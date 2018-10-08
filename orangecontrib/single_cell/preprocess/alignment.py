import numpy as np
import itertools as it
from fastdtw import fastdtw
from scipy.stats import spearmanr, scoreatpercentile
from orangecontrib.single_cell.preprocess.biweight import biweight_midcorrelation
from orangecontrib.single_cell.preprocess.cca import MultiCCA, SVDCCA


GENE_SCORING_METHODS = ("pearson", "spearman", "bicor")


def score_genes(Xs, Ws, n_metagenes=30, method=GENE_SCORING_METHODS[0]):
    """ Score genes To construct a metagene.
        Each CC component gets its own set of characteristic genes. """
    n_datasets = len(Xs)
    n_genes = Xs[0].shape[1]
    n_cc = Ws[0].shape[1]
    use_genes = dict()

    # Each CC is scored independently
    for j in range(n_cc):
        gene_scores = np.zeros((n_genes, n_datasets))
        if method == "pearson":
            for d in range(n_datasets):
                Nx = np.linalg.norm(Xs[d], axis=0)
                gene_scores[:, d] = Ws[d][:, j].T.dot(Xs[d]) / Nx
        elif method == "spearman":
            for d, g in it.product(range(n_datasets), range(n_genes)):
                gene_scores[g, d] = spearmanr(Xs[d][:, g], Ws[d][:, j])[0]
        elif method == "bicor":
            for d, g in it.product(range(n_datasets), range(n_genes)):
                gene_scores[g, d] = biweight_midcorrelation(Xs[d][:, g], Ws[d][:, j])
        else:
            raise ValueError("Parameter method must be in {pearson, spearman, bicor}.")

        # Apply sorting criteria; high is good, low is bad
        f0 = np.isnan(gene_scores).sum(axis=1)
        f1 = - np.product(np.sign(gene_scores[:, 0:1]) == np.sign(gene_scores[:, 1:]), axis=1)
        f2 = - np.min(np.absolute(gene_scores), axis=1)
        order = sorted(range(n_genes), key=lambda i: (f0[i], f1[i], f2[i]))
        use_genes[j] = order[:n_metagenes]

    return use_genes


def reference_range(x, min_p=2.5, max_p=97.5):
    """ Map values in x to its (min_p), (max_p) quantile range. """
    qmin = scoreatpercentile(x, min_p)
    qmax = scoreatpercentile(x, max_p)
    return (x - qmin) / (qmax - qmin)


def quantile_shift(x, y):
    """ Correct y for quantile shifts with reference to x. """
    shift = np.min([abs(scoreatpercentile(x, z) - scoreatpercentile(y, z))
                    for z in range(10, 100, 10)])
    return y + shift


def metagene_map(Xs, Ws, use_genes, align=True, qmin=2.5, qmax=97.5):
    """ Linear map to the meta genes.
        Assumption: Xs (and Ws) are sorted on descending number of rows. """
    n_datasets = len(Xs)
    n_cc = Ws[0].shape[1]
    Phis = [np.zeros((X.shape[0], n_cc)) for X in Xs]

    # Project to the selected set of metagenes
    for j, genes in use_genes.items():
        for X, W, Phi in zip(Xs, Ws, Phis):
            Phi[:, j] = W[:, j].T.dot(X[:, genes]).dot(X[:, genes].T).T

        # Align all datasets to the largest one
        if align:
            Phis[0][:, j] = reference_range(Phis[0][:, j], min_p=qmin, max_p=qmax)
            for d in range(1, n_datasets):
                Phis[d][:, j] = quantile_shift(Phis[0][:, j],
                                               reference_range(Phis[d][:, j]))

    return Phis


def shared_correlation(Xs, Ws, use_genes):
    """ Measure shared correlation strength for each CC component in each dataset. """
    n_datasets = len(Xs)
    n_cc = Ws[0].shape[1]
    C = np.zeros((n_datasets, n_cc))
    for ci in range(n_cc):
        for di, (X, W) in enumerate(zip(Xs, Ws)):
            genes = use_genes[ci]
            C[di, ci] = np.nanmean([abs(biweight_midcorrelation(X[:, g], W[:, ci])) for g in genes])
    return C


def duplicated(x):
    """ Mirror R function. """
    z = np.zeros((len(x),), dtype=bool)
    count = dict()
    for i in range(len(x)):
        z[i] = x[i] in count
        count[x[i]] = True
    return z


def align(Ws):
    """ Align data sets with dynamic time warping.
        Largest of the datasets (first in Ws) is chosen as the reference. """
    n_cc = Ws[0].shape[1]
    n_datasets = len(Ws)

    # Process each component independently
    for d in range(1, n_datasets):
        for j in range(n_cc):
            # Alignment makes sense only if the two vectors are sorted,
            # as there is no other logical correspondence between the cells in the two datasets.
            wl = Ws[0][:, j]
            ws = Ws[d][:, j]
            assert len(wl) >= len(ws)
            ars = np.argsort(ws)
            arl = np.argsort(wl)
            _, path = fastdtw(wl[arl], ws[ars])

            # Map smaller data set to larger
            P = np.array(path)
            keep = np.logical_not(duplicated(P[:, 1]))
            Pd = P[keep]

            # New values must be read from the larger vector (wl),
            # otherwise indices do not make sense
            Ws[d][ars[Pd[:, 1]], j] = Ws[0][arl[Pd[:, 0]], j]

    return Ws


class SeuratAlignmentModel:
    """
    Implementation of the Seurat multi data set alignment method,
    based on Canonical Correlation Analysis (CCA)
    and Dynamic Time Warping (DTW).
    """

    def __init__(self, n_components=20, n_metagenes=30, gene_scoring=GENE_SCORING_METHODS[0], random_state=None):
        """
        :param n_components: Number of components of final projection.
        :param n_metagenes: Number of characteristic genes for each CCA vector.
        :param gene_scoring: Correlation method to select metagenes. One of {'pearson', 'spearman'}.
        :param random_state: Random seed.
        """
        if gene_scoring not in GENE_SCORING_METHODS:
            raise ValueError("Parameter gene_scoring must be in %s." % str(GENE_SCORING_METHODS))
        self.gene_scoring = gene_scoring
        self.n_components = n_components
        self.n_metagenes = n_metagenes
        self.random_state = random_state
        self.Ws = None
        self.use_genes = None
        self.correlations = None
        self.shared_correlations = None

    def fit_transform(self, X, y):
        """
        Fit model parameters and return transformed data.
        :param X: Numeric matrix.
        :param y: Class (data set) labels.
        :return:
        """
        # Datasets are sorted by size for later alignment
        ys = sorted(set(y), key=lambda yi: sum(y == yi), reverse=True)
        index_map = dict([(yi, np.where(y == yi)[0]) for yi in ys])
        Xs = [X[index_map[yi]] for yi in ys]

        # Step 1: Run canonical correlation analysis
        if len(Xs) == 2:
            cca_model = SVDCCA(n_components=self.n_components,
                               random_state=self.random_state,
                               standardize=True)
            Ws = cca_model.fit_transform(*Xs)
            self.correlations = cca_model.correlations
        elif len(Xs) > 2:
            cca_model = MultiCCA(n_components=self.n_components,
                                 random_state=self.random_state,
                                 standardize=True)
            Ws = cca_model.fit_transform(Xs)
            self.correlations = cca_model.correlations
        else:
            raise ValueError("At least two data sets (indicated by y) are required!")

        # Step 2: Select metagenes and map to reference range
        self.use_genes = score_genes(Xs, Ws, n_metagenes=self.n_metagenes, method=self.gene_scoring)
        Phis = metagene_map(Xs, Ws, self.use_genes)

        # Step 2a: Compute shared correlations
        self.shared_correlations = shared_correlation(Xs, Ws, use_genes=self.use_genes)

        # Step 3: Run dynamic time warping
        Us = align(Phis)

        # Aggregate results according to y
        Z = np.zeros((X.shape[0], self.n_components))
        for U, yi in zip(Us, ys):
            Z[index_map[yi]] = U

        return Z

    def fit(self, X, y):
        """
        Fit model parameters.
        :param X: Numeric matrix.
        :param y: Class (data set) labels.
        :return: Ws: Model.
        """
        # Datasets are sorted by size for later alignment
        ys = sorted(set(y), key=lambda yi: sum(y == yi), reverse=True)
        index_map = dict([(yi, np.where(y == yi)[0]) for yi in ys])
        Xs = [X[index_map[yi]] for yi in ys]

        # Step 1: Run canonical correlation analysis
        if len(Xs) == 2:
            cca_model = SVDCCA(n_components=self.n_components,
                               random_state=self.random_state,
                               standardize=True)
            Ws = cca_model.fit_transform(*Xs)
            self.correlations = cca_model.correlations
        elif len(Xs) > 2:
            cca_model = MultiCCA(n_components=self.n_components,
                                 random_state=self.random_state,
                                 standardize=True)
            Ws = cca_model.fit_transform(Xs)
            self.correlations = cca_model.correlations
        else:
            raise ValueError("At least two data sets (indicated by y) are required!")

        # Step 2: Select metagenes and map to reference range
        self.use_genes = score_genes(Xs, Ws, n_metagenes=self.n_metagenes, method=self.gene_scoring)

        # Step 2a: Compute shared correlations
        self.shared_correlations = shared_correlation(Xs, Ws, use_genes=self.use_genes)

        self.Ws = Ws
        return Ws

    def transform(self, X, y, normalize=True, quantile=2.5, dtw=True, Ws=None):
        """
        Transform given data based on a model.
        :param X: Numeric matrix.
        :param y: Class (data set) labels.
        :param normalize: Boolean, whether to apply quantile normalization.
        :param quantile: float (between 0 and 49(, on which quantile will be transformed data normalized.
        :param dtw: Boolean, wheter to apply Dynamic time warping.
        :param Ws: Pretrained model.
        :return: Transformed data.
        """
        if Ws is None:
            Ws = self.Ws

        # Datasets are sorted by size for later alignment
        ys = sorted(set(y), key=lambda yi: sum(y == yi), reverse=True)
        index_map = dict([(yi, np.where(y == yi)[0]) for yi in ys])
        Xs = [X[index_map[yi]] for yi in ys]

        # Step 2: Map to reference range
        qmin = quantile
        qmax = 100 - quantile
        Us = metagene_map(Xs, Ws, self.use_genes, normalize, qmin, qmax)

        # Step 3: Run dynamic time warping
        if dtw:
            Us = align(Us)

        # Aggregate results according to y
        Z = np.zeros((X.shape[0], self.n_components))
        for U, yi in zip(Us, ys):
            Z[index_map[yi]] = U

        return Z

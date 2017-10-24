import numpy as np
import scipy.stats as spstats
from Orange.regression import Learner, Model, SklLearner
from Orange.data import Table
from statsmodels.regression.quantile_regression import QuantReg
from sklearn.preprocessing import PolynomialFeatures


class QuantRegLearner(Learner):
    """Generate a quantile regression modeland learn a prediction model

    Parameters
    ----------
    quantile : float in [0, 1]
        target quantile
    preprocessors : List[Preprocessor]
        preprocessors to be applied on the data before learning
    """
    name = 'quantile regression'
    preprocessors = SklLearner.preprocessors

    def __init__(self, quantile=0.5, preprocessors=None):
        super().__init__(preprocessors=preprocessors)
        self.quantile = quantile

    def fit(self, X, Y, W=None):
        qr = QuantReg(Y, X)
        model = qr.fit(self.quantile)
        params = model.params
        return QuantRegModel(model, params)


class QuantRegModel(Model):
    def __init__(self, model, params):
        self.model = model
        self.params = params
        super().__init__()

    def predict(self, X):
        return self.model.predict(X)

    def __str__(self):
        return 'QuantRegModel {}'.format(self.model)


class ScNormModel:
    """
    Single cell RNA-seq normalization based on Quantile normalization.
    Genes are split in K groups based on quantile linear regression slopes.
    The data is corrected such that the expression does not depend on sequencing depth.

    Bacher, Rhonda, et al. "SCnorm: robust normalization of single-cell RNA-seq data."
    Nature Methods 14.6 (2017): 584-586.
    """

    def __init__(self, p_subgroup=0.25, K=10, n_genes=None):
        """
        :param p_subgroup: Proportion of genes within a subgroup.
        :param n_genes: Number of genes to use for fitting of slopes.
        :param K: number of groups.
        """
        self.p_subgroup = p_subgroup
        self.n_genes = n_genes
        self.K = K

        # Fixed hyperparameters
        self.q_range = np.linspace(0.05, 0.95, 19)

        # Construct a dict of group-specific models,
        # characterized by gene median expression
        self.group_models = dict()


    def fit(self, X, Y=None):
        """
        Determine groups of genes based on quantile regression slopes.
        :param X: Cell-gene expression matrix with raw counts.
        :param Y: Unused.
        :return:
        """

        # Remove genes expressed in less than two cells
        gen_sum = np.array((X > 0), dtype=int).sum(axis=0).ravel()
        X = X[:, gen_sum > 1]

        # Calculate sequencing depth for each cell in the training set
        n, m = X.shape
        cell_sum = X.sum(axis=1).reshape((n, 1))

        # Filter genes for speed (use groups based on medians)
        X[X == 0] = np.nan
        med_exprs = np.nanmedian(X, axis=0)

        # Partitioning of genes into equally-sized groups based on percentile of median expression
        percentiles = np.array(list(map(lambda e: spstats.percentileofscore(med_exprs, e), med_exprs)))
        groups = np.array(self.K * percentiles / 100.0, dtype=int)
        groups[groups == self.K] = self.K - 1

        # Process each group independently
        # If a subset of genes is specified, use an equal number of genes closest to
        # the group median
        for k in range(self.K):
            group_cols = np.where(groups == k)[0]
            group_med_expr = med_exprs[group_cols]

            # Sub-sample group if needed ; n_group genes closest to group median are taken
            group_protyps = group_cols
            if self.n_genes is not None:
                n_group = max(int(self.n_genes / self.K), 3)
                dist = np.absolute(group_med_expr - np.median(group_med_expr))
                group_protyps = group_cols[np.argsort(dist)[:n_group]]

            # Compute slopes for selected genes
            group_slopes = np.zeros((len(group_protyps),))
            for ji, j in enumerate(group_protyps):
                Y = X[:, j]
                rows = np.where(np.isnan(Y) == False)[0]
                assert len(rows) >= 2
                Dp = PolynomialFeatures(degree=1).fit_transform(np.log(cell_sum[rows]))
                depth_model = QuantRegLearner().fit(Dp, np.log(Y[rows]), W=None)
                group_slopes[ji] = depth_model.params[1]

            # Convert to long format
            X_group = X[:, group_protyps]
            use_rows, use_cols = np.where(np.logical_not(np.isnan(X_group)))
            D_sub = cell_sum[use_rows].reshape((len(use_rows), 1))
            Y_sub = X_group[(use_rows, use_cols)].reshape((len(use_rows), 1))

            # Compute a model for the whole group and choose one
            # with best matching to the group median slope
            group_med_slope = np.median(group_slopes)
            submodels = dict()
            for q in self.q_range:
                Dp = PolynomialFeatures(degree=1).fit_transform(np.log(D_sub))
                model = QuantRegLearner(quantile=q).fit(Dp, np.log(Y_sub), W=None)
                model_dist = np.absolute(group_med_slope - model.params[1])
                submodels[q] = model_dist, model

            q_best, (_, model_best) = sorted(submodels.items(), key=lambda tup: tup[1])[0]
            group_quantile = spstats.scoreatpercentile(np.log(Y_sub), 100.0 * q_best)
            group_median = np.median(group_med_expr)
            self.group_models[k] = (model_best, group_median, group_quantile)


    def group_genes(self, data):
        """
        Map genes to groups.
        :param data:
        :return:
        """
        # Temporary use nan to compute non-zero median without copying the data
        X_new = data.X
        X_new[np.where(X_new == 0)] = np.nan

        group_keys = sorted(self.group_models.keys())
        group_expr_med = np.array(list(map(lambda ky: self.group_models[ky][1], group_keys)))

        gen_expr_med = np.nanmedian(X_new, axis=0)
        gen_groups = np.array(list(map(lambda dm:
                                        group_keys[np.argmin((group_expr_med - dm) ** 2)],
                                        gen_expr_med)))

        X_new[np.where(np.isnan(X_new))] = 0
        return gen_groups


    def transform(self, data):
        """
        Map new data to groups based on median expressions.
        Correct values using previously fit models.
        :param data: Orange.data.Table with cells in rows and genes in columns.
        :return: Orange.data.Table with normalized values.
        """
        X_new = data.X.copy()
        cell_sum = data.X.sum(axis=1)
        group_keys = sorted(self.group_models.keys())
        gen_groups = self.group_genes(data)

        for gky in group_keys:
            # Extract values for group
            group_cols = np.where(gen_groups == gky)[0]
            use_rows, use_cols = inxs = np.nonzero(data.X[:, group_cols])
            D_sub = cell_sum[use_rows].reshape((len(use_rows), 1))

            # Predict new values
            model, _, quantile = self.group_models[gky]
            Dp = PolynomialFeatures(degree=1).fit_transform(np.log(D_sub))
            predicted = model.predict(Dp)
            size_factors = np.exp(predicted) / np.exp(quantile)

            # Assign to X_new (cannot assign values to double indexing X[...][...])
            X_tmp = X_new[:, group_cols].copy()
            X_tmp[inxs] = X_tmp[inxs] / size_factors
            X_new[:, group_cols] = X_tmp

        return Table.from_numpy(domain=data.domain,
                                X=X_new,
                                Y=data.Y,
                                metas=data.metas,
                                W=data.W)




def test():
    n_genes = 1000
    test_file = "/Users/martin/Dev/data/singlecell/bacher2017/H1_data.csv"
    data = Table(test_file)[:, :n_genes]

    model = ScNormModel(K=5)
    model.fit(data.X)
    normed_data = model.transform(data)

    xs = np.log(data.X.ravel() + 1)
    ys = np.log(normed_data.X.ravel() + 1)

    print("Difference: %f" % np.linalg.norm(data.X - normed_data.X))
    print("Sp. correlation: %f (%f)" % spstats.spearmanr(xs, ys))

if __name__ == "__main__":
    test()
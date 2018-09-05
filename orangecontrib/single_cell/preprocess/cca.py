import itertools as it
import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition.truncated_svd import TruncatedSVD

import Orange.statistics.util as ut


def _standardize(X, tol=1e-8):
    """ Standardize matrix X avoiding NaNs and small values. """
    s = ut.std(X, axis=0)
    m = np.zeros(*s.shape)
    valid = s > tol
    m[np.logical_not(valid)] = 0
    m[valid] = 1.0 / s[valid]
    return (X - X.mean(axis=0)) * m


class SVDCCA:
    """ Two data set CCA approximated via SVD. Follows
        Butler, A. et al. (2018) "Integrating single-cell transcriptomic data across different conditions,
        technologies, and species",  Nature Biotechnology, (July 2017). doi: 10.1038/nbt.4096.
    """

    def __init__(self, n_components, random_state=None, standardize=True):
        self.n_components = n_components
        self.standardize = standardize
        self.random_state = random_state
        self.correlations = None

    def fit_transform(self, X, Y):
        if self.standardize:
            X = _standardize(X)
            Y = _standardize(Y)
        K = X.dot(Y.T)
        model = TruncatedSVD(n_components=self.n_components,
                             random_state=self.random_state)
        U = model.fit_transform(K)
        U = U / np.linalg.norm(U, axis=0)
        V = model.components_.T
        self.correlations = np.array([pearsonr(u.dot(X), v.dot(Y))[0]
                                      for u, v in zip(U.T, V.T)])
        return U, V


class MultiCCA:

    """ Approximation to penalized CCA for multiple datasets.
        Follows Seurat MultiCCA method.
    """

    def __init__(self, n_components, max_iter=25, tol=1e-3, random_state=None, standardize=True):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        self.standardize = standardize
        self.correlations = None

    @staticmethod
    def _objective(Xs, ws):
        """ CCA objective function. """
        c = 0
        for i, j in it.combinations(range(len(Xs)), 2):
            c += ws[i].T.dot(Xs[i]).dot(Xs[j].T.dot(ws[j]))
        return c

    def fit_transform(self, Xs):
        """
        Optimize each CC components and return per-data set projections.
        :param Xs: List of matrices with the same number of columns.
        :return: CCA subspace.
        """
        p = len(Xs)
        Ws = [np.zeros((X.shape[0], self.n_components)) for X in Xs]
        if self.standardize:
            Xs = list(map(_standardize, Xs))
        Ws_init = [TruncatedSVD(n_components=self.n_components,
                                random_state=self.random_state).fit_transform(X) for X in Xs]
        correlations = np.zeros((self.n_components,))


        # Optimize each CC component individually
        for cc in range(self.n_components):
            w_cur = [Wi[:, cc] / np.linalg.norm(Wi[:, cc]) for Wi in Ws_init]
            for itr in range(self.max_iter):
                o1 = self._objective(Xs, w_cur)
                for i in range(p):
                    wi = 0
                    for j in range(p):
                        if i == j:
                            continue
                        wj = w_cur[j]
                        Dj = np.diag(np.diagonal(Ws[i].T.dot(Xs[i]).dot(Xs[j].T.dot(Ws[j]))))
                        wi += Xs[i].dot((Xs[j].T.dot(wj))) - Ws[i].dot(Dj).dot(Ws[j].T).dot(wj)
                    w_cur[i] = wi / np.linalg.norm(wi)
                o2 = self._objective(Xs, w_cur)
                if abs(o2-o1) / abs(o1) < self.tol:
                    break
            for i in range(p):
                Ws[i][:, cc] = w_cur[i]

            # Compute average correlations
            n_pairs = p * (p - 1) / 2
            for i, j in it.combinations(range(p), 2):
                wi = Ws[i][:, cc].T.dot(Xs[i])
                wj = Ws[j][:, cc].T.dot(Xs[j])
                correlations[cc] += pearsonr(wi, wj)[0] / n_pairs

        # Orientate vectors
        s = np.sign(Ws[0][0, :])
        for i in range(p):
            Ws[i] = Ws[i] * s

        self.correlations = correlations
        return Ws

import unittest
import numpy as np
from orangecontrib.single_cell.preprocess.cca import MultiCCA, SVDCCA
from scipy.stats import pearsonr


class TestSVDCCA(unittest.TestCase):

    def test_random_data(self):
        """ Assert approximately decreasing explained variance of components."""
        np.random.seed(42)
        n1 = 200
        n2 = 300
        m = 100
        k = 20
        X = np.random.rand(n1, m)
        Y = np.random.rand(n2, m)
        model = SVDCCA(n_components=k, random_state=42)
        U, V = model.fit_transform(X, Y)
        slope, pvalue = pearsonr(model.correlations, np.arange(k))
        assert slope < 0
        assert pvalue < 0.05
        assert np.linalg.norm(U.T.dot(U) - np.eye(k)) < 1e-5
        assert np.linalg.norm(V.T.dot(V) - np.eye(k)) < 1e-5


class TestMultiCCA(unittest.TestCase):

    def test_random_data(self):
        """ Test approximately decreasing explained variance of components."""
        np.random.seed(42)
        ns = [10, 12, 14]
        m = 100
        k = 10
        Xs = [np.random.rand(n, m) for n in ns]
        model = MultiCCA(n_components=k, random_state=42)
        Ws = model.fit_transform(Xs)
        slope, pvalue = pearsonr(model.correlations, np.arange(k))
        assert slope < 0
        assert pvalue < 0.05
        for W in Ws:
            # Normal but not orthogonal
            assert np.linalg.norm(np.linalg.norm(W, axis=0) - np.ones((k,))) < 1e-5

if __name__ == "__main__":
    unittest.main()

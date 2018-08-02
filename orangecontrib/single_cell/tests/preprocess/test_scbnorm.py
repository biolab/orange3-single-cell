import numpy as np
import unittest
from itertools import product
from scipy.stats import pearsonr
from Orange.data import Table, ContinuousVariable, DiscreteVariable, Domain
from orangecontrib.single_cell.preprocess.scbnorm import *


class ScBatchNormalize(unittest.TestCase):

    def setUp(self):
        """ Set up datasets. Insert 0.5% random zeros. """
        np.random.seed(42)
        n = 100
        m = 200
        bias = 10
        self.noise = 0.1 * np.random.randn(n, 1)
        Z = np.hstack([np.random.randn(n, 1),
                       np.random.choice([0, 1], size=n).reshape(n, 1)])
        Y = np.random.randn(n, 1)
        W = np.random.randn(2, m)
        X_log = INV_LINKS[LINK_LOG](bias + Z.dot(W) + self.noise)
        X_lin = INV_LINKS[LINK_IDENTITY](bias + Z.dot(W) + self.noise)
        X_lin_neg = INV_LINKS[LINK_IDENTITY](Z.dot(W))

        self.zeros = (np.random.choice(range(n), size=100, replace=True),
                      np.random.choice(range(m), size=100, replace=True))
        X_log[self.zeros] = 0
        X_lin[self.zeros] = 0

        self.data_lin = Table.from_numpy(X=X_lin,
                                         metas=Z,
                                         Y=Y,
                                         domain=Domain(
                                            class_vars=[DiscreteVariable("class", values=["c1", "c2"])],
                                            attributes=[ContinuousVariable("X%d" % i) for i in range(m)],
                                            metas=[ContinuousVariable("Z0"),
                                                   DiscreteVariable("Z1", values=["a", "b"])]))

        self.data_log = Table.from_numpy(X=X_log,
                                         metas=Z,
                                         Y=Y,
                                         domain=Domain(
                                             class_vars=[DiscreteVariable("class", values=["c1", "c2"])],
                                             attributes=[ContinuousVariable("X%d" % i) for i in range(m)],
                                             metas=[ContinuousVariable("Z0"),
                                                    DiscreteVariable("Z1", values=["a", "b"])]))

        self.data_lin_neg = Table.from_numpy(X=X_lin_neg,
                                             metas=Z,
                                             Y=Y,
                                             domain=Domain(
                                                 class_vars=[DiscreteVariable("class", values=["c1", "c2"])],
                                                 attributes=[ContinuousVariable("X%d" % i) for i in range(m)],
                                                 metas=[ContinuousVariable("Z0"),
                                                        DiscreteVariable("Z1", values=["a", "b"])]))

    def test_multivariate_lin(self):
        """ Assert batch variables are decorrelated and the remaining factors (noise) become correlated. """
        alpha = 0.05
        model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"], link=LINK_IDENTITY)
        model.fit(self.data_lin)
        data1 = model(self.data_lin)
        assert pearsonr(data1.get_column_view("Z0")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(data1.get_column_view("Z1")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(self.noise.ravel(), data1.X[:, 0].ravel())[1] < alpha
        assert np.linalg.norm(data1.X[self.zeros]) == 0

    def test_multivariate_log(self):
        """ Assert batch variables are decorrelated and the remaining factors (noise) become correlated. """
        alpha = 0.05
        model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"], link=LINK_LOG)
        model.fit(self.data_log)
        data1 = model(self.data_log)
        assert pearsonr(data1.get_column_view("Z0")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(data1.get_column_view("Z1")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(self.noise.ravel(), data1.X[:, 0].ravel())[1] < alpha
        assert np.linalg.norm(data1.X[self.zeros]) == 0

    def test_preprocessor_lin(self):
        """ Test as a preprocessor (linear link). """
        pp = SCBatchNormalizer(batch_vars=["Z0", "Z1"], link=LINK_IDENTITY)
        data1 = pp(self.data_lin)
        data2 = self.data_lin.transform(data1.domain)
        np.testing.assert_almost_equal(data1.X, data2.X)
        assert np.linalg.norm(data1.X[self.zeros]) == 0

    def test_preprocessor_log(self):
        """ Test as a preprocessor (log link). """
        pp = SCBatchNormalizer(batch_vars=["Z0", "Z1"], link=LINK_LOG)
        data1 = pp(self.data_log)
        data2 = self.data_log.transform(data1.domain)
        np.testing.assert_almost_equal(data1.X, data2.X)
        assert np.linalg.norm(data1.X[self.zeros]) == 0

    def test_preprocessor_vs_model_log(self):
        """ Normalizer and model should return the same values (log). """
        pp = SCBatchNormalizer(batch_vars=["Z0", "Z1"], link=LINK_LOG)
        model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"], link=LINK_LOG)
        model.fit(self.data_log)
        data1 = model(self.data_log)
        data2 = pp(self.data_log)
        np.testing.assert_almost_equal(data1.X, data2.X)

    def test_preprocessor_vs_model_lin(self):
        """ Normalizer and model should return the same values (lin). """
        pp = SCBatchNormalizer(batch_vars=["Z0", "Z1"], link=LINK_IDENTITY, nonzero_only=False)
        model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"], link=LINK_IDENTITY, nonzero_only=False)
        model.fit(self.data_lin)
        data1 = model(self.data_lin)
        data2 = pp(self.data_lin)
        np.testing.assert_almost_equal(data1.X, data2.X)

    def test_nonzero_lin(self):
        """ Fit all including zeros (linear link). """
        alpha = 0.05
        model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"],
                                      nonzero_only=False,
                                      link=LINK_IDENTITY)
        model.fit(self.data_lin)
        data1 = model(self.data_lin)
        assert pearsonr(data1.get_column_view("Z0")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(data1.get_column_view("Z1")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(self.noise.ravel(), data1.X[:, 0].ravel())[1] < alpha
        assert np.linalg.norm(data1.X[self.zeros]) > 0

    def test_nonzero_log(self):
        """ Fit all including zeros (log link)"""
        def dummy_call():
            return ScBatchNormalizeModel(batch_vars=["Z0", "Z1"],
                                         nonzero_only=False,
                                         link=LINK_LOG)
        self.assertRaises(ValueError, dummy_call)

    def test_batch_class(self):
        """ Test including class variable. """
        alpha = 0.05
        pp = SCBatchNormalizer(batch_vars=("class", "Z0"))
        data1 = pp(self.data_lin)
        assert pearsonr(data1.get_column_view("class")[0], data1.X[:, 0].ravel())[1] > alpha
        assert pearsonr(data1.get_column_view("Z0")[0], data1.X[:, 0].ravel())[1] > alpha

    def test_no_vars(self):
        """ Test with no batch variables. """
        pp = SCBatchNormalizer(batch_vars=())
        data1 = pp(self.data_lin)
        assert np.linalg.norm(data1.X - self.data_lin.X) == 0

    def test_scorer_correlations(self):
        """ Matrix form correlations should equal default numpy implementations. """
        A = self.data_lin.X[:, :10]
        B = self.data_log.X[:, :10]
        C, P = ScBatchScorer.correlations(A, B)
        for i, j in product(range(A.shape[1]), range(B.shape[1])):
            c, p = pearsonr(A[:, i], B[:, j])
            assert abs(C[i, j] - c) < 1e-3
            assert abs(P[i, j] - p) < 1e-3

    def test_negative_lin(self):
        """ Assert batch variables are decorrelated and the remaining factors (noise) become correlated. """
        def dummy_call():
            model = ScBatchNormalizeModel(batch_vars=["Z0", "Z1"],
                                          nonzero_only=True,
                                          link=LINK_IDENTITY)
            model.fit(self.data_lin_neg)
        self.assertRaises(ValueError, dummy_call)
import unittest

import numpy as np
import numpy.testing as npt

from Orange.data import Table, Domain
from orangecontrib.single_cell.preprocess.scpreprocess import (
    LogarithmicScale, Binarize, NormalizeSamples, NormalizeGroups,
    Standardize, SelectMostVariableGenes
)


class TestLogarithmicScale(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(1, 13).reshape((3, 4)))

    def test_default(self):
        table = LogarithmicScale()(self.table)
        self.assertIsInstance(table, Table)
        npt.assert_array_equal(table, np.log2(1 + self.table.X))
        self.assertNotEqual(table, self.table)

    def test_options(self):
        table = LogarithmicScale(LogarithmicScale.NaturalLog)(self.table)
        npt.assert_array_equal(table, np.log(1 + self.table.X))


class TestBinarize(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(12).reshape((3, 4)))

    def test_default(self):
        table = Binarize()(self.table)
        self.assertIsInstance(table, Table)
        self.assertEqual(np.sum(table.X.ravel()[:1]), 0)
        self.assertEqual(np.sum(table.X.ravel()[1:]), 11)
        self.assertNotEqual(table, self.table)

    def test_options(self):
        table = Binarize(Binarize.Greater, 5)(self.table)
        self.assertEqual(np.sum(table.X.ravel()[:6]), 0)
        self.assertEqual(np.sum(table.X.ravel()[6:]), 6)


class TestNormalizeSamples(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(1, 7).reshape((2, 3)).T)

    def test_default(self):
        table = NormalizeSamples()(self.table)
        self.assertIsInstance(table, Table)
        self.assertNotEqual(table, self.table)
        npt.assert_array_almost_equal(
            table, [[200000, 800000],
                    [285714.286, 714285.714],
                    [333333.333, 666666.667]], decimal=3)

    def test_options(self):
        table = NormalizeSamples(NormalizeSamples.Median)(self.table)
        npt.assert_array_almost_equal(
            table, [[1.4, 5.6], [2, 5], [2.333, 4.667]], decimal=3)


class TestNormalizeGroups(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        table = Table([[1, 4, 0], [2, 5, 0], [3, 6, 1]])
        domain = Domain(table.domain.attributes[:-1],
                        metas=table.domain.attributes[-1:])
        cls.table = table.transform(domain)

    def test_default(self):
        table = NormalizeGroups(-1)(self.table)
        self.assertIsInstance(table, Table)
        self.assertNotEqual(table, self.table)
        npt.assert_array_almost_equal(
            table.X, [[83333.333, 333333.333],
                      [166666.667, 416666.667],
                      [333333.333, 666666.667]], decimal=3)

    def test_options(self):
        table = NormalizeGroups(-1, NormalizeGroups.Median)(self.table)
        npt.assert_array_equal(table.X, [[0.5, 2], [1, 2.5], [2, 4]])


class TestStandardize(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(1, 7).reshape((2, 3)).T)

    def test_default(self):
        table = Standardize()(self.table)
        self.assertIsInstance(table, Table)
        self.assertNotEqual(table, self.table)
        npt.assert_array_almost_equal(
            table, [[-1.225, -1.225], [0., 0.], [1.225, 1.225]], decimal=3)

    def test_options(self):
        table = Standardize(-1, 1)(self.table)
        npt.assert_array_equal(table, [[-1, -1], [0., 0.], [1, 1]])


class TestSelectMostVariableGenes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = Table("dermatology")
        cls.data = data.transform(Domain(data.domain.attributes[:-1],
                                         data.domain.class_vars))
        cls.variance = np.nanvar(cls.data.X, axis=0)
        cls.mean = np.nanmean(cls.data.X, axis=0)

    def test_default(self):
        table = Table(np.arange(12).reshape((3, 4)).T ** (1 / 2))
        pp_table = SelectMostVariableGenes(n_genes=2, n_groups=2)(table)
        self.assertIsInstance(pp_table, Table)
        self.assertNotEqual(pp_table, table)
        npt.assert_array_equal(pp_table, table[:, :2])

    def test_options_dispersion(self):
        attrs, cls = self.data.domain.attributes, self.data.domain.class_vars
        inds = sorted(np.argsort(self.variance / self.mean))
        data = self.data.transform(Domain(tuple(np.array(attrs)[inds]), cls))
        npt.assert_array_equal(
            SelectMostVariableGenes(method=SelectMostVariableGenes.Dispersion,
                                    n_groups=None)(self.data), data)

    def test_options_variance(self):
        attrs, cls = self.data.domain.attributes, self.data.domain.class_vars
        inds = sorted(np.argsort(self.variance))
        data = self.data.transform(Domain(tuple(np.array(attrs)[inds]), cls))
        npt.assert_array_equal(
            SelectMostVariableGenes(method=SelectMostVariableGenes.Variance,
                                    n_groups=None)(self.data), data)

    def test_options_mean(self):
        attrs, cls = self.data.domain.attributes, self.data.domain.class_vars
        inds = sorted(np.argsort(self.mean))
        data = self.data.transform(Domain(tuple(np.array(attrs)[inds]), cls))
        npt.assert_array_equal(
            SelectMostVariableGenes(method=SelectMostVariableGenes.Mean,
                                    n_groups=None)(self.data), data)

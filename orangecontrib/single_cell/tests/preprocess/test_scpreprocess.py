import warnings
import os
import unittest
import numpy as np
import numpy.testing as npt

from Orange.data import Table, Domain
from Orange.widgets.tests.utils import table_dense_sparse
from orangecontrib.single_cell.preprocess.scpreprocess import (
    LogarithmicScale, Binarize, NormalizeSamples, NormalizeGroups,
    Standardize, SelectMostVariableGenes, DropoutGeneSelection,
    DropoutWarning
)


class TestLogarithmicScale(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(1, 13).reshape((3, 4)))

    @table_dense_sparse
    def test_default(self, prepare_table):
        table = prepare_table(self.table)
        new_table = LogarithmicScale()(table)
        self.assertIsInstance(new_table, Table)
        npt.assert_array_equal(new_table, np.log2(1 + self.table.X))
        self.assertNotEqual(new_table, self.table)

    @table_dense_sparse
    def test_options(self, prepare_table):
        table = prepare_table(self.table)
        new_table = LogarithmicScale(LogarithmicScale.NaturalLog)(table)
        npt.assert_allclose(new_table, np.log(1 + self.table.X))

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        new_table = LogarithmicScale(LogarithmicScale.NaturalLog)(table)
        self.assertEqual(is_sparse, new_table.is_sparse())


class TestBinarize(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.arange(12).reshape((3, 4)))

    @table_dense_sparse
    def test_default(self, prepare_table):
        table = prepare_table(self.table)
        new_table = Binarize()(table)
        self.assertIsInstance(new_table, Table)
        self.assertEqual(new_table.X[0, 0], 0)
        self.assertEqual(np.sum(new_table[1]), 4)
        self.assertNotEqual(new_table, self.table)

    @table_dense_sparse
    def test_options(self, prepare_table):
        table = prepare_table(self.table)
        new_table = Binarize(Binarize.Greater, 5)(table)
        self.assertEqual(np.sum(new_table.X[0]), 0)
        self.assertEqual(np.sum(new_table.X[1]), 2)
        self.assertEqual(np.sum(new_table.X[2]), 4)

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        new_table = Binarize(Binarize.Greater, 5)(table)
        self.assertEqual(is_sparse, new_table.is_sparse())


class TestNormalizeSamples(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = Table(np.array([
            [1.5, 8.5, 0],
            [2, 0, 3],
            [3, np.nan, 7],
        ]))

    @table_dense_sparse
    def test_default(self, prepare_table):
        table = prepare_table(self.table)
        new_table = NormalizeSamples()(table)
        self.assertIsInstance(new_table, Table)
        self.assertNotEqual(new_table, self.table)

        new_table = new_table.to_dense()
        npt.assert_array_almost_equal(
            new_table, [[150000, 850000, 0],
                        [400000, 0, 600000],
                        [300000, np.nan, 700000]])

    @table_dense_sparse
    def test_options(self, prepare_table):
        table = prepare_table(self.table)
        new_table = NormalizeSamples(NormalizeSamples.Median)(table)
        npt.assert_array_almost_equal(
            new_table,
            [[1.5, 8.5, 0], [4, 0, 6], [3, np.nan, 7]],
        )

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        new_table = NormalizeSamples(NormalizeSamples.Median)(table)
        self.assertEqual(is_sparse, new_table.is_sparse())


class TestNormalizeGroups(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        table = Table([[1, 4, 0], [2, 5, 0], [3, 6, 1]])
        domain = Domain(table.domain.attributes[:-1],
                        metas=table.domain.attributes[-1:])
        cls.table = table.transform(domain)

    @table_dense_sparse
    def test_default(self, prepare_table):
        table = prepare_table(self.table)
        new_table = NormalizeGroups(-1)(table)
        self.assertIsInstance(new_table, Table)
        self.assertNotEqual(new_table, self.table)

        new_table = new_table.to_dense()
        npt.assert_array_almost_equal(
            new_table.X, [[83333.333, 333333.333],
                          [166666.667, 416666.667],
                          [333333.333, 666666.667]], decimal=3)

    @table_dense_sparse
    def test_options(self, prepare_table):
        table = prepare_table(self.table)
        new_table = NormalizeGroups(-1, NormalizeGroups.Median)(table)
        new_table = new_table.to_dense()
        npt.assert_array_equal(new_table.X, [[0.5, 2], [1, 2.5], [2, 4]])

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        new_table = NormalizeGroups(-1, NormalizeGroups.Median)(table)
        self.assertEqual(is_sparse, new_table.is_sparse())


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
        path = os.path.join(os.path.dirname(__file__), '..', "data")
        data = Table(os.path.join(path, 'dermatology.tab'))
        cls.table = data.transform(Domain(data.domain.attributes[:-1],
                                          data.domain.class_vars))
        cls.variance = np.nanvar(cls.table.X, axis=0)
        cls.mean = np.nanmean(cls.table.X, axis=0)

    @table_dense_sparse
    def test_default(self, prepare_table):
        table = Table(np.arange(12).reshape((3, 4)).T ** (1 / 2))
        table = prepare_table(table)
        pp_table = SelectMostVariableGenes(n_genes=2, n_groups=2)(table)
        self.assertIsInstance(pp_table, Table)
        self.assertNotEqual(pp_table, table)
        npt.assert_array_equal(pp_table, table[:, :2])

    @table_dense_sparse
    def test_options_dispersion(self, prepare_table):
        attrs, cls = self.table.domain.attributes, self.table.domain.class_vars
        inds = sorted(np.argsort(self.variance / self.mean))
        data = self.table.transform(Domain(tuple(np.array(attrs)[inds]), cls))

        table = prepare_table(data)
        filtered_data = SelectMostVariableGenes(
            method=SelectMostVariableGenes.Dispersion, n_groups=None
        )(table)
        npt.assert_array_equal(filtered_data, data)

    @table_dense_sparse
    def test_options_variance(self, prepare_table):
        attrs, cls = self.table.domain.attributes, self.table.domain.class_vars
        inds = sorted(np.argsort(self.variance))
        table = self.table.transform(Domain(tuple(np.array(attrs)[inds]), cls))
        table = prepare_table(table)

        filtered_table = SelectMostVariableGenes(
            method=SelectMostVariableGenes.Dispersion, n_groups=None
        )(table)

        npt.assert_array_equal(filtered_table, table)

    @table_dense_sparse
    def test_options_mean(self, prepare_table):
        attrs, cls = self.table.domain.attributes, self.table.domain.class_vars
        inds = sorted(np.argsort(self.mean))
        table = self.table.transform(Domain(tuple(np.array(attrs)[inds]), cls))
        table = prepare_table(table)

        filtered_table = SelectMostVariableGenes(
            method=SelectMostVariableGenes.Dispersion, n_groups=None
        )(table)

        npt.assert_array_equal(filtered_table, table)

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        filtered_table = SelectMostVariableGenes(
            method=SelectMostVariableGenes.Dispersion, n_groups=None
        )(table)
        self.assertEqual(is_sparse, filtered_table.is_sparse())


class TestDropoutGeneSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(os.path.dirname(__file__), '..', "data")
        cls.table = Table(os.path.join(path, 'dermatology.tab'))

    @table_dense_sparse
    def test_default(self, prepare_table):
        n_genes = 2
        table = prepare_table(self.table)
        filtered_table = DropoutGeneSelection(n_genes)(table)
        self.assertIsInstance(filtered_table, Table)
        self.assertEqual(len(filtered_table.domain.attributes), n_genes)

    @table_dense_sparse
    def test_detection(self, prepare_table):
        table = prepare_table(self.table)
        preprocessor = DropoutGeneSelection(n_genes=2)
        zero_rate, mean_expr = preprocessor.detection(table.X)
        self.assertEqual(len(zero_rate), len(table.domain.attributes))
        self.assertEqual(len(mean_expr), len(table.domain.attributes))
        npt.assert_allclose(
            zero_rate[:5],
            np.array([0.01092896, 0.02185792, 0.16120219, 0.32240437, 0.61202186]),
            atol=1e-7,
        )
        npt.assert_allclose(
            mean_expr[:5],
            np.array([0.9879741, 0.77491075, 0.78471751, 0.88894012, 0.58119243]),
            atol=1e-7,
        )

    @table_dense_sparse
    def test_warning(self, prepare_table):
        n_genes = 30
        table = prepare_table(self.table)
        with warnings.catch_warnings(record=True) as w:
            filtered_table = DropoutGeneSelection(n_genes)(table)
        self.assertIsInstance(filtered_table, Table)
        self.assertNotEqual(len(filtered_table.domain.attributes), n_genes)
        dropout_warnings = [warning for warning in w if
                            issubclass(warning.category, DropoutWarning)]
        self.assertEqual(len(dropout_warnings), 1)

    @table_dense_sparse
    def test_preserves_density(self, prepare_table):
        table = prepare_table(self.table)
        is_sparse = table.is_sparse()
        filtered_table = DropoutGeneSelection(3)(table)
        self.assertEqual(is_sparse, filtered_table.is_sparse())

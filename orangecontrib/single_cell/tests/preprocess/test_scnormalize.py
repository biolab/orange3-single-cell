import unittest
import numpy as np
from Orange.data.table import Table
from orangecontrib.single_cell.preprocess.scnormalize import ScNormalizeModel, ScNormalizeProjector


class ScNormalizeTest(unittest.TestCase):

    def setUp(self):
        self.table = Table("iris")

    def test_projection(self):
        # Fit with projector
        projector = ScNormalizeProjector(domain=self.table.domain)
        projection = projector.fit(self.table.X)
        proj_data = projection(self.table)
        # Fit with model
        model = ScNormalizeModel()
        model.fit(self.table.X)
        model_data = model(self.table)
        self.assertTrue(np.all(np.array(proj_data.X == model_data.X)))

    def test_new_data(self):
        inxs1 = list(range(0, len(self.table), 2))
        inxs2 = list(range(1, len(self.table), 2))
        projector = ScNormalizeProjector(domain=self.table.domain,
                                         equalize_var=None,
                                         normalize_cells=True,
                                         log_base=None)

        mi1 = np.median(self.table[inxs1].X.mean(axis=1))
        mi2 = np.median(self.table[inxs2].X.mean(axis=1))
        self.assertNotEqual(mi1, mi2)
        projection = projector.fit(self.table[inxs1].X)
        proj_data1 = projection(self.table[inxs1])
        proj_data2 = projection(self.table[inxs2])
        self.assertFalse(np.all(np.array(proj_data1.X == proj_data2.X)))
        m1 = np.median(proj_data1.X.mean(axis=1))
        m2 = np.median(proj_data2.X.mean(axis=1))
        self.assertTrue(m1 - m2 < 1)
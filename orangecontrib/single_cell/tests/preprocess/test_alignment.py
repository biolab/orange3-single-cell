import unittest
import numpy as np
from orangecontrib.single_cell.preprocess.alignment import SeuratAlignmentModel, GENE_SCORING_METHODS


class TestAlignment(unittest.TestCase):

    def setUp(self):
        m = 50
        n = 70
        p = 40
        self.X = np.random.randint(0, 10, size=(m + n) * p).reshape(((m + n), p))
        self.y = np.zeros((m+n, ))
        self.y[m:] = 1


    def test_seurat_alignment(self):
        """ Test on sample data with included scoring methods. """
        for method in GENE_SCORING_METHODS:
            model = SeuratAlignmentModel(n_components=10,
                                         n_metagenes=20,
                                         gene_scoring=method,
                                         random_state=42)
            Z = model.fit_transform(self.X, self.y)
            assert Z.shape == (self.X.shape[0], model.n_components)
            assert np.isnan(Z).sum() == 0


if __name__ == "__main__":
    unittest.main()
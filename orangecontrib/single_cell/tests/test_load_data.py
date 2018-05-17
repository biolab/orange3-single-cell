import os
import unittest

import numpy as np
import numpy.testing as npt
import pandas as pd

from orangecontrib.single_cell.widgets.load_data import (
    MtxLoader, CountLoader, Loader,
    get_data_loader
)


class TestLoadData(unittest.TestCase):
    def test_get_data_loader(self):
        dir_name = os.path.dirname(__file__)
        file_name = os.path.join(dir_name, "data/matrix.mtx")
        self.assertIsInstance(get_data_loader(file_name), MtxLoader)
        file_name = os.path.join(dir_name, "DATA_MATRIX_LOG_TPM.txt")
        self.assertIsInstance(get_data_loader(file_name), Loader)
        file_name = os.path.join(dir_name, "lib.cell.count")
        self.assertIsInstance(get_data_loader(file_name), CountLoader)

    def test_file_summary_mtx(self):
        file_name = os.path.join(os.path.dirname(__file__),
                                 "data/10x/mm10/matrix.mtx")
        loader = MtxLoader(file_name)
        self.assertEqual(loader.file_size, 105)
        self.assertEqual(loader.n_rows, 5)
        self.assertEqual(loader.n_cols, 5)

    def test_file_summary_broad(self):
        file_name = os.path.join(os.path.dirname(__file__),
                                 "data/DATA_MATRIX_LOG_TPM.txt")
        loader = Loader(file_name)
        self.assertEqual(loader.file_size, 1084)
        self.assertEqual(loader.n_rows, 10)
        self.assertEqual(loader.n_cols, 15)

    def test_file_summary_hhmi(self):
        file_name = os.path.join(os.path.dirname(__file__),
                                 "data/lib.cell.count")
        loader = CountLoader(file_name)
        self.assertEqual(loader.file_size, 428)
        self.assertEqual(loader.n_rows, 10)
        self.assertEqual(loader.n_cols, 11)

    def test_load_data_mtx(self):
        file_name = os.path.join(os.path.dirname(__file__),
                                 "data/10x/mm10/matrix.mtx")
        loader = MtxLoader(file_name)
        df = pd.read_csv(file_name, sep=" ", header=None, skiprows=[0, 1])
        attrs, X, meta_df, meta_df_index = loader._load_data()
        array = np.zeros((df.iloc[0, 0], df.iloc[0, 1]))
        for i, series in df.iterrows():
            if i == 0:
                continue
            array[series.iloc[1] - 1, series.iloc[0] - 1] = series.iloc[2]
        npt.assert_array_equal(X, array)

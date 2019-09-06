import os

import numpy as np
import numpy.testing as npt
import pandas as pd

from Orange.data import ContinuousVariable, Variable
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owloaddata import OWLoadData


class TestOWLoadData(WidgetTest):
    def setUp(self):
        Variable._clear_all_caches()
        self.widget = self.create_widget(OWLoadData)
        self._path = os.path.join(os.path.dirname(__file__), "data")

    def test_widget_no_file(self):
        self.assertEqual(self.widget.summary_label.text(), "")

    def test_open_mtx(self):
        file_name = os.path.join(self._path, "10x/mm10/matrix.mtx")
        self.widget.set_current_path(file_name)
        text = self.widget.summary_label.text()
        self.assertIn("5 rows, 5 columns", text)
        self._check_headers_and_row_labels_box((0, False, 0, False))
        self._check_input_data_structure_box((True, False, False, False))
        self._check_sample_data_box()
        values = ((True, False, "barcodes.tsv", False),
                  (True, False, "genes.tsv", False))
        self._check_annotation_files_box(values)

    def test_open_broad(self):
        file_name = os.path.join(self._path, "DATA_MATRIX_LOG_TPM.txt")
        self.widget.set_current_path(file_name)
        text = self.widget.summary_label.text()
        self.assertIn("10 rows, 15 columns", text)
        self._check_headers_and_row_labels_box((1, True, 1, True))
        self._check_input_data_structure_box((True, True, False, True))
        self._check_sample_data_box()
        values = ((False, True, "", False),
                  (False, True, "", False))
        self._check_annotation_files_box(values)

    def test_open_xls(self):
        file_name = os.path.join(self._path, "data.xlsx")
        self.widget.set_current_path(file_name)
        text = self.widget.summary_label.text()
        self.assertIn("11 rows, 15 columns", text)
        self._check_headers_and_row_labels_box((1, True, 1, True))
        self._check_input_data_structure_box((True, True, False, True))
        self._check_sample_data_box()
        values = ((False, True, "", False),
                  (False, True, "", False))
        self._check_annotation_files_box(values)

    def test_open_hhmi(self):
        file_name = os.path.join(self._path, "lib.cell.count")
        self.widget.set_current_path(file_name)
        text = self.widget.summary_label.text()
        self.assertIn("10 rows, 11 columns", text)
        self._check_headers_and_row_labels_box((1, False, 1, False))
        self._check_input_data_structure_box((True, False, False, False))
        self._check_sample_data_box()
        values = ((True, True, "lib.cell.meta", True),
                  (False, True, "", False))
        self._check_annotation_files_box(values)

    def test_load_data_mtx(self):
        file_name_mtx = os.path.join(self._path, "10x/hg19/matrix.mtx")
        self.widget.set_current_path(file_name_mtx)
        self.widget.commit()
        data = self.get_output("Data")
        file_name = os.path.join(self._path, "10x/hg19/genes.tsv")
        df = pd.read_csv(file_name, sep="\t", header=None)
        self._test_load_data_mtx_attributes(data.domain.attributes, df)
        df = pd.read_csv(file_name_mtx, sep=" ", header=None, skiprows=[0, 1])
        self._test_load_data_mtx_x(data.X, df)
        file_name = os.path.join(self._path, "10x/hg19/barcodes.tsv")
        df = pd.read_csv(file_name, sep="\t", header=None)
        self._test_load_data_metas(data.metas, df)

    def test_load_data_mtx_sample(self):
        file_name_mtx = os.path.join(self._path, "10x/mm10/matrix.mtx")
        self.widget.set_current_path(file_name_mtx)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(10)
        self.widget.set_sample_cols_p(10)
        self.widget.commit()
        data = self.get_output("Data")
        file_name = os.path.join(self._path, "10x/mm10/genes.tsv")
        df = pd.read_csv(file_name, sep="\t", header=None, skiprows=[4])
        self._test_load_data_mtx_attributes(data.domain.attributes, df)
        df = pd.read_csv(file_name_mtx, sep=" ", header=None,
                         skiprows=[0, 1, 9, 10])
        df.iloc[0, 0] = 4
        df.iloc[0, 1] = 4
        df.iloc[0, 2] = 6
        self._test_load_data_mtx_x(data.X, df)
        file_name = os.path.join(self._path, "10x/mm10/barcodes.tsv")
        df = pd.read_csv(file_name, sep="\t", header=None, skiprows=[4])
        self._test_load_data_metas(data.metas, df)

    def test_load_data_broad(self):
        file_name = os.path.join(self._path, "DATA_MATRIX_LOG_TPM.txt")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=0)
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_xls(self):
        file_name = os.path.join(self._path, "data.xlsx")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_excel(file_name, header=0, index_col=0)
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_compressed(self):
        file_name = os.path.join(self._path, "data.txt.gz")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=0)
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_pickle(self):
        file_name = os.path.join(self._path, "data.pkl")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        data = self.get_output("Data")
        file_name = os.path.join(self._path, "data.txt.gz")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=0)
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_loom(self):
        file_name = os.path.join(self._path, "data.loom")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        self.assertEqual(self.get_output("Data").X.shape, (20, 10))
        self.assertEqual(self.get_output("Data").metas.shape, (20, 1))

    def test_load_data_broad_sample(self):
        file_name = os.path.join(self._path, "DATA_MATRIX_LOG_TPM.txt")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(40)
        self.widget.set_sample_cols_p(60)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_csv(
            file_name, header=0, sep="\t", index_col=0,
            usecols=[0, 1, 2, 3, 5, 7, 12, 14], skiprows=[4, 5, 6, 7, 9, 10])
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_xls_sample(self):
        file_name = os.path.join(self._path, "data.xlsx")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(40)
        self.widget.set_sample_cols_p(60)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_excel(
            file_name, header=0, index_col=0,
            usecols=[0, 1, 2, 3, 5, 7, 12, 14], skiprows=[4, 5, 6, 7, 9, 10])
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        self._test_load_data_broad_metas(data.metas, df)

    def test_load_data_hhmi(self):
        file_name = os.path.join(self._path, "lib.cell.count")
        self.widget.set_current_path(file_name)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=0)
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        file_name = os.path.join(self._path, "lib.cell.meta")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=None)
        df.iloc[:, [1, 3]] = 0
        df.iloc[:, 6] = [0, 2, 0, 2, 0, 2, 0, 1, 0, 1]
        self._test_load_data_metas(data.metas, df)

    def test_load_data_hhmi_sample(self):
        file_name = os.path.join(self._path, "lib.cell.count")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(40)
        self.widget.set_sample_cols_p(60)
        self.widget.commit()
        data = self.get_output("Data")
        df = pd.read_csv(
            file_name, header=0, sep="\t", index_col=0,
            usecols=[0, 1, 2, 3, 5, 6], skiprows=[6, 8, 9, 10])
        self._test_load_data_attributes(data.domain.attributes, df)
        self._test_load_data_x(data.X, df)
        file_name = os.path.join(self._path, "lib.cell.meta")
        df = pd.read_csv(file_name, header=0, sep="\t", index_col=None,
                         skiprows=[4, 6, 8, 9, 10])
        df.iloc[:, [1, 3]] = 0
        df.iloc[:, 5] = [1, 0, 1, 1, 1]
        df.iloc[:, 6] = [0, 1, 0, 0, 0]
        self._test_load_data_metas(data.metas, df)

    def test_load_data_loom_sample(self):
        file_name = os.path.join(self._path, "data.loom")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(40)
        self.widget.set_sample_cols_p(60)
        self.widget.commit()
        data = self.get_output("Data")
        X = np.array([[4, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0], [0, 7, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 5, 0], [9, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0]])
        npt.assert_array_equal(data.X, X)
        metas = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 13, 16]]).T
        npt.assert_array_equal(data.metas, metas)
        self.assertListEqual([attr.name for attr in data.domain.attributes],
                             ["0", "1", "2", "3", "6", "7"])

    def test_load_data_h5ad_sparse_sample(self):
        file_name = os.path.join(self._path, "data_sparse.h5ad")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(30)
        self.widget.set_sample_cols_p(10)
        self.widget.commit()
        data = self.get_output("Data")
        X = np.array([[0, 2, 1, 0], [1, 0, 0, 0], [2, 1, 2, 1], [1, 1, 0, 0],
                      [1, 0, 0, 0], [0, 1, 0, 2], [1, 0, 0, 1], [0, 2, 0, 1],
                      [0, 1, 1, 0]])
        npt.assert_array_equal(data.X, X)
        metas = np.array([[0, 0, 10], [1, 0, 10], [2, 0, 15], [3, 0, 10],
                          [6, 1, 11], [7, 1, 15], [8, 1, 11], [9, 1, 20],
                          [14, 2, 12]])
        npt.assert_array_equal(data.metas, metas)
        self.assertListEqual([attr.name for attr in data.domain.attributes],
                             ["Gene 0", "Gene 1", "Gene 2", "Gene 3"])

    def test_load_data_h5ad_dense_sample(self):
        file_name = os.path.join(self._path, "data_dense.h5ad")
        self.widget.set_current_path(file_name)
        self.widget.sample_rows_cb.setChecked(True)
        self.widget.sample_cols_cb.setChecked(True)
        self.widget.set_sample_rows_p(30)
        self.widget.set_sample_cols_p(10)
        self.widget.commit()
        data = self.get_output("Data")
        X = np.array([[0, 2, 1, 0], [1, 0, 0, 0], [2, 1, 2, 1], [1, 1, 0, 0],
                      [1, 0, 0, 0], [0, 1, 0, 2], [1, 0, 0, 1], [0, 2, 0, 1],
                      [0, 1, 1, 0]])
        npt.assert_array_equal(data.X, X)
        metas = np.array([[0, 0, 10], [1, 0, 10], [2, 0, 15], [3, 0, 10],
                          [6, 1, 11], [7, 1, 15], [8, 1, 11], [9, 1, 20],
                          [14, 2, 12]])
        npt.assert_array_equal(data.metas, metas)
        self.assertListEqual([attr.name for attr in data.domain.attributes],
                             ["Gene 0", "Gene 1", "Gene 2", "Gene 3"])

    def test_not_enough_headers(self):
        file_name = os.path.join(self._path, "DATA_MATRIX_LOG_TPM.txt")
        self.widget.set_current_path(file_name)
        self.widget.set_header_rows_count(0)
        self.widget.commit()
        self.assertTrue(self.widget.Error.inadequate_headers.is_shown())
        self.widget.set_header_cols_count(0)
        self.widget.commit()
        self.assertTrue(self.widget.Error.inadequate_headers.is_shown())
        self.widget.set_header_rows_count(1)
        self.widget.set_header_cols_count(1)
        self.widget.commit()
        self.assertFalse(self.widget.Error.inadequate_headers.is_shown())

    def _check_headers_and_row_labels_box(self, values):
        self.assertEqual(self.widget.header_rows_spin.value(), values[0])
        self.assertEqual(self.widget.header_rows_spin.isEnabled(), values[1])
        self.assertEqual(self.widget.header_cols_spin.value(), values[2])
        self.assertEqual(self.widget.header_cols_spin.isEnabled(), values[3])

    def _check_input_data_structure_box(self, values):
        rbs = self.widget.data_struct_box.children()[1].buttons
        self.assertEqual(rbs[0].isChecked(), values[0])
        self.assertEqual(rbs[0].isEnabled(), values[1])
        self.assertEqual(rbs[1].isChecked(), values[2])
        self.assertEqual(rbs[1].isEnabled(), values[3])

    def _check_sample_data_box(self):
        self.assertFalse(self.widget.sample_rows_cb.isChecked())
        self.assertTrue(self.widget.sample_rows_cb.isEnabled())
        self.assertFalse(self.widget.sample_cols_cb.isChecked())
        self.assertTrue(self.widget.sample_cols_cb.isEnabled())

    def _check_annotation_files_box(self, v):
        self.assertEqual(self.widget.row_annotations_cb.isChecked(), v[0][0])
        self.assertEqual(self.widget.row_annotations_cb.isEnabled(), v[0][1])
        text = self.widget.row_annotations_combo.currentText()
        self.assertEqual(text, v[0][2])
        self.assertEqual(self.widget.row_annotations_combo.isEnabled(), v[0][3])
        self.assertEqual(self.widget.col_annotations_cb.isChecked(), v[1][0])
        self.assertEqual(self.widget.col_annotations_cb.isEnabled(), v[1][1])
        text = self.widget.col_annotations_combo.currentText()
        self.assertEqual(text, v[1][2])
        self.assertEqual(self.widget.col_annotations_combo.isEnabled(), v[1][3])

    def _test_load_data_attributes(self, attributes, df):
        for i, attribute in enumerate(attributes):
            self.assertIsInstance(attribute, ContinuousVariable)
            self.assertEqual(attribute.name, df.index[i])
            self.assertDictEqual(attribute.attributes, {})

    def _test_load_data_x(self, x, df):
        npt.assert_array_equal(x, df.values.T)

    def _test_load_data_metas(self, metas, df):
        npt.assert_array_equal(metas, df.values)

    def _test_load_data_mtx_attributes(self, attributes, df):
        for i, attribute in enumerate(attributes):
            self.assertIsInstance(attribute, ContinuousVariable)
            self.assertEqual(attribute.name, df.iloc[i, 0])
            self.assertEqual(attribute.attributes["Id"], df.iloc[i, 0])
            self.assertEqual(attribute.attributes["Gene"], df.iloc[i, 1])

    def _test_load_data_mtx_x(self, x, df):
        array = np.zeros((df.iloc[0, 1], df.iloc[0, 0]))
        for i, series in df.iterrows():
            if i == 0:
                continue
            array[series.iloc[1] - 1, series.iloc[0] - 1] = series.iloc[2]
        npt.assert_array_equal(x, array)

    def _test_load_data_broad_metas(self, metas, data_frame):
        npt.assert_array_equal(metas.flatten(), data_frame.columns.values)

    # def test_valid_output(self):
    #     widget = self.create_widget(OWDataTable)
    #     self.widget.set_current_path(os.path.join(self._path, "data.xlsx"))
    #     self.widget.set_header_rows_count(3)
    #     self.widget.set_header_cols_count(2)
    #     self.widget.sample_rows_cb.setChecked(True)
    #     self.widget.sample_cols_cb.setChecked(True)
    #     self.widget.set_sample_rows_p(40)
    #     self.widget.set_sample_cols_p(60)
    #     self.widget.commit()
    #     data = self.get_output("Data")
    #     self.assertIsNotNone(data)
    #     self.assertEqual(data.X.shape, (3, 7))
    #     widget.set_dataset(data)

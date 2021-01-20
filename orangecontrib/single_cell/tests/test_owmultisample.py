import os
import unittest
from unittest.mock import patch, Mock
from pathlib import PurePath

import numpy.testing as npt

from AnyQt.QtCore import Qt, QMimeData, QUrl, QPoint
from AnyQt.QtGui import QDropEvent

from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owmultisample import OWMultiSample


class TestOWMultiSample(WidgetTest):
    def create_widget(self, cls, stored_settings=None,
                      reset_default_settings=True):
        if reset_default_settings:
            self.reset_default_settings(cls)
        widget = cls.__new__(cls, signal_manager=self.signal_manager,
                             stored_settings=stored_settings)
        widget.__init__()
        self.process_events()
        self.widgets.append(widget)
        return widget

    def setUp(self):
        self.widget = self.create_widget(
            OWMultiSample, stored_settings={"auto_commit": True}
        )
        self._path = path = os.path.join(os.path.dirname(__file__), "data")
        self.file_name_1 = os.path.join(path, "10x/hg19/matrix.mtx")
        self.file_name_2 = os.path.join(path, "10x/mm10/matrix.mtx")
        self.widget.set_current_path(self.file_name_1)
        self.widget.set_current_path(self.file_name_2)
        model = self.widget.view.model()
        model.item(0).setCheckState(True)
        model.item(1).setCheckState(True)
        self.widget.commit()

    def test_minimum_size(self):
        pass

    def test_load_samples(self):
        self.assertEqual(self.widget.view.model().rowCount(), 2)

    def test_concatenate_intersection_mtx(self):
        concatenated_data = self.get_output("Data")
        domain = concatenated_data.domain
        self.assertEqual(len(concatenated_data), 11)
        self.assertEqual(len(domain.attributes), 1)
        self.assertEqual(len(domain.metas), 2)

    def test_concatenate_union_mtx(self):
        self.widget.controls.output_type.buttons[1].click()
        concatenated_data = self.get_output("Data")
        domain = concatenated_data.domain
        self.assertEqual(len(concatenated_data), 11)
        self.assertEqual(len(domain.attributes), 8)
        self.assertEqual(len(domain.metas), 2)
        self.assertTrue(all(
            [list(attr.attributes.keys()) == ["Id", "Gene"]
             for attr in domain.attributes]
        ))

    def test_settings(self):
        self.widget.saveSettings()
        widget = self.create_widget(
            OWMultiSample, reset_default_settings=False
        )
        self.assertEqual(widget.view.model().rowCount(), 2)
        self.assertIsNotNone(self.get_output("Data", widget))
        npt.assert_array_equal(
            self.get_output("Data").X, self.get_output("Data", widget).X
        )

    def test_drop_sample(self):
        path = os.path.join(self._path, "lib.cell.count")
        data = QMimeData()
        data.setUrls([QUrl.fromLocalFile(path)])
        event = QDropEvent(
            QPoint(10, 10), Qt.MoveAction, data,
            Qt.NoButton, Qt.NoModifier, QDropEvent.Drop)
        self.widget.view.dropEvent(event)
        self.assertEqual(self.widget.view.model().rowCount(), 3)

    @patch("Orange.widgets.widget.OWWidget.workflowEnv",
           Mock(return_value={"basedir": os.path.join(os.path.dirname(__file__), "data")}))
    def test_load_workflow(self):
        base_path = os.path.join("010-Showcase-LoadDataM-Data", "batch-A")
        base_name = str(PurePath(os.path.join(base_path, "matrix.mtx")))
        base_row_name = str(PurePath(os.path.join(base_path, "barcodes.tsv")))
        base_col_name = str(PurePath(os.path.join(base_path, "genes.tsv")))

        # store settings - simulate saving workflow
        w1 = self.create_widget(OWMultiSample)
        w1.set_current_path(os.path.join(os.path.join(os.path.dirname(__file__), "data"), str(base_name)))
        w1.write_settings()
        settings = self.widget.settingsHandler.pack_data(w1)

        # check other widget settings - simulate opening workflow
        w2 = self.create_widget(OWMultiSample, stored_settings=settings)
        loader = w2._data_loader

        self.assertEqual(str(PurePath(w2._recent[0].relpath)), base_name)
        self.assertEqual(w2._recent[0].prefix, "basedir")
        self.assertEqual(str(PurePath(loader.recent_path.relpath)), base_name)
        self.assertEqual(loader.recent_path.prefix, "basedir")
        self.assertEqual(str(PurePath(loader.row_annotation_file.relpath)), base_row_name)
        self.assertEqual(loader.row_annotation_file.prefix, "basedir")
        self.assertEqual(str(PurePath(loader.col_annotation_file.relpath)), base_col_name)
        self.assertEqual(loader.col_annotation_file.prefix, "basedir")


if __name__ == "__main__":
    unittest.main()

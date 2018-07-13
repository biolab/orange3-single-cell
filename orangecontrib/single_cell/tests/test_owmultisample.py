import os

import numpy.testing as npt
from AnyQt.QtCore import Qt, QMimeData, QUrl, QPoint
from AnyQt.QtGui import QDropEvent

from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owmultisample import OWMultiSample


class TestOWMultiSample(WidgetTest):
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

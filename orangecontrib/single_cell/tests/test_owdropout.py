import unittest
from unittest.mock import Mock, patch

from pyqtgraph import PlotCurveItem

from AnyQt.QtCore import Qt, QPoint
from AnyQt.QtTest import QTest

from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.widget import OWWidget

from orangecontrib.single_cell.widgets.owdropout import OWDropout, DropoutGraph


class DummyWidget(OWWidget):
    name = "dummy"


class DummyPoint:
    xi = iter(range(10))
    yi = iter(range(10))

    def x(self):
        return next(self.xi)

    def y(self):
        return next(self.yi) * 2


class TestDropoutGraph(WidgetTest):
    def setUp(self):
        self.parent = DummyWidget()
        self.graph = DropoutGraph(self.parent)
        self.results = results = Mock()
        results.decay = 1
        results.x_offset = 0.1
        results.y_offset = 0.1
        results.mean_expr = [1, 2, 3]
        results.zero_rate = [1, 2, 3]
        results.threshold = 0

    def test_on_curve(self):
        self.assertFalse(self.graph._on_curve(1, 0.5))
        self.graph.set_data(self.results)
        self.assertTrue(self.graph._on_curve(1, 0.5))
        self.assertFalse(self.graph._on_curve(1, 1))

    @patch("orangecontrib.single_cell.widgets.owdropout.DropoutGraph."
           "_on_curve", Mock(return_value=True))
    @patch.object(PlotCurveItem, "mapFromScene",
                  Mock(return_value=DummyPoint()))
    def test_cursor_event(self):
        graph = self.graph
        slot = Mock()
        graph.curve_moved.connect(slot)
        graph.set_data(self.results)
        graph.is_curve_movable = True
        graph.cursor_event(Mock())
        QTest.mousePress(graph.viewport(), Qt.LeftButton, pos=QPoint(1, 1))
        # QTest.mouseMove(g.viewport()) does not work; set flag manually
        graph._state = 3
        graph.cursor_event(Mock())
        slot.assert_called_once_with(1, 2)


class TestOWDropout(WidgetTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")

    def setUp(self):
        self.widget = self.create_widget(OWDropout)

    def test_controls_enabled(self):
        controls = self.widget.controls
        buttons = controls.filter_type.buttons
        self.assertEqual(self.widget.filter_type, 0)
        self.assertTrue(controls.n_genes.isEnabled())
        self.assertFalse(controls.decay.isEnabled())
        self.assertFalse(controls.x_offset.isEnabled())
        self.assertFalse(controls.y_offset.isEnabled())

        buttons[1].click()
        self.assertEqual(self.widget.filter_type, 1)
        self.assertFalse(controls.n_genes.isEnabled())
        self.assertTrue(controls.decay.isEnabled())
        self.assertTrue(controls.x_offset.isEnabled())
        self.assertTrue(controls.y_offset.isEnabled())

        buttons[2].click()
        self.assertEqual(self.widget.filter_type, 2)
        self.assertFalse(controls.n_genes.isEnabled())
        self.assertFalse(controls.decay.isEnabled())
        self.assertFalse(controls.x_offset.isEnabled())
        self.assertFalse(controls.y_offset.isEnabled())

    def test_manual_move_enabled(self):
        controls = self.widget.controls
        buttons = controls.filter_type.buttons
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertFalse(self.widget.graph.is_curve_movable)
        buttons[2].click()
        self.assertTrue(self.widget.graph.is_curve_movable)
        self.send_signal(self.widget.Inputs.data, None)
        self.assertTrue(self.widget.graph.is_curve_movable)
        buttons[1].click()
        self.assertFalse(self.widget.graph.is_curve_movable)
        buttons[2].click()
        self.assertTrue(self.widget.graph.is_curve_movable)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertTrue(self.widget.graph.is_curve_movable)

    def test_output(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        output = self.get_output(self.widget.Outputs.data)
        self.assertEqual(len(output.domain.attributes), 497)
        self.assertTrue(self.widget.Warning.less_selected.is_shown())
        self.send_signal(self.widget.Inputs.data, None)
        self.assertIsNone(self.get_output(self.widget.Outputs.data))
        self.assertFalse(self.widget.Warning.less_selected.is_shown())

    def test_info(self):
        self.assertEqual(self.widget.info_label.text(), "No data on input.")
        self.send_signal(self.widget.Inputs.data, self.data)
        text = "Data with 1000 cells and 1000 genes\n497 genes in selection"
        self.assertEqual(self.widget.info_label.text(), text)
        self.send_signal(self.widget.Inputs.data, self.data[:1])
        text = "Data with 1 cell and 1000 genes\n0 genes in selection"
        self.assertEqual(self.widget.info_label.text(), text)
        self.send_signal(self.widget.Inputs.data, None)
        self.assertEqual(self.widget.info_label.text(), "No data on input.")

    def test_n_genes_changed(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.controls.n_genes.setValue(200)
        output = self.get_output(self.widget.Outputs.data)
        self.assertEqual(len(output.domain.attributes), 200)

    @unittest.skip("Skip due to numpy.linalg.LinAlgError: Singular matrix")
    def test_send_report(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.report_button.click()
        self.widget.controls.filter_type.buttons[1].click()
        self.widget.report_button.click()
        self.send_signal(self.widget.Inputs.data, None)
        self.widget.report_button.click()


if __name__ == "__main__":
    unittest.main()

import unittest
from unittest.mock import Mock, patch

from pyqtgraph import PlotCurveItem

from AnyQt.QtCore import Qt, QPoint
from AnyQt.QtTest import QTest

from Orange.data import Table
from Orange.data.filter import Values, FilterString
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.widget import OWWidget

from orangecontrib.bioinformatics.utils import serverfiles

from orangecontrib.single_cell.widgets.owdropout import OWDropout, \
    DropoutGraph, FilterType


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
        self.data = Table("iris")[:3]

    def test_on_curve(self):
        self.assertFalse(self.graph._on_curve(1, 0.5))
        self.graph.set_data(self.results, self.data, None)
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
        graph.set_data(self.results, self.data, None)
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
        genes_path = serverfiles.localpath_download("marker_genes",
                                                    "panglao_gene_markers.tab")
        filter_ = FilterString("Organism", FilterString.Equal, "Human")
        cls.genes = Values([filter_])(Table(genes_path))
        cls.iris = Table("iris")

    def setUp(self):
        self.widget = self.create_widget(OWDropout)

    def test_input_data(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.send_signal(self.widget.Inputs.data, self.iris)
        self.send_signal(self.widget.Inputs.data, self.genes)
        self.send_signal(self.widget.Inputs.data, None)

    @patch("orangecontrib.single_cell.widgets.owdropout."
           "DropoutGraph.update_markers")
    @patch("orangecontrib.single_cell.widgets.owdropout."
           "DropoutGraph.set_data")
    def test_input_genes(self, set_data: Mock, update_markers: Mock):
        # genes
        self.send_signal(self.widget.Inputs.genes, self.genes)
        update_markers.assert_called_once_with(None, self.genes)
        update_markers.reset_mock()
        set_data.assert_not_called()

        # data, genes
        self.send_signal(self.widget.Inputs.data, self.data)
        update_markers.assert_not_called()
        update_markers.reset_mock()
        set_data.assert_called_once()

        # data
        self.send_signal(self.widget.Inputs.genes, None)
        update_markers.assert_called_once_with(self.data, None)
        update_markers.reset_mock()
        set_data.assert_called_once()

        # data, genes
        self.send_signal(self.widget.Inputs.genes, self.iris)
        self.assertTrue(self.widget.Warning.missing_entrez_id.is_shown())
        update_markers.assert_called_once_with(self.data, None)
        update_markers.reset_mock()
        set_data.assert_called_once()

        # data
        self.send_signal(self.widget.Inputs.genes, None)
        self.assertFalse(self.widget.Warning.missing_entrez_id.is_shown())
        update_markers.assert_called_once_with(self.data, None)
        set_data.assert_called_once()

    def test_controls_enabled(self):
        controls = self.widget.controls
        buttons = controls.filter_type.buttons
        self.assertEqual(self.widget.filter_type, FilterType.ByNumber)
        self.assertTrue(controls.n_genes.isEnabled())
        self.assertFalse(controls.decay.isEnabled())
        self.assertFalse(controls.x_offset.isEnabled())
        self.assertFalse(controls.y_offset.isEnabled())

        buttons[1].click()
        self.assertEqual(self.widget.filter_type, FilterType.ByEquation)
        self.assertFalse(controls.n_genes.isEnabled())
        self.assertTrue(controls.decay.isEnabled())
        self.assertTrue(controls.x_offset.isEnabled())
        self.assertTrue(controls.y_offset.isEnabled())

    def test_manual_move(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.graph.curve_moved.emit(1, 2)
        self.assertEqual(self.widget.filter_type, FilterType.ByEquation)
        controls = self.widget.controls
        self.assertTrue(controls.filter_type.buttons[1].clicked)
        self.assertFalse(controls.n_genes.isEnabled())
        self.assertTrue(controls.decay.isEnabled())
        self.assertTrue(controls.x_offset.isEnabled())
        self.assertTrue(controls.y_offset.isEnabled())

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


if __name__ == "__main__":
    unittest.main()

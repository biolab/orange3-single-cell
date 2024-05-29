import unittest
from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest
from orangecontrib.single_cell.widgets.owspatialautocorrelation import OWSpatialAutocorrelation
from Orange.widgets.utils.concurrent import TaskState


class TestOWSpatialAutocorrelation(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSpatialAutocorrelation)
        
    def test_input_data(self):
        data = Table("iris")
        self.widget.set_data(data)
        self.assertEqual(self.widget.data, data)

    def test_feature_selection(self):
        data = Table("iris")
        self.widget.set_data(data)
        self.widget.feature_x_combo.setCurrentIndex(1)
        self.widget.feature_y_combo.setCurrentIndex(2)
        self.assertEqual(self.widget.feature_x, "sepal width")
        self.assertEqual(self.widget.feature_y, "petal length")

    def test_method_selection(self):
        self.widget.method = "Geary C"
        self.assertEqual(self.widget.method, "Geary C")

    def test_k_neighbors_input(self):
        self.widget.k_input.setText("10")
        self.widget._on_k_changed()
        self.assertEqual(self.widget.k_neighbors, 10)

    def test_auto_commit(self):
        self.widget.auto_commit = False
        self.assertFalse(self.widget.auto_commit)

    def test_calculate(self):
        data = Table("iris")
        self.widget.set_data(data)
        self.widget.feature_x_combo.setCurrentIndex(0)
        self.widget.feature_y_combo.setCurrentIndex(1)
        self.widget.k_input.setText("5")
        self.widget._on_k_changed()
        self.widget.calculate(TaskState())
        self.assertIsNotNone(self.widget.adjacency_matrix)


if __name__ == "__main__":
    unittest.main()
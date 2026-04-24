from AnyQt.QtCore import Qt
from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest
from orangecontrib.single_cell.widgets.owharmony import OWHarmony


class TestOWBatchNorm(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWHarmony)

    def test_batch_var_change(self):
        data = Table("iris")
        w = self.widget
        self.send_signal(w.Inputs.data, data)

        self.widget.model.item(0).setCheckState(Qt.CheckState.Checked)
        self.assertEqual(w.batch_vars, [data.domain["iris"]])
        output = self.get_output(w.Outputs.data)
        self.assertNotEqual(output.X[0, 0], data.X[0, 0])

        self.widget.model.item(0).setCheckState(Qt.CheckState.Unchecked)
        self.assertEqual(w.batch_vars, [])
        output = self.get_output(w.Outputs.data)
        self.assertEqual(output.X[0, 0], data.X[0, 0])

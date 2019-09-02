import unittest

from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owdropout import OWDropout


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

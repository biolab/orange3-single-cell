import numpy as np

from AnyQt.QtCore import Qt

from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.tests.utils import simulate

from orangecontrib.single_cell.widgets.owbatchnorm import (
    OWBatchNorm, LinkMethod
)


class TestOWBatchNorm(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWBatchNorm)

    def test_valid_input(self):
        """Check batch normalization for data set with continuous attributes"""
        data = Table("iris")
        self.send_signal(self.widget.Inputs.data, data)
        self.assertEqual(self.widget.model.rowCount(), 1)
        self.widget.model.item(0).setCheckState(Qt.Checked)
        output = self.get_output(self.widget.Outputs.data)
        self.assertIsInstance(output, Table)
        self.assertFalse((output.X == data.X).any())

    def test_discrete_attributes_input(self):
        """Check batch normalization for data set with discrete attributes"""
        self.send_signal(self.widget.Inputs.data, Table("heart_disease"))
        self.assertIsNone(self.get_output(self.widget.Outputs.data))
        self.assertTrue(self.widget.Error.discrete_attributes.is_shown())
        self.send_signal(self.widget.Inputs.data, None)
        self.assertFalse(self.widget.Error.discrete_attributes.is_shown())
        self.assertIsNone(self.get_output(self.widget.Outputs.data))
        self.assertEqual(self.widget.info_label.text(), "No data on input.")

    def test_missing_values_input(self):
        """Check batch normalization for data set with unknown values """
        data = Table("iris")
        data[0, 3] = np.nan
        self.send_signal(self.widget.Inputs.data, data)
        self.assertIsInstance(self.get_output(self.widget.Outputs.data), Table)
        self.assertTrue(self.widget.Warning.missing_values.is_shown())
        self.send_signal(self.widget.Inputs.data, None)
        self.assertFalse(self.widget.Warning.missing_values.is_shown())
        self.assertIsNone(self.get_output(self.widget.Outputs.data))

    def test_negative_values_input_log_link(self):
        """Check batch normalization for data set with negative values"""
        data = Table("iris")
        data[0, 3] = -1
        self.send_signal(self.widget.Inputs.data, data)
        self.widget.model.item(0).setCheckState(Qt.Checked)
        link_method = self.widget.controls.link_method
        simulate.combobox_activate_index(link_method, LinkMethod.LOG_LINK)

        self.assertTrue(self.widget.Warning.negative_values.is_shown())
        output = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(output.X, data.X)

        self.send_signal(self.widget.Inputs.data, None)
        self.assertFalse(self.widget.Warning.negative_values.is_shown())

    def test_negative_values_input_id_link(self):
        """Check batch normalization for data set with negative values"""
        data = Table("iris")
        data[0, 3] = -1
        self.send_signal(self.widget.Inputs.data, data)
        self.widget.model.item(0).setCheckState(Qt.Checked)
        link_method = self.widget.controls.link_method
        simulate.combobox_activate_index(link_method, LinkMethod.IDENTITY_LINK)

        self.widget.skip_zeros_check.setChecked(True)
        self.assertTrue(self.widget.Warning.negative_values.is_shown())
        output = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(output.X, data.X)

        self.widget.skip_zeros_check.setChecked(False)
        self.assertFalse(self.widget.Warning.negative_values.is_shown())
        output = self.get_output(self.widget.Outputs.data)
        self.assertFalse((output.X == data.X).any())

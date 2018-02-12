import unittest
import numpy as np

from orangecontrib.single_cell.widgets.ownormalization import OWNormalization
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWNormalization(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWNormalization)

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]

        domain = Domain(attributes, class_vars=class_var)

        self.data = Table(domain, [[0, 1, 3, 7, 15, 'STG1']])
        self.expected_base2_transform = np.array([[0., 1., 2., 3., 4.]])

    def test_no_input(self):
        self.send_signal(self.widget.Inputs.data, None)

        self.assertEqual(self.widget.info.text(), 'No data on input')
        self.assertIsNone(self.widget.data)

    def test_log_transformation(self):
        # test if input is OK
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        # test transformation results
        self.widget.log_base = 2  # ensure base2
        self.widget.commit()
        np.testing.assert_array_equal(self.get_output(self.widget.Outputs.data).X,
                                      self.expected_base2_transform)


if __name__ == '__main__':
    unittest.main()

import unittest
import numpy as np

from orangecontrib.single_cell.widgets.owscoregenes import OWRank
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWRank(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWRank)
        self.widget.selected_methods = set()

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(2)]
        domain = Domain(attributes, class_vars=class_var)

        self.data = Table(domain, [[8, -4, 'STG1'],
                                   [9, -2, 'STG1'],
                                   [10, 0, 'STG2'],
                                   [11, 2, 'STG2'],
                                   [12, 4, 'STG2']])

        self.expected_mean = np.array([[10.], [0.]])
        self.expected_variance = np.array([[2.], [8.]])
        self.expected_dispersion = np.array([[0.2], [8.]])

        # GeneName1(second column) -> if mean is zero, we set it to 1
        self.expected_cv = np.array([[0.141], [2.828]])

    def test_mean(self):
        self.widget.selected_methods.add('Mean')

        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        out_data = self.get_output(self.widget.Outputs.scores)
        self.assertTrue((out_data.X == self.expected_mean).all())

    def test_variance(self):
        self.widget.selected_methods.add('Variance')

        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        out_data = self.get_output(self.widget.Outputs.scores)
        self.assertTrue((out_data.X == self.expected_variance).all())

    def test_dispersion(self):
        self.widget.selected_methods.add('Dispersion')

        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        out_data = self.get_output(self.widget.Outputs.scores)
        self.assertTrue((out_data.X == self.expected_dispersion).all())

    def test_coef_of_variation(self):
        self.widget.selected_methods.add('Coef. of variation (CV)')

        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        out_data = self.get_output(self.widget.Outputs.scores)
        self.assertTrue((np.around(out_data.X, decimals=3) == self.expected_cv).all())


if __name__ == '__main__':
    unittest.main()
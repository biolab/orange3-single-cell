import unittest
import numpy as np

from orangecontrib.single_cell.widgets.owclusterstatistics import OWClusterStatistics
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWClusterStatistics(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWClusterStatistics)

        class_var = DiscreteVariable('Cluster', values=['C1', 'C2', 'C3'])
        domain = Domain([ContinuousVariable('GeneName')], class_vars=class_var)

        self.data = Table(domain, [[1, 'C1'],
                                   [1, 'C1'],
                                   [1, 'C1'],
                                   [2, 'C2'],
                                   [2, 'C2'],
                                   [3, 'C3']])

        self.expected_cluster_sizes = np.array([3, 2, 1])
        self.expected_sum = np.array([[3.], [4.], [3.]])

    def test_no_data(self):
        # input data
        self.send_signal(self.widget.Inputs.data, None)

        self.assertIsNone(self.widget.data)
        self.assertEqual(self.widget.info.text(), 'No data on input')

    def test_no_discrete_var(self):
        # input data
        self.send_signal(self.widget.Inputs.data, Table([[1, 2, 3]]))
        self.assertIsNotNone(self.widget.data)
        # test discrete var
        self.assertTrue(self.widget.Warning.no_discrete_attributes)

    def test_cluster_statistics(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data)

        # we have 3 clusters in data
        self.assertTrue(len(out_data) == 3)

        # check cluster sizes
        size_column = out_data.get_column_view(out_data.columns.size)[0]
        self.assertTrue((self.expected_cluster_sizes == size_column).all())

        # check sum
        self.assertTrue((out_data.X == self.expected_sum).all())


if __name__ == '__main__':
    unittest.main()

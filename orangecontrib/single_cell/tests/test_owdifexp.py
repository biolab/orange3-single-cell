import unittest

from orangecontrib.single_cell.widgets.owdifexp import OWDifferentialExpression
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWDifferentialExpression(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWDifferentialExpression)

        class_var = DiscreteVariable('Cluster', values=['C1', 'C2', 'C3'])
        attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]
        domain = Domain(attributes, class_vars=class_var)

        self.data = Table(domain, [[2, 2, 2, 2, 2, 'C1'],
                                   [4, 4, 4, 4, 4, 'C1'],
                                   [6, 6, 6, 6, 6, 'C2'],
                                   [8, 8, 8, 8, 8, 'C2']])

    def test_dif_expression(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.differential_expression)

        # check if correctly calculated
        # we expect genes from C1 to have values less than 0
        # we expect genes from C2 to have values greater than 0
        self.assertTrue((0 > out_data.X[:, 0]).all())
        self.assertTrue((0 < out_data.X[:, 1]).all())


if __name__ == '__main__':
    unittest.main()

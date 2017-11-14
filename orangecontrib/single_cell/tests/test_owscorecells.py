import unittest
import numpy as np

from orangecontrib.single_cell.widgets.owscorecells import OWScoreCells
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, StringVariable, Domain, Table


class TestOWScoreCells(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWScoreCells)

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]
        domain = Domain(attributes, class_vars=class_var)

        self.data = Table(domain, [[1, 2, 3, 4, 5, 'STG1'],
                                   [4, 4, 4, 4, 4, 'STG1'],
                                   [2, 3, 1, 1, 1, 'STG1'],
                                   [-1, 0, 1, 0, 0, 'STG2']])

        self.genes = Table(Domain([], metas=[StringVariable('Name')]), [[gene] for gene in self.data.domain.attributes])
        self.expected_score_values = np.array([[5.], [4.], [3.], [1.]])

    def test_score_cells(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        # input genes
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.assertIsNotNone(self.widget.genes)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data)
        self.assertTrue((self.expected_score_values == out_data.metas).all())

    def test_no_input_genes(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.data)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data)
        self.assertIsNone(out_data)


if __name__ == '__main__':
    unittest.main()

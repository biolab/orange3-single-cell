import unittest
import numpy as np

from orangecontrib.single_cell.widgets.owscorecells import OWScoreCells
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, StringVariable, Domain, Table

from orangecontrib.bioinformatics.widgets.utils.data import (
    GENE_ID_ATTRIBUTE,
    GENE_ID_COLUMN,
    GENE_AS_ATTRIBUTE_NAME,
    TAX_ID,
)


class TestOWScoreCells(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWScoreCells)

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]

        for idx, attr in enumerate(attributes):
            attr.attributes['Entrez ID'] = idx

        domain = Domain(attributes, class_vars=class_var)

        self.data = Table.from_numpy(
            domain,
            [[1, 2, 3, 4, 5], [4, 4, 4, 4, 4], [2, 3, 1, 1, 1], [-1, 0, 1, 0, 0]],
            [[0], [0], [0], [1]]
        )
        self.data.attributes[TAX_ID] = '9606'
        self.data.attributes[GENE_AS_ATTRIBUTE_NAME] = True
        self.data.attributes[GENE_ID_ATTRIBUTE] = 'Entrez ID'

        self.genes = Table.from_list(
            Domain([], metas=[StringVariable('Entrez ID')]),
            [[attr.attributes.get('Entrez ID')] for attr in self.data.domain.attributes],
        )

        self.genes.attributes[TAX_ID] = '9606'
        self.genes.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        self.genes.attributes[GENE_ID_COLUMN] = 'Entrez ID'

        self.expected_score_values = np.array([[3.0], [4.0], [1.6], [0]])

    def test_score_cells(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.input_data)
        self.assertIsNotNone(self.widget.input_genes)

        # marker data
        self.send_signal(self.widget.Inputs.genes, self.genes)

        self.assertIsNotNone(self.widget.marker_data)
        self.assertIsNotNone(self.widget.marker_genes)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data)
        self.assertTrue((self.expected_score_values == out_data.metas).all())

    def test_no_input_genes(self):
        # input data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNotNone(self.widget.input_data)

        # output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data)
        self.assertTrue(np.all(out_data.get_column_view('Score')[0] == 1.0))


if __name__ == '__main__':
    unittest.main()

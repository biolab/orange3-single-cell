import unittest

from Orange.data import Table, Domain, StringVariable
from Orange.widgets.tests.base import WidgetTest
from orangecontrib.bioinformatics.ncbi.gene import NCBI_ID
from orangecontrib.bioinformatics.widgets.utils.data import GENE_AS_ATTRIBUTE_NAME, GENE_ID_ATTRIBUTE, GENE_ID_COLUMN

from orangecontrib.single_cell.widgets.owclusteranalysis import OWClusterAnalysis


class TestOWClusterAnalysis(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWClusterAnalysis)

        self.data_table = Table("iris")
        self.data_table.attributes[GENE_AS_ATTRIBUTE_NAME] = True
        self.data_table.attributes[GENE_ID_ATTRIBUTE] = NCBI_ID
        for i, var in enumerate(self.data_table.domain.attributes):
            var.attributes[NCBI_ID] = str(i)

        domain = Domain(self.data_table.domain.attributes[0:2])
        self.genes_as_attributes = self.data_table.transform(domain)

        domain = Domain([], metas=[StringVariable("Gene ID")])
        self.genes_as_rows = Table.from_list(domain, [["1"], ["2"]])
        self.genes_as_rows.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        self.genes_as_rows.attributes[GENE_ID_COLUMN] = "Gene ID"

    def tearDown(self):
        self.widget.onDeleteWidget()
        super().tearDown()

    def test_input(self):
        """
        Test combinations of input and whether outputs appear at all.
        """
        # data input
        self.send_signal(self.widget.Inputs.data, self.data_table)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertIsNone(self.get_output(self.widget.Outputs.selected_data))
        self.assertIsNotNone(self.get_output(self.widget.Outputs.annotated_data))
        self.assertIsNotNone(self.get_output(self.widget.Outputs.contingency))

        # remove data input
        self.send_signal(self.widget.Inputs.data, None)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertIsNone(self.get_output(self.widget.Outputs.selected_data))
        self.assertIsNone(self.get_output(self.widget.Outputs.annotated_data))
        self.assertIsNone(self.get_output(self.widget.Outputs.contingency))

        # only genes
        self.send_signal(self.widget.Inputs.genes, self.genes_as_attributes)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertIsNone(self.get_output(self.widget.Outputs.selected_data))
        self.assertIsNone(self.get_output(self.widget.Outputs.annotated_data))
        self.assertIsNone(self.get_output(self.widget.Outputs.contingency))

        # data input and genes input
        self.send_signal(self.widget.Inputs.data, self.data_table)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertIsNone(self.get_output(self.widget.Outputs.selected_data))
        self.assertIsNotNone(self.get_output(self.widget.Outputs.annotated_data))
        self.assertIsNotNone(self.get_output(self.widget.Outputs.contingency))

    def test_gene_list(self):
        """
        Test whether genes input constructs a correct contingency table.
        """
        self.send_signal(self.widget.Inputs.data, self.data_table)
        self.process_events(until=lambda: self.widget._task is None)

        self.send_signal(self.widget.Inputs.genes, self.genes_as_attributes)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertEqual(len(self.genes_as_attributes.domain.attributes),
                         len(self.get_output(self.widget.Outputs.contingency).domain.attributes))

        self.send_signal(self.widget.Inputs.genes, self.genes_as_rows)
        self.process_events(until=lambda: self.widget._task is None)
        self.assertEqual(len(self.genes_as_rows),
                         len(self.get_output(self.widget.Outputs.contingency).domain.attributes))

    def test_filtering(self):
        """
        Test whether output data gets filtered correctly.
        """
        self.send_signal(self.widget.Inputs.data, self.data_table)
        self.process_events(until=lambda: self.widget._task is None)

        self.widget.tableview.set_selection({})
        self.assertIsNone(self.get_output(self.widget.Outputs.selected_data))

        self.widget.selection = {(0, 0)}
        self.widget.commit()
        self.assertEqual(50, len(self.get_output(self.widget.Outputs.selected_data)))
        self.assertEqual(1, self.get_output(self.widget.Outputs.selected_data).X.shape[1])

        self.widget.selection = {(0, 0), (0, 1)}
        self.widget.commit()
        self.assertEqual(50, len(self.get_output(self.widget.Outputs.selected_data)))
        self.assertEqual(2, self.get_output(self.widget.Outputs.selected_data).X.shape[1])

        self.widget.selection = {(0, 0), (1, 0)}
        self.widget.commit()
        self.assertEqual(100, len(self.get_output(self.widget.Outputs.selected_data)))
        self.assertEqual(1, self.get_output(self.widget.Outputs.selected_data).X.shape[1])


if __name__ == "__main__":
    unittest.main()

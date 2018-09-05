import unittest

from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owmarkergenes import OWMarkerGenes


class TestOWMarkerGenes(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMarkerGenes)

    def test_set_source(self):
        w = self.widget

        # reset combobox
        w.group_cb.clear()

        no_meta_data = Table.from_list(Domain([ContinuousVariable('test')]),
                                       [[1], [2], [3], [4]])
        w.set_source(no_meta_data)

        # data with no meta attributes does not set combobox values
        self.assertTrue(w.group_cb.count() == 0)

    def test_source_not_empty(self):
        w = self.widget
        self.assertIsNotNone(w.source)

    def test_row_selection(self):
        w = self.widget

        w.selected_group = w.group_cb.itemText(w.group_index)
        w.selected_genes = [('9246', 'Dendritic cell', 'Butler et al. (2018)'), ('9246', 'CD4', 'Butler et al. (2018)')]
        w.set_selection()

        out = self.get_output(self.widget.Outputs.genes)

        self.assertIsNotNone(out)
        self.assertIsInstance(out, Table)

        # check if rows were selected
        self.assertTrue(len(out) == len(w.selected_genes))

        # change group
        w.set_group_index(w.group_index + 1)
        self.assertNotEqual(w.group_cb.itemText(w.group_index - 1), w.selected_group)

    def test_read_data(self):
        self.widget.URL = "https://docs.google.com/spreadsheets/test_url"
        self.widget.read_data()
        self.assertNotEqual(self.widget.Warning.local_data.is_shown(),
                            self.widget.Error.file_not_found.is_shown())
        self.assertEqual(self.widget.Warning.local_data.is_shown(),
                         self.widget.source is not None)


if __name__ == '__main__':
    unittest.main()

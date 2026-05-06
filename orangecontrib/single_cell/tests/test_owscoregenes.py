import unittest

from AnyQt.QtCore import Qt
from AnyQt.QtWidgets import QStackedWidget, QCheckBox

from orangecontrib.single_cell.widgets.owscoregenes import OWRank
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWScoreGenes(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWRank)

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        attrs = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]
        domain = Domain(attributes=attrs, class_vars=class_var, metas=[])
        self.data = Table.from_numpy(
            domain,
            [[1, 2, 3, 4, 5], [4, 4, 4, 4, 4], [2, 3, 1, 1, 1], [-1, 0, 1, 0, 0]],
            [[0], [0], [0], [1]],
        )
    def test_data_no_attributes(self):
        """
        Do not crash when data has not attributes. Show an error.
        GH-51
        """
        w = self.widget
        self.assertFalse(w.Error.no_attributes.is_shown())
        # move attributes to metas
        data = self.data.transform(Domain([], self.data.domain.class_var, self.data.domain.attributes))
        for is_shown, input_data in ((True, data), (False, None)):
            self.send_signal(w.Inputs.data, input_data)
            self.assertEqual(is_shown, w.Error.no_attributes.is_shown())

    def test_store_settings(self):
        w = self.widget
        self.send_signal(w.Inputs.data, self.data)
        w.nSelected = 2
        w.setSelectionMethod(OWRank.SelectNBest)
        self.assertIsInstance(w.selected_attrs, list)
        self.assertEqual(len(w.selected_attrs), w.nSelected)
        self.assertIsInstance(w.selected_attrs[0], ContinuousVariable)
        w.setSelectionMethod(OWRank.SelectManual)
        settings = w.settingsHandler.pack_data(w)
        w1 = self.create_widget(OWRank, stored_settings=settings)
        self.send_signal(w1.Inputs.data, self.data)
        self.assertEqual(w1.selectionMethod, OWRank.SelectManual)
        self.assertSequenceEqual(w.selected_attrs, w1.selected_attrs)

    def test_update_selected_scores(self):
        w = self.widget
        self.send_signal(w.Inputs.data, self.data)
        stack = w.findChild(QStackedWidget)
        box = stack.widget(2)  # for unsupervised
        cb: QCheckBox = box.findChild(QCheckBox, "Mean")
        cb.setCheckState(Qt.CheckState.Unchecked)
        self.assertNotIn("Mean", w.selected_methods)
        cb.setCheckState(Qt.CheckState.Checked)
        self.assertIn("Mean", w.selected_methods)


if __name__ == '__main__':
    unittest.main()

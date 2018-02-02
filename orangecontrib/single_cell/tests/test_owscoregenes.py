import unittest

from orangecontrib.single_cell.widgets.owscoregenes import OWRank
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWScoreGenes(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWRank)

    def test_data_no_attributes(self):
        """
        Do not crash when data has not attributes. Show an error.
        GH-51
        """
        w = self.widget

        class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        metas = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]
        domain = Domain(attributes=[], class_vars=class_var, metas=metas)

        data = Table(
            domain,
            [['STG1', 1, 2, 3, 4, 5],
             ['STG1', 4, 4, 4, 4, 4],
             ['STG1', 2, 3, 1, 1, 1],
             ['STG2', -1, 0, 1, 0, 0]])

        self.assertFalse(w.Error.no_attributes.is_shown())
        for is_shown, input_data in ((True, data), (False, None)):
            self.send_signal(w.Inputs.data, input_data)
            self.assertEqual(is_shown, w.Error.no_attributes.is_shown())


if __name__ == '__main__':
    unittest.main()

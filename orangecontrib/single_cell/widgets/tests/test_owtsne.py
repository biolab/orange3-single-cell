import unittest

from orangecontrib.single_cell.widgets.owtsne import OWtSNE
from Orange.widgets.tests.base import WidgetTest
from Orange.data import DiscreteVariable, ContinuousVariable, Domain, Table


class TestOWtSNE(WidgetTest):

    def setUp(self):
        self.widget = self.create_widget(OWtSNE)

        self.class_var = DiscreteVariable('Stage name', values=['STG1', 'STG2'])
        self.attributes = [ContinuousVariable('GeneName' + str(i)) for i in range(5)]
        self.domain = Domain(self.attributes, class_vars=self.class_var)

        self.data = None

    def test_wrong_input(self):
        # no data
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNone(self.widget.data)

        # no rows
        self.data = Table(self.domain, [[1, 2, 3, 4, 5, 'STG1']])
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertIsNone(self.widget.data)

    def test_input(self):
        self.data = Table(self.domain, [[1, 1, 1, 1, 1, 'STG1'],
                                        [2, 2, 2, 2, 2, 'STG1'],
                                        [4, 4, 4, 4, 4, 'STG2'],
                                        [5, 5, 5, 5, 5, 'STG2']])

        self.send_signal(self.widget.Inputs.data, self.data)


if __name__ == '__main__':
    unittest.main()

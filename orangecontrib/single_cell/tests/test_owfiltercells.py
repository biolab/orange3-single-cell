import numpy as np

import Orange.data
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owfiltercells import OWFilterCells


class TestOWFilterCells(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWFilterCells)

    def tearDown(self):
        self.widget.onDeleteWidget()
        self.widget.deleteLater()
        del self.widget

    def test(self):
        domain = Orange.data.Domain(
            [Orange.data.ContinuousVariable("A1"),
             Orange.data.ContinuousVariable("A2"),
             Orange.data.ContinuousVariable("A3")]
        )
        data = Orange.data.Table.from_list(
            domain,
            [[0, 0, 0],
             [0, 0, 1],
             [0, 1, 1],
             [1, 1, 1]]
        )
        self.send_signal(self.widget.Inputs.data, data)

        out = self.get_output(self.widget.Outputs.data)
        self.assertTrue(len(out) == len(data))

        self.widget.limit_lower = 2
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        self.assertTrue(len(out), 2)
        np.testing.assert_array_equal(
            np.count_nonzero(out.X, axis=1),
            [2, 3]
        )
        self.send_signal(self.widget.Inputs.data, None)
        self.assertIsNone(self.get_output(self.widget.Outputs.data))


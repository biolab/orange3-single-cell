import numpy as np
import numpy.testing

import Orange.data
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owfilter import OWFilter, log1p, expm1


class TestOWFilterCells(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWFilter)

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
             [0, 1, 2],
             [1, 2, 3]]
        )
        self.send_signal(self.widget.Inputs.data, data)

        out = self.get_output(self.widget.Outputs.data)
        self.assertTrue(len(out) == len(data))

        self.widget.limit_lower = 2
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        self.assertEqual(len(out), 2)
        np.testing.assert_array_equal(
            np.count_nonzero(out.X, axis=1),
            [2, 3]
        )
        self.widget.set_filter_type(OWFilter.Genes)
        self.widget.limit_lower = 2
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            np.count_nonzero(out.X, axis=0),
            [2, 3]
        )
        self.widget.set_filter_type(OWFilter.Data)
        self.widget.limit_lower = 2
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            out.X,
            [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 2],
             [0, 2, 3]]
        )

        self.widget.limit_upper = 2
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            out.X,
            [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 2],
             [0, 2, 0]]
        )
        self.widget.set_upper_limit_enabled(False)
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            out.X,
            [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 2],
             [0, 2, 3]]
        )
        self.widget.set_lower_limit_enabled(False)
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            out.X,
            [[0, 0, 0],
             [0, 0, 1],
             [0, 1, 2],
             [1, 2, 3]]
        )
        self.widget.set_upper_limit_enabled(True)
        self.widget.commit()
        out = self.get_output(self.widget.Outputs.data)
        np.testing.assert_array_equal(
            out.X,
            [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 2],
             [0, 2, 3]]
        )
        #  For paint code coverage
        self.widget.grab()

        self.send_signal(self.widget.Inputs.data, None)
        self.assertIsNone(self.get_output(self.widget.Outputs.data))

    def test_utils(self):
        # test log1p and its inverse
        x = 10 ** np.arange(9)
        np.testing.assert_allclose(log1p(x), np.log10(1 + x))
        np.testing.assert_allclose(expm1(log1p(x)), x)
        x = 1 / (10 ** np.arange(1, 7))
        np.testing.assert_allclose(log1p(x), np.log10(1 + x))
        np.testing.assert_allclose(expm1(log1p(x)), x)


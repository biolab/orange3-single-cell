import unittest

from orangecontrib.single_cell.widgets.owscdatasets import OWscDataSets
from Orange.widgets.tests.base import WidgetTest


class TestOWscDataSets(WidgetTest):

    def test_widget_setup(self):
        self.widget = self.create_widget(OWscDataSets)


if __name__ == '__main__':
    unittest.main()

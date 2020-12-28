import unittest
from unittest.mock import patch

from orangecontrib.single_cell.widgets.owscdatasets import OWscDataSets
from Orange.widgets.tests.base import WidgetTest


class TestOWscDataSets(WidgetTest):

    def test_widget_setup(self):
        with patch.object(OWscDataSets, "list_remote", staticmethod(lambda: [])):
            self.widget = self.create_widget(OWscDataSets)
            self.wait_until_finished(self.widget)


if __name__ == '__main__':
    unittest.main()

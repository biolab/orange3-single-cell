import unittest
from unittest.mock import patch

import pkg_resources
from packaging import version

from orangecontrib.single_cell.widgets.owscdatasets import OWscDataSets
from Orange.widgets.tests.base import WidgetTest


class TestOWscDataSets(WidgetTest):
    def test_widget_setup(self):
        orange_version = pkg_resources.get_distribution("orange3").version
        if version.parse(orange_version) > version.parse("3.27.0"):
            with patch("Orange.widgets.data.owdatasets.list_remote", lambda _: {}):
                self.widget = self.create_widget(OWscDataSets)
                self.wait_until_finished(self.widget)
        else:
            with patch.object(OWscDataSets, "list_remote", staticmethod(lambda: [])):
                self.widget = self.create_widget(OWscDataSets)
                self.wait_until_finished(self.widget)


if __name__ == '__main__':
    unittest.main()

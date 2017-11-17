import sys

from AnyQt.QtWidgets import QApplication

import Orange.widgets.data.owdatasets


class OWscDataSets(Orange.widgets.data.owdatasets.OWDataSets):
    name = "Single Cell Datasets"
    description = "Load a data set from an online repository"
    icon = "icons/DataSets.svg"
    priority = 20

    INDEX_URL = "http://datasets.orange.biolab.si/sc/"
    DATASET_DIR = "sc-datasets"


def main(args=None):
    if args is None:
        args = sys.argv

    app = QApplication(list(args))
    w = OWscDataSets()
    w.show()
    w.raise_()
    rv = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return rv


if __name__ == "__main__":
    sys.exit(main())

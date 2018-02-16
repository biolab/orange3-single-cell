import sys

from AnyQt.QtWidgets import QApplication

from Orange.widgets.data import owdatasets


def hide_attributes(data):
    """Sets hidden=True for all attributes in large datasets"""
    if data is not None and len(data.domain.attributes) > 1000:
        for att in data.domain.attributes:
            att.attributes["hidden"] = True
    return data


class OWscDataSets(owdatasets.OWDataSets):
    name = "Single Cell Datasets"
    description = "Load a data set from an online repository"
    icon = "icons/SingleCellDatasets.svg"
    priority = 150

    INDEX_URL = "http://datasets.orange.biolab.si/sc/"
    DATASET_DIR = "sc-datasets"

    @staticmethod
    def load_data(path):
        data = owdatasets.OWDataSets.load_data(path)
        return hide_attributes(data)


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

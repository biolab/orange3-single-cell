import os
import sys
import Orange.widgets.data.owdatasets

from AnyQt.QtWidgets import QApplication
from AnyQt.QtGui import QStandardItemModel, QStandardItem
from AnyQt.QtCore import Qt

from Orange.widgets.gui import IndicatorItemDelegate

from orangecontrib.bioinformatics.ncbi.taxonomy import shortname, common_taxids
from orangecontrib.bioinformatics.ncbi.gene import ENTREZ_ID
from orangecontrib.bioinformatics.widgets.utils.data import GENE_AS_ATTRIBUTE_NAME, TAX_ID, GENE_ID_ATTRIBUTE


class OWscDataSets(Orange.widgets.data.owdatasets.OWDataSets):
    name = "Single Cell Datasets"
    description = "Load a data set from an online repository"
    icon = "icons/SingleCellDatasets.svg"
    priority = 150

    INDEX_URL = "https://datasets.biolab.si/sc/"
    DATASET_DIR = "sc-datasets"

    HEADER_SCHEMA = [
        ['islocal',      {'label': ''}],
        ['title',        {'label': 'Title'}],
        ['size',         {'label': 'Size'}],
        ['instances',    {'label': 'Cells'}],
        ['num_of_genes', {'label': 'Genes'}],
        ['taxid',        {'label': 'Organism'}],
        ['target',       {'label': 'Target'}],
        ['tags',         {'label': 'Tags'}]
    ]

    def load_data(self, path):
        info = self.selected_dataset()
        data = Orange.data.Table(path)
        data.attributes[TAX_ID] = info.taxid
        data.attributes[GENE_AS_ATTRIBUTE_NAME] = True  # Will all data sets have gene names in columns?
        data.attributes[GENE_ID_ATTRIBUTE] = ENTREZ_ID
        return data

    def assign_delegates(self):
        self.view.setItemDelegateForColumn(
            self.Header.islocal, IndicatorItemDelegate(self, role=Qt.DisplayRole)
        )
        self.view.setItemDelegateForColumn(
            self.Header.size,
            Orange.widgets.data.owdatasets.SizeDelegate(self))

        self.view.setItemDelegateForColumn(
            self.Header.instances,
            Orange.widgets.data.owdatasets.NumericalDelegate(self)
        )
        self.view.setItemDelegateForColumn(
            self.Header.num_of_genes,
            Orange.widgets.data.owdatasets.NumericalDelegate(self)
        )

    def create_model(self):
        allkeys = set(self.allinfo_local)

        if self.allinfo_remote is not None:
            allkeys = allkeys | set(self.allinfo_remote)

        allkeys = sorted(allkeys)

        model = QStandardItemModel(self)
        model.setHorizontalHeaderLabels(self._header_labels)

        current_index = -1
        for i, file_path in enumerate(allkeys):
            data_info = self._parse_info(file_path)
            row = []

            for info_tag, header_setting in self.HEADER_SCHEMA:
                item = QStandardItem()

                try:
                    data = data_info.__getattribute__(info_tag)
                except AttributeError:
                    # unknown tag in JSON
                    data = ''

                # first column indicating cached data sets
                if info_tag == 'islocal':
                    item.setData(' ' if data else '', Qt.DisplayRole)
                    item.setData(data_info, Qt.UserRole)

                else:
                    # parse taxid to common name
                    if info_tag == 'taxid' and data in common_taxids():
                        data = shortname(data)[0].title()

                    if info_tag == 'tags':
                        data = ', '.join(data) if data else ''

                    item.setData(data, Qt.DisplayRole)

                # set icon to Target column
                if info_tag == 'target' and data:
                    item.setIcon(Orange.widgets.data.owdatasets.variable_icon(data))

                row.append(item)
            model.appendRow(row)

            if os.path.join(*file_path) == self.selected_id:
                current_index = i

        return model, current_index


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


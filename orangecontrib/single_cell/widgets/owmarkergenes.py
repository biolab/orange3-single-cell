import sys
import os


from typing import List
from enum import EnumMeta
from collections import defaultdict
from functools import partial


from AnyQt.QtCore import (
    Qt, QSize, QSortFilterProxyModel, QModelIndex,
    QItemSelection, QItemSelectionModel, QItemSelectionRange
)
from AnyQt.QtGui import QFont, QColor
from AnyQt.QtWidgets import QTreeView, QLineEdit

from Orange.data import MISSING_VALUES, Table
from Orange.data.io import UrlReader
from Orange.misc.environ import data_dir
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import TableModel


from orangecontrib.bioinformatics.widgets.utils.data import (
    GENE_AS_ATTRIBUTE_NAME, TAX_ID, GENE_ID_COLUMN
)


def local_cache_path(path):
    return os.path.join(data_dir(), path)


class FilterProxyModel(QSortFilterProxyModel):
    def filterAcceptsRow(self, source_row, source_parent):
        # type: (int, QModelIndex) -> bool
        model = self.sourceModel()
        if isinstance(model, LinkedTableModel):
            table = model.source
            domain = table.domain
            index = model.mapToSourceRows(source_row)
            row = table[index]
            regexp = self.filterRegExp()
            return any(regexp.indexIn(text) != -1 for text in
                       (str(row[var]) if var.is_string else var.str_val(row[var])
                        for var in domain.variables + domain.metas
                        if var.is_string or var.is_discrete))
        else:
            return True


class HeaderIndex(EnumMeta):
    NAME = 0
    GENE = 1
    CELL_TYPE = 2
    FUNCTION = 3
    REFERENCE = 4
    URL = 5


HeaderLabels = {HeaderIndex.GENE: 'Entrez ID',
                HeaderIndex.REFERENCE: 'Reference',
                HeaderIndex.URL: 'URL'}


class LinkedTableModel(TableModel):

    ClassVar, Meta, Attribute = range(3)

    ColorForRole = {
        ClassVar: None,
        Meta: None,
        Attribute: None,
    }

    def __init__(self, data, parent):
        TableModel.__init__(self, data[:, data.domain.metas[:-1]], parent)
        self._data = data
        self._roleData = {Qt.DisplayRole: self.source}
        self._roleData = partial(
            defaultdict,
            partial(defaultdict,
                    partial(defaultdict, lambda: None)))(self._roleData)
        self.set_column_links()

    def set_column_links(self):
        domain = self._data.domain
        ref_col = domain.metas.index(domain[HeaderLabels[HeaderIndex.REFERENCE]])
        font = QFont()
        font.setUnderline(True)
        color = QColor(Qt.blue)
        for i, row in enumerate(self._data):
            link = row[HeaderLabels[HeaderIndex.URL]].value
            if len(link):
                self._roleData[gui.LinkRole][i][ref_col] = link
                self._roleData[Qt.FontRole][i][ref_col] = font
                self._roleData[Qt.ForegroundRole][i][ref_col] = color

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            cell_data = super().data(index, role)
            return "" if cell_data in MISSING_VALUES else cell_data
        elif role in (gui.LinkRole, Qt.FontRole, Qt.ForegroundRole):
            row, col = index.row(), index.column()
            return self._roleData[role][row][col]
        return super().data(index, role)


class MarkerGroupContextHandler(settings.ContextHandler):

    def __init__(self):
        super().__init__()

    def match(self, context, group, *args):
        if not context.group == group:
            return self.NO_MATCH

        return self.PERFECT_MATCH

    def new_context(self, group):
        context = super().new_context()
        context.group = group
        return context


class OWMarkerGenes(widget.OWWidget):
    name = "Marker Genes"
    icon = 'icons/MarkerGenes.svg'
    priority = 170

    URL = "https://docs.google.com/spreadsheets/d/1ik-Ju5F-" \
          "wcsFhjszM7pRZXTTfSJl_EQOn3uNIzYfAY4/edit#gid=0"
    DIR_NAME = "datasets/sc/"
    FILE_NAME = DIR_NAME + "markers.tab"

    class Outputs:
        genes = widget.Output("Genes", Table)

    want_main_area = False

    selected_group = settings.Setting("")  # type: str
    filter_text = settings.Setting('')     # type: str
    header_state = settings.Setting(b'')   # type: bytes

    settingsHandler = MarkerGroupContextHandler()
    selected_genes = settings.ContextSetting([])  # type: List[tuple]

    class Error(widget.OWWidget.Error):
        file_not_found = widget.Msg("File not found.")

    class Warning(widget.OWWidget.Warning):
        local_data = widget.Msg("File not found, local cached data is shown.")

    def __init__(self):
        super().__init__()
        self.source = None
        self.group_index = -1
        self.group_cb = gui.comboBox(self.controlArea, self, "group_index")
        self.group_cb.activated[int].connect(self.set_group_index)

        # TODO: to avoid this, marker genes table should have 'tax_id' column
        self.map_group_to_taxid = {
            'Human': '9606',
            'Mouse': '10090'
        }

        filter_line_edit = gui.lineEdit(
            self.controlArea, self, "filter_text"
        )  # type: QLineEdit
        filter_line_edit.setPlaceholderText("Filter...")
        filter_line_edit.textEdited.connect(self.set_filter_str)

        self.view = view = QTreeView(
            rootIsDecorated=False,
            uniformRowHeights=True,
            selectionMode=QTreeView.ExtendedSelection,
            sortingEnabled=True,
        )

        self.proxy_model = FilterProxyModel(
            self, filterCaseSensitivity=Qt.CaseInsensitive,
        )
        view.setModel(self.proxy_model)
        view.selectionModel().selectionChanged.connect(
            self._on_selection_changed
        )
        view.viewport().setMouseTracking(True)
        self.controlArea.layout().addWidget(view)

        self.read_data()
        if self.header_state:
            view.header().restoreState(self.header_state)

    def read_data(self):
        try:
            data = UrlReader(self.URL).read()
        except Exception:
            data = self._read_cached_data()
            if data is None:
                return
        else:
            dir_path = local_cache_path(self.DIR_NAME)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
            data.save(local_cache_path(self.FILE_NAME))
        self.set_source(data)
        self.commit()

    def _read_cached_data(self):
        data = None
        try:
            data = Table(local_cache_path(self.FILE_NAME))
        except OSError:
            self.Error.file_not_found()
        else:
            self.Warning.local_data()
        return data

    def set_selection(self):
        selected = self.selected_rows()

        if len(selected):
            header_count = self.view.header().count() - 1

            if self.view.model().rowCount() <= selected[-1]:
                return

            selection = QItemSelection()

            for row_index in selected:
                selection.append(
                    QItemSelectionRange(
                        self.view.model().index(row_index, 0),
                        self.view.model().index(row_index, header_count)
                    )
                )

            self.view.selectionModel().select(
                selection, QItemSelectionModel.ClearAndSelect)

    def set_source(self, data):
        # type: (Table) -> None
        """
        Set the source data from which to fetch the output

        The output is a subset filtered on the first meta column (group)
        """
        self.source = data
        domain = data.domain

        if domain.metas:
            group = domain.metas[0]
            groupcol, _ = data.get_column_view(group)

            if group.is_string:
                group_values = list(map(str, unique(groupcol)))
            elif group.is_discrete:
                group_values = group.values
            else:
                raise TypeError("Invalid column type")
            try:
                idx = group_values.index(self.selected_group)
            except ValueError:
                idx = -1

            self.group_cb.clear()
            self.group_cb.addItems(group_values)
            if idx != -1:
                self.group_index = idx
                self.selected_group = group_values[idx]
            elif group_values:
                self.group_index = min(max(self.group_index, 0),
                                       len(group_values) - 1)

            self.set_group_index(self.group_index)
            self._setup()

    def set_group_index(self, group_index):
        self.closeContext()
        self.group_index = group_index
        self.selected_group = self.group_cb.itemText(group_index)
        self._setup()

    def set_filter_str(self, string):
        if string != self.filter_text:
            self.filter_text = string
        proxy = self.view.model()
        assert isinstance(proxy, QSortFilterProxyModel)
        proxy.setFilterFixedString(string)

    def _setup(self):
        self.closeContext()
        data = self.source
        group = data.domain.metas[0]
        gvec = data.get_column_view(group)[0]
        if group.is_string:
            mask = gvec == self.group_cb.itemData(self.group_index,
                                                  Qt.DisplayRole)
        else:
            mask = gvec == self.group_index

        data = data[mask]
        rest = data[:, data.domain.metas[1:]]
        model = LinkedTableModel(rest, parent=self)
        ref_col = rest.domain.metas.index(rest.domain[HeaderLabels[HeaderIndex.REFERENCE]])
        self.view.setItemDelegateForColumn(
            ref_col, gui.LinkStyledItemDelegate(self.view))

        if self.proxy_model.sourceModel():
            self.proxy_model.sourceModel().deleteLater()
        self.proxy_model.setSourceModel(model)

        self.openContext(self.selected_group)
        self.set_filter_str(self.filter_text)
        self.set_selection()

        self.commit()

    def _on_selection_changed(self):
        self.commit()

    def selected_rows(self):
        """ Return row index for selected genes
        """
        if not self.selected_genes:
            return []

        model = self.view.model()
        return [row_index for row_index in range(model.rowCount())
                if (model.index(row_index, HeaderIndex.GENE).data(),
                    model.index(row_index, HeaderIndex.CELL_TYPE).data(),
                    model.index(row_index, HeaderIndex.REFERENCE).data())
                in self.selected_genes]

    def commit(self):
        model = self.view.model()
        assert isinstance(model, QSortFilterProxyModel)
        table = model.sourceModel()
        assert isinstance(table, LinkedTableModel)
        rows = [model.mapToSource(mi).row()
                for mi in self.view.selectionModel().selectedRows(0)]

        if rows:
            rows = table.mapToSourceRows(rows)
            output = table.source[rows]
        else:
            output = table.source

        gene_id = self.view.selectionModel().selectedRows(HeaderIndex.GENE)
        cell_type = self.view.selectionModel().selectedRows(HeaderIndex.CELL_TYPE)
        ref = self.view.selectionModel().selectedRows(HeaderIndex.REFERENCE)

        self.selected_genes = [(entrez.data(), cell.data(), ref.data())
                               for entrez, cell, ref in zip(gene_id, cell_type, ref)]

        # always false for marker genes data tables in single cell
        output.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        # set taxonomy id in data.attributes
        output.attributes[TAX_ID] = self.map_group_to_taxid.get(self.selected_group, '')
        # set column id flag
        output.attributes[GENE_ID_COLUMN] = HeaderLabels[HeaderIndex.GENE]
        output.name = "Marker Genes"

        self.Outputs.genes.send(output)

    def closeEvent(self, event):
        self.header_state = bytes(self.view.header().saveState())
        super().closeEvent(event)

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(600, 500))


def unique(iterable):
    seen = set()

    def observed(el):
        observed = el in seen
        seen.add(el)
        return observed

    return (el for el in iterable if not observed(el))


def main(argv=None):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(argv or sys.argv)
    w = OWMarkerGenes()
    w.show()
    w.activateWindow()
    rv = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return rv


if __name__ == "__main__":
    sys.exit(main())

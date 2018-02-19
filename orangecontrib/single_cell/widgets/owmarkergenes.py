import sys
import os

from AnyQt.QtCore import Qt, QSize, QSortFilterProxyModel, QModelIndex
from AnyQt.QtWidgets import QTreeView, QLineEdit

import Orange.data

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import TableModel


def resource_path(path):
    return os.path.join(os.path.dirname(__file__), path)


class FilterProxyModel(QSortFilterProxyModel):
    def filterAcceptsRow(self, source_row, source_parent):
        # type: (int, QModelIndex) -> bool
        model = self.sourceModel()
        if isinstance(model, TableModel):
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


class OWMarkerGenes(widget.OWWidget):
    name = "Marker Genes"
    icon = 'icons/MarkerGenes.svg'
    priority = 158

    class Outputs:
        genes = widget.Output("Genes", Orange.data.Table)

    want_main_area = False

    selected_group = settings.Setting("")  # type: str
    header_state = settings.Setting(b'')   # type: bytes

    def __init__(self):
        super().__init__()
        self.source = None
        self.group_index = -1
        self.filter_text = ""
        self.group_cb = gui.comboBox(self.controlArea, self, "group_index")
        self.group_cb.activated[int].connect(self.set_group_index)

        filter = gui.lineEdit(
            self.controlArea, self, "filter_text"
        )  # type: QLineEdit
        filter.setPlaceholderText("Filter...")
        filter.textEdited.connect(self.set_filter_str)
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
        self.controlArea.layout().addWidget(view)

        # NEEDS updating/
        self.set_source(Orange.data.Table(resource_path("data/markers.tab")))
        if self.header_state:
            view.header().restoreState(self.header_state)
        self.commit()

    def set_source(self, data):
        # type: (Orange.data.Table) -> None
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
            self._setup()

    def set_group_index(self, group_index):
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
        model = TableModel(rest, parent=self)
        if self.proxy_model.sourceModel():
            self.proxy_model.sourceModel().deleteLater()
        self.proxy_model.setSourceModel(model)
        self.commit()

    def _on_selection_changed(self):
        self.commit()

    def commit(self):
        model = self.view.model()
        assert isinstance(model, QSortFilterProxyModel)
        table = model.sourceModel()
        assert isinstance(table, TableModel)
        rows = [model.mapToSource(mi).row()
                for mi in self.view.selectionModel().selectedRows(0)]
        if rows:
            rows = table.mapToSourceRows(rows)
            output = table.source[rows]
        else:
            output = table.source
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

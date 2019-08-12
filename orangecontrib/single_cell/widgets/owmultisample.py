import os
import numbers
from collections import namedtuple
from typing import Dict, Tuple, List, Optional
from serverfiles import sizeformat

from AnyQt.QtCore import (
    Qt, QItemSelectionModel, pyqtSlot as Slot, QModelIndex, Signal
)
from AnyQt.QtGui import (
    QStandardItemModel, QStandardItem, QIcon, QBrush, QColor
)
from AnyQt.QtWidgets import (
    QStyledItemDelegate, QTreeView, QSizePolicy,
    QStyle, QHeaderView, QLineEdit, QFileDialog
)
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Msg

from orangecontrib.single_cell.widgets.load_data import Concatenate, Loader
import orangecontrib.single_cell.widgets.owloaddata as owloaddata

LoaderObjectRole = next(gui.OrangeUserRole)


class FileDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, str):
            option.text = os.path.split(data)[1]


class SizeDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, numbers.Integral):
            option.text = sizeformat(int(data))
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter


class NumericalDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, numbers.Number):
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter


class SparsityDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, numbers.Number):
            option.text = "~" + str(round(data * 100)) + " %"
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter


class RecentLineEdit(QLineEdit):
    def __init__(self, parent):
        super().__init__(parent)
        self.setText("Recent files")
        self.setReadOnly(True)

    def mouseReleaseEvent(self, event):
        self.parent().showPopup()


class MultiSampleModel(QStandardItemModel):
    def supportedDropActions(self):
        return Qt.MoveAction

    def supportedDragActions(self):
        return Qt.MoveAction

    def dropMimeData(self, data, action, row, column, parent):
        column = 0
        if parent.row() != -1:
            row = parent.row()
        return super().dropMimeData(data, action, row, column, QModelIndex())


class MultiSampleTreeView(QTreeView):
    drop_finished = Signal()
    data_changed = Signal()

    def __init__(self, parent):
        self._parent = parent
        super().__init__(
            sortingEnabled=True,
            selectionMode=QTreeView.SingleSelection,
            selectionBehavior=QTreeView.SelectRows,
            alternatingRowColors=True
        )
        self.setRootIsDecorated(False)
        self.setDragEnabled(True)
        self.setDragDropMode(QTreeView.InternalMove)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

        else:
            super().dragEnterEvent(event)

    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                event.acceptProposedAction()
                if os.path.isfile(url.toLocalFile()):
                    self._parent.set_current_path(url.toLocalFile())
            self._parent.select_last_item()
            self.drop_finished.emit()
        else:
            super().dropEvent(event)

    def startDrag(self, action):
        super().startDrag(action)
        self.drop_finished.emit()

    def commitData(self, event):
        super().commitData(event)
        self.data_changed.emit()


class OWMultiSample(owloaddata.OWLoadData):
    name = "Load Data"
    description = "Load samples for multi-sample analysis"
    icon = "icons/LoadData.svg"
    priority = 110

    class Information(owloaddata.OWLoadData.Information):
        file_in_list = Msg("File {} already in the list.")

    HEADER_SCHEMA = (
        ("selected", ""),
        ("name", "File"),
        ("source", "Source Name"),
        ("cells", "Cells"),
        ("genes", "Genes"),
        ("size", "Size"),
        ("sparsity", "Sparsity"),
        ("remove", "")
    )
    OUTPUT_TYPES = ["Intersection", "Union"]

    # override settings
    _recent_row_annotations = []
    _recent_col_annotations = []
    _cells_in_rows = False
    _col_annotations_enabled = False
    _row_annotations_enabled = False
    _last_path = ""
    _header_rows_count = 1
    _header_cols_count = 1
    _sample_rows_enabled = False
    _sample_cols_enabled = False
    _sample_cols_p = 10.0
    _sample_rows_p = 10.0

    samples = Setting([])  # type: List[Tuple[str, str, bool]
    loaders = Setting({})  # type: Dict[str, Loader]
    output_type = Setting(Concatenate.INTERSECTION)
    auto_commit = Setting(False)

    resizing_enabled = True

    def __init__(self):

        self._header_labels = hls = [label for _, label in self.HEADER_SCHEMA]
        header = namedtuple("header", [tag for tag, _ in self.HEADER_SCHEMA])
        self._Header = header(*[index for index, _ in enumerate(hls)])

        self.view = MultiSampleTreeView(self)
        self._view_setup()
        self._assign_delegates()

        super().__init__()
        file_label = self.file_layout.itemAtPosition(0, 0).widget()
        file_label.setText("Add file:")

        self.recent_combo.setLineEdit(RecentLineEdit(self))
        self.recent_combo.setCurrentText("Recent files")
        self.__recent_combo_view_resize()

        open_button = self.file_layout.itemAtPosition(0, 2).widget()
        open_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        open_button.setText("Browse...")
        open_button.setIcon(QIcon())

        self.summary_label.hide()
        self.load_data_button.hide()
        self.controlArea.layout().removeItem(self.file_layout)
        self.controlArea.setSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.MinimumExpanding)

        main_box = gui.hBox(self.left_side, margin=4)
        left_box = gui.vBox(main_box, spacing=10)
        left_box.layout().addLayout(self.file_layout)
        left_box.layout().addWidget(self.view)
        right_box = gui.vBox(main_box)
        right_box.layout().addWidget(self.controlArea)
        output_box = gui.hBox(left_box)
        box = gui.widgetBox(output_box, True)
        gui.radioButtons(
            box, self, "output_type", self.OUTPUT_TYPES,
            label="Output genes:", orientation=Qt.Horizontal,
            callback=self._output_type_changed
        )
        gui.rubber(output_box)
        cb_box = gui.auto_commit(output_box, self, "auto_commit", "Load Data")
        cb_box.setMaximumWidth(250)
        self.left_side.layout().addWidget(main_box)
        self.read_settings()
        self.view.setFocus()

        self.cells_in_rows_changed.connect(self._update_view_cells_genes)

    def _update_view_cells_genes(self):
        row = self.view.currentIndex().row()
        if row > -1:
            cell_item = self.view.model().item(row, self._Header.cells)
            cell_item.setData(self._data_loader.n_cells, Qt.DisplayRole)
            gene_item = self.view.model().item(row, self._Header.genes)
            gene_item.setData(self._data_loader.n_genes, Qt.DisplayRole)
            self.repaint()

    def _view_clicked(self, index):
        if index.column() == self._Header.remove:
            self.remove_item(index)
        elif index.column() == self._Header.selected:
            self.commit()

    def _selection_changed(self):
        rows = self.view.selectionModel().selectedRows(0)
        current = rows[0] if rows else None  # type: Optional[QModelIndex]
        if current is not None:
            self.set_current_loader(
                current.data(LoaderObjectRole),
                current.sibling(
                    current.row(), self._Header.name
                ).data(Qt.DisplayRole)
            )

    def _output_type_changed(self):
        self.commit()

    def _assign_delegates(self):
        self.view.setItemDelegateForColumn(
            self._Header.name, FileDelegate(self)
        )
        self.view.setItemDelegateForColumn(
            self._Header.cells, NumericalDelegate(self)
        )
        self.view.setItemDelegateForColumn(
            self._Header.genes, NumericalDelegate(self)
        )
        self.view.setItemDelegateForColumn(
            self._Header.size, SizeDelegate(self)
        )
        self.view.setItemDelegateForColumn(
            self._Header.sparsity, SparsityDelegate(self)
        )

    def _view_setup(self):
        model = MultiSampleModel(self)
        model.setHorizontalHeaderLabels(self._header_labels)
        self.view.setModel(model)
        self.view.clicked[QModelIndex].connect(self._view_clicked)
        self.view.data_changed.connect(lambda: self.commit())
        self.view.drop_finished.connect(lambda: self.commit())
        self.view.header().setStretchLastSection(False)
        self.view.header().setSectionResizeMode(
            self._Header.name, QHeaderView.Stretch
        )
        self.view.selectionModel().selectionChanged.connect(
            self._selection_changed
        )
        self.view.resizeColumnToContents(self._Header.selected)
        self.view.setColumnWidth(self._Header.remove, 20)
        self._resize_columns_to_contents()

    def _resize_columns_to_contents(self):
        self.view.resizeColumnToContents(self._Header.cells)
        self.view.resizeColumnToContents(self._Header.genes)
        self.view.resizeColumnToContents(self._Header.size)
        self.view.resizeColumnToContents(self._Header.sparsity)

    def add_item(self, checked, source_name):
        path = self.current_path()
        if not path:
            return

        self.Information.file_in_list.clear()
        if len(self.view.model().findItems(path, Qt.MatchFixedString, 1)):
            self.Information.file_in_list(path)
            return

        exists = os.path.exists(path)
        self.view.model().appendRow([
            self.__selected_item(checked, exists),
            self.__file_item(exists),
            self.__source_item(source_name),
            self.__cells_item(),
            self.__genes_item(),
            self.__size_item(),
            self.__sparsity_item(),
            self.__remove_item()
        ])
        self._resize_columns_to_contents()
        self.commit()

    def __selected_item(self, checked, exists):
        item = QStandardItem()
        item.setData(self._data_loader, LoaderObjectRole)
        item.setCheckable(exists)
        item.setCheckState(Qt.Checked if checked and exists else Qt.Unchecked)
        item.setEditable(False)
        return item

    def __file_item(self, exists):
        item = QStandardItem()
        item.setData(self.current_path(), Qt.DisplayRole)
        item.setToolTip(self.current_path())
        item.setEditable(False)
        if not exists:
            item.setToolTip("Missing from file system")
            item.setForeground(QBrush(QColor(Qt.red)))
        return item

    def __source_item(self, source):
        return QStandardItem(source or os.path.split(self.current_path())[1])

    def __cells_item(self):
        item = QStandardItem()
        item.setData(self._data_loader.n_cells, Qt.DisplayRole)
        item.setEditable(False)
        return item

    def __genes_item(self):
        item = QStandardItem()
        item.setData(self._data_loader.n_genes, Qt.DisplayRole)
        item.setEditable(False)
        return item

    def __size_item(self):
        item = QStandardItem()
        item.setData(self._data_loader.file_size, Qt.DisplayRole)
        item.setEditable(False)
        return item

    def __sparsity_item(self):
        item = QStandardItem()
        item.setData(self._data_loader.sparsity, Qt.DisplayRole)
        item.setEditable(False)
        return item

    def __remove_item(self):
        item = QStandardItem()
        item.setIcon(
            self.style().standardIcon(QStyle.SP_DockWidgetCloseButton))
        item.setEditable(False)
        return item

    def select_item(self, index):
        if not isinstance(index, QModelIndex):
            index = self.view.model().index(index, 0)
        self.view.selectionModel().setCurrentIndex(
            index,
            QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        self.view.setFocus()
        self.repaint()

    def remove_item(self, index):
        self.view.model().removeRows(index.row(), 1)
        self.commit()
        if not self.view.model().rowCount():
            self.set_current_loader(Loader(), "")

    def select_last_item(self):
        n_rows = self.view.model().rowCount()
        if n_rows == 0:
            return
        n_cols = self.view.model().columnCount()
        self.select_item(self.view.model().index(n_rows - 1, n_cols - 1))

    def set_current_loader(self, loader, path):
        self._data_loader = loader
        self._current_path = path
        self.setup_gui()

    def set_current_path(self, path, checked=True, source_name=""):
        super().set_current_path(path)
        self.add_item(checked, source_name)
        self.recent_combo.setCurrentText("Recent files")

    def _select_recent(self, index):
        super()._select_recent(index)
        self.select_last_item()

    @Slot()
    def browse(self):
        dlg = QFileDialog(
            self, acceptMode=QFileDialog.AcceptOpen,
            fileMode=QFileDialog.ExistingFiles,
        )
        filters = owloaddata.Formats
        dlg.setNameFilters(filters)
        if filters:
            dlg.selectNameFilter(filters[-1])
        if dlg.exec_() == QFileDialog.Accepted:
            for filename in dlg.selectedFiles():
                self.set_current_path(filename)
                self.select_last_item()
        self.__recent_combo_view_resize()

    def __recent_combo_view_resize(self):
        self.recent_combo.view().setMinimumHeight(
            self.recent_combo.view().sizeHintForRow(0) *
            self.recent_model.rowCount())

    def _invalidate(self):
        self.commit()

    def commit(self):
        if not self.view or not self.view.model():
            return

        self.clear_messages()
        self.write_settings()

        data_collection = []
        for index, (path, checked, source_name) in enumerate(self.samples):
            if checked:
                loader = self.loaders.get(path)
                table = loader()
                if any(loader.errors.values()):
                    self.select_item(index)
                    self.show_error_messages()
                    break
                else:
                    data_collection.append((table, source_name))

        data = Concatenate.concatenate(self.output_type, data_collection)
        self.Outputs.data.send(data)

    def write_settings(self):
        self.samples = []
        self.loaders = {}
        model = self.view.model()
        for i in range(model.rowCount()):
            item = model.item(i)
            checked, loader = item.checkState(), item.data(LoaderObjectRole)
            path = model.item(i, self._Header.name).text()
            source_name = model.item(i, self._Header.source).text()
            self.samples.append((path, checked, source_name))
            loader.recent_path = self.relocate_path(loader.recent_path)
            loader.row_annotation_file = self.relocate_path(
                loader.row_annotation_file)
            loader.col_annotation_file = self.relocate_path(
                loader.col_annotation_file)
            self.loaders[path] = loader

    def saveSettings(self):
        self.write_settings()
        super().saveSettings()

    def _saveState(self):
        super()._saveState()
        self.write_settings()

    def read_settings(self):
        samples = self.samples.copy()
        loaders = self.loaders.copy()
        for path, checked, source_name in samples:
            self._current_path = path
            loader = loaders.get(path)
            if loader is not None:
                loader = loader.copy()
                resolved_path = self.resolve_path(loader.recent_path)
                loader.recent_path = resolved_path
                loader.row_annotation_file = self.resolve_path(
                    loader.row_annotation_file)
                loader.col_annotation_file = self.resolve_path(
                    loader.col_annotation_file)
                self._current_path = resolved_path.abspath
                self.set_current_loader(loader, resolved_path.abspath)
                self.add_item(checked, source_name)
        self.select_last_item()
        self.recent_combo.setCurrentText("Recent files")

    def send_report(self):
        model = self.view.model()
        paths = [os.path.split(model.item(i, self._Header.name).text())[1] for
                 i in range(model.rowCount()) if model.item(i).checkState()]
        self.report_items(
            "Concatenated samples",
            [("Files", ", ".join(paths) if paths else "No files."),
             ("Output genes", self.OUTPUT_TYPES[self.output_type])])


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication

    app = QApplication([])
    w = OWMultiSample()
    w.show()
    w.raise_()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()

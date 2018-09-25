import os
import sys

from typing import List

from serverfiles import sizeformat

from AnyQt.QtCore import Qt, QFileInfo, QTimer, Signal
from AnyQt.QtGui import QStandardItemModel, QStandardItem
from AnyQt.QtWidgets import (
    QSizePolicy, QGridLayout, QHBoxLayout, QFormLayout,
    QLabel, QComboBox, QSpinBox, QCheckBox, QPushButton,
    QStyle, QApplication, QFileDialog, QFileIconProvider,
    QWidget
)
from AnyQt.QtCore import pyqtSlot as Slot

from Orange.data import Table
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.filedialogs import RecentPath
from Orange.widgets.utils.buttons import VariableTextPushButton

from orangecontrib.single_cell.widgets.load_data import get_data_loader, Loader


Formats = [
    "Count file (*.count)",
    "Tab separated file (*.tsv *.tab)",
    "Comma separated file (*.csv)",
    "10x gene-barcode matrix (matrix.mtx)",
    "Microsoft Excel spreadsheet (*.xls *.xlsx)",
    "Pickled Python object file (*.pkl *.pickle)",
    "All files (*.*)"
]

AnnotationFormats = [
    "Meta file (*.meta)",
    "Tab separated file (*.tsv *.tab)",
    "Comma separated file (*.csv)",
    "All files (*.*)",
]


def RecentPath_asqstandarditem(pathitem):
    # type: (RecentPath) -> QStandardItem
    icon_provider = QFileIconProvider()
    # basename of a normalized name (strip right path component separators)
    basename = os.path.basename(os.path.normpath(pathitem.abspath))
    item = QStandardItem(
        icon_provider.icon(QFileInfo(pathitem.abspath)),
        basename
    )
    item.setToolTip(pathitem.abspath)
    item.setData(pathitem, Qt.UserRole)
    return item


def init_recent_paths_model(model, paths, relpaths=[]):
    # type: (QStandardItemModel, List[RecentPath]) -> None
    for pathitem in paths:
        item = RecentPath_asqstandarditem(pathitem)
        abspath = pathitem.search(searchpaths=relpaths)
        if not (abspath and os.path.exists(abspath)):
            item.setEnabled(False)
            item.setSelectable(False)
            item.setToolTip(item.toolTip() + " (Missing from file system)")
        model.appendRow(item)


def insert_recent_path(model, path, relpaths=[]):
    # type: (QStandardItemModel, RecentPath) -> int
    index = -1
    for i in range(model.rowCount()):
        item = model.item(i, 0)
        pathitem = item.data(Qt.UserRole)
        if isinstance(pathitem, RecentPath) and \
                samepath(pathitem.abspath, path.abspath):
            index = i
            break
    if index != -1:
        item = model.takeRow(index)
    else:
        item = RecentPath_asqstandarditem(path)
    model.insertRow(0, item)
    return 0


def samepath(p1, p2):
    return (os.path.normpath(os.path.normcase(p1)) ==
            os.path.normpath(os.path.normcase(p2)))


class RunaroundSettingsHandler(settings.SettingsHandler):
    def pack_data(self, widget):
        widget._saveState()
        return super().pack_data(widget)


class OWLoadData(widget.OWWidget):
    name = ""
    icon = "icons/LoadData.svg"
    priority = 110

    class Outputs:
        data = widget.Output("Data", Table)

    class Information(widget.OWWidget.Information):
        modified = widget.Msg(
            "Uncommited changes\nPress 'Load data' to submit changes"
        )

    class Warning(widget.OWWidget.Warning):
        sampling_in_effect = widget.Msg("Sampling is in effect.")

    class Error(widget.OWWidget.Error):
        row_annotation_mismatch = widget.Msg(
            "Row annotation length mismatch\n"
            "Expected {} rows got {}"
        )
        col_annotation_mismatch = widget.Msg(
            "Column annotation length mismatch\n"
            "Expected {} rows got {}"
        )
        inadequate_headers = widget.Msg(
            "Headers and Row Labels error"
        )
        reading_error = widget.Msg(
            "Cannot read data using given parameters."
        )

    _recent = settings.Setting([])  # type: List[RecentPath]
    _recent_row_annotations = settings.Setting([])  # type: List[RecentPath]
    _recent_col_annotations = settings.Setting([])  # type: List[RecentPath]
    _cells_in_rows = settings.Setting(False)
    _col_annotations_enabled = settings.Setting(False)
    _row_annotations_enabled = settings.Setting(False)

    _last_path = settings.Setting("")  # type: str
    _header_rows_count = settings.Setting(1)  # type: int
    _header_cols_count = settings.Setting(1)  # type: int
    _sample_rows_enabled = settings.Setting(False)  # type: bool
    _sample_cols_enabled = settings.Setting(False)  # type: bool
    _sample_cols_p = settings.Setting(10.0)  # type: bool
    _sample_rows_p = settings.Setting(10.0)  # type: bool

    settingsHandler = RunaroundSettingsHandler()

    want_main_area = False
    resizing_enabled = False

    cells_in_rows_changed = Signal()

    def __init__(self):
        super().__init__()
        self._current_path = ""
        self._data_loader = Loader()
        icon_open_dir = self.style().standardIcon(QStyle.SP_DirOpenIcon)

        # Top grid with file selection combo box
        self.file_layout = grid = QGridLayout()
        lb = QLabel("File:")
        lb.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.recent_combo = cb = QComboBox(
            sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLengthWithIcon,
            minimumContentsLength=20,
            toolTip="Select a recent file"
        )
        self.recent_model = cb.model()  # type: QStandardItemModel
        self.recent_combo.activated[int].connect(self._select_recent)

        browse = QPushButton("...", autoDefault=False, icon=icon_open_dir,
                             clicked=self.browse)

        # reload = QPushButton("Reload", autoDefault=False, icon=icon_reload)

        grid.addWidget(lb, 0, 0, Qt.AlignVCenter)
        grid.addWidget(cb, 0, 1)
        grid.addWidget(browse, 0, 2)
        # grid.addWidget(reload, 0, 3)

        self.summary_label = label = QLabel("", self)
        label.ensurePolished()
        f = label.font()
        if f.pointSizeF() != -1:
            f.setPointSizeF(f.pointSizeF() * 5 / 6)
        else:
            f.setPixelSize(f.pixelSize() * 5 / 6)
        label.setFont(f)
        grid.addWidget(label, 1, 1, 1, 3)

        self.controlArea.layout().addLayout(grid)

        box = gui.widgetBox(
            self.controlArea, "Headers and Row Labels", spacing=-1
        )
        hl = QHBoxLayout()
        hl.setContentsMargins(0, 0, 0, 0)
        self.header_rows_spin = spin = QSpinBox(
            box, minimum=0, maximum=3, value=self._header_rows_count,
            keyboardTracking=False
        )
        spin.valueChanged.connect(self.set_header_rows_count)
        hl.addWidget(QLabel("Data starts with", box))
        hl.addWidget(self.header_rows_spin)
        hl.addWidget(QLabel("header row(s)", box))
        hl.addStretch(10)
        box.layout().addLayout(hl)

        hl = QHBoxLayout()
        hl.setContentsMargins(0, 0, 0, 0)
        self.header_cols_spin = spin = QSpinBox(
            box, minimum=0, maximum=3, value=self._header_cols_count,
            keyboardTracking=False
        )
        spin.valueChanged.connect(self.set_header_cols_count)

        hl.addWidget(QLabel("First", box))
        hl.addWidget(self.header_cols_spin)
        hl.addWidget(QLabel("column(s) are row labels", box))
        hl.addStretch(10)
        box.layout().addLayout(hl)

        self.data_struct_box = box = gui.widgetBox(
            self.controlArea, "Input Data Structure"
        )
        gui.radioButtons(
            box, self, "_cells_in_rows",
            ["Genes in rows, cells in columns",
             "Cells in rows, genes in columns"],
            callback=self._cells_in_rows_changed
        )

        box = gui.widgetBox(
            self.controlArea, "Sample Data", spacing=-1)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        box.layout().addLayout(grid)

        self.sample_rows_cb = cb = QCheckBox(checked=self._sample_rows_enabled)
        self.sample_rows_p_spin = spin = QSpinBox(
            minimum=0, maximum=100, value=self._sample_rows_p
        )
        spin.valueChanged.connect(self.set_sample_rows_p)
        suffix = QLabel("% of cells")
        cb.toggled.connect(self.set_sample_rows_enabled)

        grid.addWidget(cb, 0, 0)
        grid.addWidget(spin, 0, 1)
        grid.addWidget(suffix, 0, 2)

        self.sample_cols_cb = cb = QCheckBox(checked=self._sample_cols_enabled)
        self.sample_cols_p_spin = spin = QSpinBox(
            minimum=0, maximum=100, value=self._sample_cols_p
        )
        spin.valueChanged.connect(self.set_sample_cols_p)
        suffix = QLabel("% of genes")
        cb.toggled.connect(self.set_sample_cols_enabled)

        grid.addWidget(cb, 1, 0)
        grid.addWidget(spin, 1, 1)
        grid.addWidget(suffix, 1, 2)
        grid.setColumnStretch(3, 10)

        self.annotation_files_box = box = gui.widgetBox(
            self.controlArea, "Cell && Gene Annotation Files"
        )
        form = QFormLayout(
            formAlignment=Qt.AlignLeft,
            rowWrapPolicy=QFormLayout.WrapAllRows,
        )
        box.layout().addLayout(form)

        self.row_annotations_cb = cb = QCheckBox(
            "Cell annotations", checked=self._row_annotations_enabled
        )
        self._row_annotations_w = w = QWidget(enabled=self._row_annotations_enabled)
        cb.toggled.connect(self.set_row_annotations_enabled)
        cb.toggled.connect(w.setEnabled)
        hl = QHBoxLayout()
        hl.setContentsMargins(0, 0, 0, 0)
        w.setLayout(hl)
        self.row_annotations_combo = QComboBox(
            sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLengthWithIcon,
            minimumContentsLength=18
        )
        self.row_annotations_combo.activated.connect(
            self._row_annotations_combo_changed
        )
        hl.addWidget(self.row_annotations_combo)
        hl.addWidget(QPushButton("...", box, autoDefault=False,
                                 icon=icon_open_dir,
                                 clicked=self.browse_row_annotations))
        # hl.addWidget(QPushButton("Reload", box, autoDefault=False,
        #                          icon=icon_reload))
        form.addRow(cb, w)

        self.col_annotations_cb = cb = QCheckBox(
            "Gene annotations", checked=self._col_annotations_enabled
        )
        self._col_annotations_w = w = QWidget(enabled=self._col_annotations_enabled)
        cb.toggled.connect(self.set_col_annotations_enabled)
        cb.toggled.connect(w.setEnabled)
        hl = QHBoxLayout()
        hl.setContentsMargins(0, 0, 0, 0)
        w.setLayout(hl)
        self.col_annotations_combo = QComboBox(
            sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLengthWithIcon,
            minimumContentsLength=18
        )
        self.col_annotations_combo.activated.connect(
            self._col_annotations_combo_changed
        )
        hl.addWidget(self.col_annotations_combo)
        hl.addWidget(QPushButton("...", box, autoDefault=False,
                                 icon=icon_open_dir,
                                 clicked=self.browse_col_annotations))
        # hl.addWidget(QPushButton("Reload", box, autoDefault=False,
        #                          icon=icon_reload))
        form.addRow(cb, w)

        self.controlArea.layout().addStretch(10)
        self.load_data_button = button = VariableTextPushButton(
            "Load data", autoDefault=True, textChoiceList=["Load data", "Reload"]
        )
        self.load_data_button.setAutoDefault(True)
        button.clicked.connect(self.commit, Qt.QueuedConnection)
        self.controlArea.layout().addWidget(button, alignment=Qt.AlignRight)

        init_recent_paths_model(
            self.recent_model,
            [self.resolve_path(p) for p in self._recent],
        )
        init_recent_paths_model(
            self.row_annotations_combo.model(),
            [self.resolve_path(p) for p in self._recent_row_annotations]
        )
        init_recent_paths_model(
            self.col_annotations_combo.model(),
            [self.resolve_path(p) for p in self._recent_col_annotations]
        )
        self._update_summary()
        self._update_warning()

        if self._last_path != "" and os.path.exists(self._last_path):
            QTimer.singleShot(
                0, lambda: self.set_current_path(self._last_path)
            )
        else:
            self.recent_combo.setCurrentIndex(-1)

    def resolve_path(self, path):
        basedir = self.workflowEnv().get("basedir", None)
        if not basedir or not path:
            return path
        return path.resolve([("basedir", basedir)]) or path

    def _cells_in_rows_changed(self):
        self._data_loader.transposed = not self._cells_in_rows
        self._invalidate()
        self.cells_in_rows_changed.emit()

    def _row_annotations_combo_changed(self):
        path = self.row_annotations_combo.currentData(Qt.UserRole)
        if isinstance(path, RecentPath) and os.path.exists(path.abspath):
            self._data_loader.row_annotation_file = path  # type: RecentPath
        else:
            self._data_loader.row_annotation_file = None
        self._invalidate()

    def _col_annotations_combo_changed(self):
        path = self.col_annotations_combo.currentData(Qt.UserRole)
        if isinstance(path, RecentPath) and os.path.exists(path.abspath):
            self._data_loader.col_annotation_file = path  # type: RecentPath
        else:
            self._data_loader.col_annotation_file = None
        self._invalidate()

    def _update_warning(self):
        if (self._sample_rows_enabled and self._sample_rows_p < 100) or \
                (self._sample_cols_enabled and self._sample_cols_p < 100):
            self.Warning.sampling_in_effect()
        else:
            self.Warning.sampling_in_effect.clear()

    def set_sample_rows_enabled(self, enabled, commit=True):
        if self._sample_rows_enabled != enabled:
            self._sample_rows_enabled = enabled
            self.sample_rows_cb.setChecked(enabled)
            self._update_warning()
            self._data_loader.sample_rows_enabled = enabled
            if commit:
                self._invalidate()

    def set_sample_cols_enabled(self, enabled, commit=True):
        if self._sample_cols_enabled != enabled:
            self._sample_cols_enabled = enabled
            self.sample_cols_cb.setChecked(enabled)
            self._update_warning()
            self._data_loader.sample_cols_enabled = enabled
            if commit:
                self._invalidate()

    def set_sample_rows_p(self, p, commit=True):
        if self._sample_rows_p != p:
            self._sample_rows_p = p
            self._update_warning()
            self.sample_rows_p_spin.setValue(p)
            self._data_loader.sample_rows_p = p
            if commit:
                self._invalidate()

    def set_sample_cols_p(self, p, commit=True):
        if self._sample_cols_p != p:
            self._sample_cols_p = p
            self._update_warning()
            self.sample_cols_p_spin.setValue(p)
            self._data_loader.sample_cols_p = p
            if commit:
                self._invalidate()

    def set_header_rows_count(self, n, commit=True):
        if self._header_rows_count != n:
            self._header_rows_count = n
            self.header_rows_spin.setValue(n)
            self._data_loader.header_rows_count = n
            if commit:
                self._invalidate()

    def set_header_cols_count(self, n, commit=True):
        if self._header_cols_count != n:
            self._header_cols_count = n
            self.header_cols_spin.setValue(n)
            self._data_loader.header_cols_count = n
            if commit:
                self._invalidate()

    def set_row_annotations_enabled(self, enabled, commit=True):
        if self._row_annotations_enabled != enabled:
            self._row_annotations_enabled = enabled
            self.row_annotations_cb.setChecked(enabled)
            self._data_loader.row_annotations_enabled = enabled
            if commit:
                self._invalidate()

    def set_col_annotations_enabled(self, enabled, commit=True):
        if self._col_annotations_enabled != enabled:
            self._col_annotations_enabled = enabled
            self.col_annotations_cb.setChecked(enabled)
            self._data_loader.col_annotations_enabled = enabled
            if commit:
                self._invalidate()

    def set_current_path(self, path):
        if samepath(self._current_path, path):
            return

        model = self.recent_model
        index = -1
        pathitem = None
        for i in range(model.rowCount()):
            item = model.item(i)
            data = item.data(Qt.UserRole) if item is not None else None
            if isinstance(data, RecentPath) and samepath(path, data.abspath):
                index, pathitem = i, data
                break

        rpaths = []

        if pathitem is None:
            assert index == -1
            pathitem = RecentPath.create(path, rpaths)

        if index != -1:
            item = model.takeRow(index)
        else:
            item = RecentPath_asqstandarditem(pathitem)

        model.insertRow(0, item)
        self._current_path = path
        self.recent_combo.setCurrentIndex(0)

        self._data_loader = get_data_loader(path)
        self._update_summary()
        self.setup_gui()
        self._invalidate()

    def setup_gui(self):
        """ Use loader predefined values. If the value is None, set
        loader's parameter to widget's setting value.
        """
        loader = self._data_loader
        if loader.header_rows_count is not None:
            self.set_header_rows_count(loader.header_rows_count, False)
        else:
            loader.header_rows_count = self._header_rows_count

        if loader.header_cols_count is not None:
            self.set_header_cols_count(loader.header_cols_count, False)
        else:
            loader.header_cols_count = self._header_cols_count

        if loader.transposed is not None:
            self._cells_in_rows = not loader.transposed
        else:
            loader.transposed = not self._cells_in_rows

        if loader.sample_rows_enabled is not None:
            self.set_sample_rows_enabled(loader.sample_rows_enabled, False)
        else:
            loader.sample_rows_enabled = self._sample_rows_enabled

        if loader.sample_cols_enabled is not None:
            self.set_sample_cols_enabled(loader.sample_cols_enabled, False)
        else:
            loader.sample_cols_enabled = self._sample_cols_enabled

        if loader.sample_rows_p is not None:
            self.set_sample_rows_p(loader.sample_rows_p, False)
        else:
            loader.sample_rows_p = self._sample_rows_p

        if loader.sample_cols_p is not None:
            self.set_sample_cols_p(loader.sample_cols_p, False)
        else:
            loader.sample_cols_p = self._sample_cols_p

        if loader.row_annotation_file is not None:
            index = insert_recent_path(
                self.row_annotations_combo.model(),
                self.resolve_path(loader.row_annotation_file)
            )
            self.row_annotations_combo.setCurrentIndex(index)
            self.set_row_annotations_enabled(
                loader.row_annotations_enabled, False
            )
        else:
            self.row_annotations_combo.setCurrentIndex(-1)
            self.set_row_annotations_enabled(False, False)

        if loader.col_annotation_file is not None:
            index = insert_recent_path(
                self.col_annotations_combo.model(),
                self.resolve_path(loader.col_annotation_file)
            )
            self.col_annotations_combo.setCurrentIndex(index)
            self.set_col_annotations_enabled(
                loader.col_annotations_enabled, False)
        else:
            self.col_annotations_combo.setCurrentIndex(-1)
            self.set_col_annotations_enabled(False, False)

        self.header_rows_spin.setEnabled(loader.FIXED_FORMAT)
        self.header_cols_spin.setEnabled(loader.FIXED_FORMAT)
        self.data_struct_box.setEnabled(loader.FIXED_FORMAT)

        self.annotation_files_box.setEnabled(loader.ENABLE_ANNOTATIONS)

    def _update_summary(self):
        size = self._data_loader.file_size
        ncols = self._data_loader.n_cols
        nrows = self._data_loader.n_rows
        text = []
        if size is not None:
            text += [sizeformat(size)]
        if nrows is not None:
            text += ["{:n} rows".format(nrows)]
        if nrows is not None:
            text += ["{:n} columns".format(ncols)]

        self.summary_label.setText(", ".join(text))

    def current_path(self):
        return self._current_path

    def _select_recent(self, index):
        # type: (int) -> None
        # select a file from the recent list (entered via combo box `activate`)
        assert 0 <= index < self.recent_model.rowCount()
        item = self.recent_model.item(index)
        pathitem = item.data(Qt.UserRole)
        assert isinstance(pathitem, RecentPath)
        self.set_current_path(pathitem.abspath)

    @Slot()
    def browse(self):
        dlg = QFileDialog(self)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setFileMode(QFileDialog.ExistingFile)

        filters = Formats
        dlg.setNameFilters(filters)
        if filters:
            dlg.selectNameFilter(filters[0])
        if dlg.exec_() == QFileDialog.Accepted:
            filename = dlg.selectedFiles()[0]
            self.set_current_path(filename)

    @Slot()
    def browse_row_annotations(self):
        dlg = QFileDialog(
            self, acceptMode=QFileDialog.AcceptOpen,
            fileMode=QFileDialog.ExistingFile
        )
        filters = AnnotationFormats
        dlg.setNameFilters(filters)

        if filters:
            dlg.selectNameFilter(filters[0])
        if dlg.exec_() == QFileDialog.Accepted:
            filename = dlg.selectedFiles()[0]
            m = self.row_annotations_combo.model()  # type: QStandardItemModel
            pathitem = RecentPath.create(filename, [])
            index = insert_recent_path(m, pathitem)
            self.row_annotations_combo.setCurrentIndex(index)
            self._invalidate()

    @Slot()
    def browse_col_annotations(self):
        dlg = QFileDialog(
            self, acceptMode=QFileDialog.AcceptOpen,
            fileMode=QFileDialog.ExistingFile
        )
        filters = AnnotationFormats
        dlg.setNameFilters(filters)

        if filters:
            dlg.selectNameFilter(filters[0])
        if dlg.exec_() == QFileDialog.Accepted:
            filename = dlg.selectedFiles()[0]
            m = self.col_annotations_combo.model()  # type: QStandardItemModel
            pathitem = RecentPath.create(filename, [])
            index = insert_recent_path(m, pathitem)
            self.col_annotations_combo.setCurrentIndex(index)
            self._invalidate()

    def _invalidate(self):
        self.set_modified(True)

    def set_modified(self, modified):
        if modified:
            text = "Load data"
        else:
            text = "Reload"

        self.load_data_button.setText(text)
        self.load_data_button.setAutoDefault(modified)
        # Setting autoDefault once also sets default which persists even after
        # settings autoDefault back to False??
        self.load_data_button.setDefault(modified)
        self.Information.modified(shown=modified)

    def commit(self):
        path = self._current_path
        if not path:
            return
        self.Outputs.data.send(self._data_loader())
        self.show_error_messages()
        self.set_modified(False)

    def show_error_messages(self):
        self.Error.row_annotation_mismatch.clear()
        self.Error.col_annotation_mismatch.clear()
        self.Error.inadequate_headers.clear()
        errors = self._data_loader.errors
        if len(errors["row_annot_mismatch"]):
            self.Error.row_annotation_mismatch(*errors["row_annot_mismatch"])
        if len(errors["col_annot_mismatch"]):
            self.Error.col_annotation_mismatch(*errors["col_annot_mismatch"])
        if len(errors["inadequate_headers"]):
            self.Error.inadequate_headers(*errors["inadequate_headers"])
        if len(errors["reading_error"]):
            self.Error.reading_error()

    def onDeleteWidget(self):
        super().onDeleteWidget()

    def _saveState(self):
        maxitems = 15

        def dataiter(model, role=Qt.UserRole):
            return (model.data(model.index(i, 0), role)
                    for i in range(model.rowCount()))

        def recent_paths(model):
            return [self.relocate_path(el) for el in dataiter(model)
                    if isinstance(el, RecentPath)][:maxitems]

        self._recent = recent_paths(self.recent_model)
        self._recent_row_annotations = recent_paths(
            self.row_annotations_combo.model())
        self._recent_col_annotations = recent_paths(
            self.col_annotations_combo.model())
        self._last_path = self._current_path

    def relocate_path(self, path):
        basedir = self.workflowEnv().get("basedir", None)
        if not basedir or not path:
            return path
        return RecentPath.create(path.abspath, [("basedir", basedir)])

    def saveSettings(self):
        self._saveState()
        super().saveSettings()


def main(argv=None):
    app = QApplication(argv or [])
    w = OWLoadData()
    w.show()
    w.raise_()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

import os
import sys
import csv
from itertools import chain
from typing import List
from types import SimpleNamespace

from serverfiles import sizeformat

import scipy.io
import scipy.sparse
import numpy as np

import pandas as pd

from AnyQt.QtCore import Qt, QFileInfo
from AnyQt.QtGui import QStandardItemModel, QStandardItem
from AnyQt.QtWidgets import (
    QSizePolicy, QGridLayout, QHBoxLayout, QFormLayout,
    QLabel, QComboBox, QSpinBox, QCheckBox, QPushButton,
    QStyle, QApplication, QFileDialog, QFileIconProvider,
    QWidget
)
from AnyQt.QtCore import pyqtSlot as Slot

import Orange.data

from Orange.data import ContinuousVariable, StringVariable

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.filedialogs import RecentPath
from Orange.widgets.utils.buttons import VariableTextPushButton


class Options(SimpleNamespace):
    format = None

    header_rows_count = None
    header_cols_count = None

    transposed = None

    sample_rows_p = None
    sample_cols_p = None

    # supplementary column/rows annotation files
    row_annotation_file = None
    column_annotation_file = None


def infer_options(path):
    dirname, basename = os.path.split(path)
    basename_no_ext, ext = os.path.splitext(basename)

    options = Options()
    if ext == ".mtx":
        genes_path = os.path.join(dirname, "genes.tsv")
        if os.path.isfile(genes_path):
            options.column_annotation_file = genes_path
        barcodes_path = os.path.join(dirname, "barcodes.tsv")
        if os.path.isfile(barcodes_path):
            options.row_annotation_file = barcodes_path
        options.transposed = True
    elif ext == ".count":
        meta_path = os.path.join(dirname, basename_no_ext + ".meta")
        if os.path.isfile(meta_path):
            options.row_annotation_file = meta_path
        options.transposed = True

    return options


Formats = [
    "Count file (*.count)",
    "Tab separated file (*.tsv)",
    "Comma separated file (*.csv)",
    "10x results - Matrix Market Exchange format (*.mtx)",
    "Any tab separated file (*.*)"
]

AnnotationFormats = [
    "Meta file (*.meta)",
    "Tab separated file (*.tsv)",
    "Comma separated file (*.csv)",
    "Any tab separated file (*.*)",
]


def separator_from_filename(path):
    path, ext = os.path.splitext(path)
    if ext == ".csv":
        return ","
    elif ext == ".tsv":
        return "\t"
    elif ext == ".count" or ext == ".meta":
        return "\t"
    else:  # assume tab separated
        return "\t"


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
    name = "Load Data"
    icon = "icons/LoadData.svg"

    class Outputs:
        data = widget.Output("Data", Orange.data.Table)

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

    _recent = settings.Setting([])  # type: List[str]
    _recent_row_annotations = settings.Setting([])  # type: List[str]
    _recent_col_annotations = settings.Setting([])  # type: List[str]
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

    def __init__(self):
        super().__init__()
        self._current_path = ""
        icon_open_dir = self.style().standardIcon(QStyle.SP_DirOpenIcon)

        # Top grid with file selection combo box
        grid = QGridLayout()
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
            callback=self._invalidate
        )

        box = gui.widgetBox(
            self.controlArea, "Sample Data", spacing=-1)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        box.layout().addLayout(grid)

        self.sample_rows_cb = cb = QCheckBox(checked=self._sample_rows_enabled)

        spin = QSpinBox(
            minimum=0, maximum=100, value=self._sample_rows_p,
            enabled=self._sample_rows_enabled
        )
        spin.valueChanged.connect(self.set_sample_rows_p)
        suffix = QLabel("% of cells", enabled=self._sample_rows_enabled)
        cb.toggled.connect(self.set_sample_rows_enabled)
        cb.toggled.connect(spin.setEnabled)
        cb.toggled.connect(suffix.setEnabled)

        grid.addWidget(cb, 0, 0)
        grid.addWidget(spin, 0, 1)
        grid.addWidget(suffix, 0, 2)

        self.sample_cols_cb = cb = QCheckBox(checked=self._sample_cols_enabled)
        spin = QSpinBox(
            minimum=0, maximum=100, value=self._sample_cols_p,
            enabled=self._sample_cols_enabled
        )
        spin.valueChanged.connect(self.set_sample_cols_p)
        suffix = QLabel("% of genes", enabled=self._sample_cols_enabled)
        cb.toggled.connect(self.set_sample_cols_enabled)
        cb.toggled.connect(spin.setEnabled)
        cb.toggled.connect(suffix.setEnabled)

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
        self.row_annotations_combo.activated.connect(self._invalidate)
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
        self.col_annotations_combo.activated.connect(self._invalidate)
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
            [RecentPath.create(p, []) for p in self._recent],
        )
        init_recent_paths_model(
            self.row_annotations_combo.model(),
            [RecentPath.create(p, []) for p in self._recent_row_annotations]
        )
        init_recent_paths_model(
            self.col_annotations_combo.model(),
            [RecentPath.create(p, []) for p in self._recent_col_annotations]
        )
        self._update_summary()
        self._update_warning()

        if self._last_path != "" and os.path.exists(self._last_path):
            self.set_current_path(self._last_path)
        else:
            self.recent_combo.setCurrentIndex(-1)

    def _update_warning(self):
        if (self._sample_rows_enabled and self._sample_rows_p < 100) or \
                (self._sample_cols_enabled and self._sample_cols_p < 100):
            self.Warning.sampling_in_effect()
        else:
            self.Warning.sampling_in_effect.clear()

    def set_sample_rows_enabled(self, enabled):
        if self._sample_rows_enabled != enabled:
            self._sample_rows_enabled = enabled
            self._update_warning()
            self._invalidate()

    def set_sample_cols_enabled(self, enabled):
        if self._sample_cols_enabled != enabled:
            self._sample_cols_enabled = enabled
            self._update_warning()
            self._invalidate()

    def set_sample_rows_p(self, p):
        if self._sample_rows_p != p:
            self._sample_rows_p = p
            self._update_warning()
            self._invalidate()

    def set_sample_cols_p(self, p):
        if self._sample_cols_p != p:
            self._sample_cols_p = p
            self._update_warning()
            self._invalidate()

    def set_header_rows_count(self, n):
        if self._header_rows_count != n:
            self._header_rows_count = n
            self.header_rows_spin.setValue(n)
            self._invalidate()

    def set_header_cols_count(self, n):
        if self._header_cols_count != n:
            self._header_cols_count = n
            self.header_cols_spin.setValue(n)
            self._invalidate()

    def set_row_annotations_enabled(self, enabled):
        if self._row_annotations_enabled != enabled:
            self._row_annotations_enabled = enabled
            self.row_annotations_cb.setChecked(enabled)
            self._invalidate()

    def set_col_annotations_enabled(self, enabled):
        if self._col_annotations_enabled != enabled:
            self._col_annotations_enabled = enabled
            self.col_annotations_cb.setChecked(enabled)
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

        opts = infer_options(path)

        if path.endswith(".count"):
            self.set_header_rows_count(1)
            self.set_header_cols_count(1)
            fixed_format = False
        elif path.endswith(".mtx"):
            self.set_header_rows_count(0)
            self.set_header_cols_count(0)
            fixed_format = False
            self._cells_in_rows = False
        else:
            fixed_format = True

        if opts.transposed is not None:
            self._cells_in_rows = not opts.transposed

        self.data_struct_box.setEnabled(fixed_format)
        self.header_rows_spin.setEnabled(fixed_format)
        self.header_cols_spin.setEnabled(fixed_format)

        model.insertRow(0, item)
        self._current_path = path
        self.recent_combo.setCurrentIndex(0)
        self._update_summary()

        if opts.row_annotation_file is not None:
            index = insert_recent_path(
                self.row_annotations_combo.model(),
                RecentPath.create(opts.row_annotation_file, [])
            )
            self.row_annotations_combo.setCurrentIndex(index)
            self.set_row_annotations_enabled(True)
        else:
            self.row_annotations_combo.setCurrentIndex(-1)

        if opts.column_annotation_file is not None:
            index = insert_recent_path(
                self.col_annotations_combo.model(),
                RecentPath.create(opts.column_annotation_file, [])
            )
            self.col_annotations_combo.setCurrentIndex(index)
            self.set_col_annotations_enabled(True)
        else:
            self.col_annotations_combo.setCurrentIndex(-1)

        self._invalidate()

    def _update_summary(self):
        path = self._current_path

        size = None
        ncols = None
        nrows = None

        try:
            st = os.stat(path)
        except OSError:
            pass
        else:
            size = st.st_size

        if os.path.splitext(path)[1] == ".mtx":
            try:
                with open(path, "rb") as f:
                    nrows, ncols = scipy.io.mminfo(f)[:2]
            except OSError:
                pass
            except ValueError:
                pass
        else:
            try:
                with open(path, "rt", encoding="latin-1") as f:
                    sep = separator_from_filename(path)
                    ncols = len(next(csv.reader(f, delimiter=sep)))
                    nrows = sum(1 for _ in f)
            except OSError:
                pass
            except StopIteration:
                pass

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
            f = dlg.selectedNameFilter()
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
            f = dlg.selectedNameFilter()
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
            f = dlg.selectedNameFilter()
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

        transpose = not self._cells_in_rows
        row_annot = self.row_annotations_combo.currentData(Qt.UserRole)
        col_annot = self.col_annotations_combo.currentData(Qt.UserRole)

        if self._row_annotations_enabled and \
                isinstance(row_annot, RecentPath) and \
                os.path.exists(row_annot.abspath):
            row_annot = row_annot.abspath  # type: str
        else:
            row_annot = None

        if self._col_annotations_enabled and \
                isinstance(col_annot, RecentPath) and \
                os.path.exists(col_annot.abspath):
            col_annot = col_annot.abspath  # type: str
        else:
            col_annot = None

        meta_parts = []  # type: List[pd.DataFrame]
        attrs = []  # type: List[ContinuousVariable]
        metas = []  # type: List[StringVariable]

        rstate = np.random.RandomState(0x667)

        skip_row = skip_col = None
        if self._sample_cols_enabled:
            p = self._sample_cols_p
            if p < 100:
                def skip_col(i, p=p):
                    return i > 3 and rstate.uniform(0, 100) > p

        if self._sample_rows_enabled:
            p = self._sample_rows_p
            if p < 100:
                def skip_row(i, p=p):
                    return i > 3 and rstate.uniform(0, 100) > p

        header_rows = self._header_rows_count
        header_rows_indices = []
        if header_rows == 0:
            header_rows = None
        elif header_rows == 1:
            header_rows = 0
            header_rows_indices = [0]
        else:
            header_rows = list(range(header_rows))
            header_rows_indices = header_rows

        header_cols = self._header_cols_count
        header_cols_indices = []
        if header_cols == 0:
            header_cols = None
        elif header_cols == 1:
            header_cols = 0
            header_cols_indices = [0]
        else:
            header_cols = list(range(header_cols))
            header_cols_indices = header_cols

        if transpose:
            _skip_row, _skip_col = skip_col, skip_row
        else:
            _skip_col, _skip_row = skip_col, skip_row

        _userows = _usecols = None
        userows_mask = usecols_mask = None

        if _skip_col is not None:
            ncols = pd.read_csv(
                path, sep=separator_from_filename(path), index_col=None,
                nrows=1).shape[1]
            usecols_mask = np.array([
                not _skip_col(i) or i in header_cols_indices
                for i in range(ncols)
                ], dtype=bool)
            _usecols = np.flatnonzero(usecols_mask)

        if _skip_row is not None:
            userows_mask = []  # record the used rows
            def _skip_row(i, test=_skip_row):
                r = test(i)
                userows_mask.append(r)
                return r

        meta_df_index = None
        row_annot_header = 0
        row_annot_columns = None
        col_annot_header = 0
        col_annot_columns = None

        if os.path.splitext(path)[1] == ".mtx":
            # 10x cellranger output
            X = scipy.io.mmread(path)
            assert isinstance(X, scipy.sparse.coo_matrix)
            if transpose:
                X = X.T
            if _skip_row is not None:
                userows_mask = np.array(
                    [not _skip_row(i) for i in range(X.shape[0])]
                )
                X = X.tocsr()[np.flatnonzero(userows_mask)]
            if _skip_col is not None:
                usecols_mask = np.array(
                    [not _skip_col(i) for i in range(X.shape[1])]
                )
                X = X.tocsc()[:, np.flatnonzero(usecols_mask)]
            X = X.todense(order="F")
            if userows_mask is not None:
                meta_df = pd.DataFrame({}, index=np.flatnonzero(userows_mask))
            else:
                meta_df = pd.DataFrame({}, index=pd.RangeIndex(X.shape[0]))

            meta_df_index = meta_df.index

            row_annot_header = None
            row_annot_columns = ["Barcodes"]
            col_annot_header = None
            col_annot_columns = ["Id", "Gene"]
            leading_cols = leading_rows = 0
        else:
            df = pd.read_csv(
                path, sep=separator_from_filename(path),
                index_col=header_cols, header=header_rows,
                skiprows=_skip_row, usecols=_usecols
            )

            if _skip_row is not None:
                userows_mask = np.array(userows_mask, dtype=bool)

            if transpose:
                df = df.transpose()
                userows_mask, usecols_mask = usecols_mask, userows_mask
                leading_rows = len(header_cols_indices)
                leading_cols = len(header_rows_indices)
            else:
                leading_rows = len(header_rows_indices)
                leading_cols = len(header_cols_indices)

            X = df.values
            attrs = [ContinuousVariable.make(str(g)) for g in df.columns]

            meta_df = df.iloc[:, :0]  # Take the index # type: pd.DataFrame
            meta_df_index = df.index
            meta_parts = (meta_df, )

        self.Error.row_annotation_mismatch.clear()
        self.Error.col_annotation_mismatch.clear()

        if row_annot is not None:
            row_annot_df = pd.read_csv(
                row_annot, sep=separator_from_filename(row_annot),
                header=row_annot_header, names=row_annot_columns,
                index_col=None
            )
            if userows_mask is not None:
                # NOTE: we account for column header/ row index
                expected = len(userows_mask) - leading_rows
            else:
                expected = X.shape[0]
            if len(row_annot_df) != expected:
                self.Error.row_annotation_mismatch(expected, len(row_annot_df))
                row_annot_df = None

            if row_annot_df is not None and userows_mask is not None:
                # use the same sample indices
                indices = np.flatnonzero(userows_mask[leading_rows:])
                row_annot_df = row_annot_df.iloc[indices]
                # if path.endswith(".count") and row_annot.endswith('.meta'):
                #     assert np.all(row_annot_df.iloc[:, 0] == df.index)

            if row_annot_df is not None and meta_df_index is not None:
                # Try to match the leading columns with the meta_df_index.
                # If found then drop the columns (or index if the level does
                # not have a name but the annotation col does)
                drop_cols = []
                drop_index_level = []
                for i in range(meta_df_index.nlevels):
                    meta_df_level = meta_df_index.get_level_values(i)
                    if np.all(row_annot_df.iloc[:, i] == meta_df_level):
                        if meta_df_level.name is None:
                            drop_index_level.append(i)
                        elif meta_df_level.name == row_annot_df.columns[i].name:
                            drop_cols.append(i)

                if drop_cols:
                    row_annot_df = row_annot_df.drop(columns=drop_cols)

                if drop_index_level:
                    for i in reversed(drop_index_level):
                        if isinstance(meta_df.index, pd.MultiIndex):
                            meta_df_index = meta_df_index.droplevel(i)
                        else:
                            assert i == 0
                            meta_df_index = pd.RangeIndex(meta_df_index.size)
                    meta_df = pd.DataFrame({}, index=meta_df_index)

            if row_annot_df is not None:
                meta_parts = (meta_df, row_annot_df)

        if col_annot is not None:
            col_annot_df = pd.read_csv(
                col_annot, sep=separator_from_filename(col_annot),
                header=col_annot_header, names=col_annot_columns,
                index_col=None
            )
            if usecols_mask is not None:
                expected = len(usecols_mask) - leading_cols
            else:
                expected = X.shape[1]
            if len(col_annot_df) != expected:
                self.Error.col_annotation_mismatch(expected, len(col_annot_df))
                col_annot_df = None
            if col_annot_df is not None and usecols_mask is not None:
                indices = np.flatnonzero(usecols_mask[leading_cols:])
                col_annot_df = col_annot_df.iloc[indices]

            if col_annot_df is not None:
                assert len(col_annot_df) == X.shape[1]
                if not attrs and X.shape[1]:  # No column names yet
                    attrs = [ContinuousVariable.make(str(v))
                             for v in col_annot_df.iloc[:, 0]]
                names = [str(c) for c in col_annot_df.columns]
                for var, values in zip(attrs, col_annot_df.values):
                    var.attributes.update({n: v for n, v in zip(names, values)})

        if meta_parts:
            meta_parts = [df_.reset_index() if not df_.index.is_integer()
                          else df_ for df_ in meta_parts]
            metas = [StringVariable.make(name)
                     for name in chain(*(_.columns for _ in meta_parts))]
            M = np.hstack(tuple(df_.values for df_ in meta_parts))
        else:
            metas = None
            M = None

        if not attrs and X.shape[1]:
            attrs = Orange.data.Domain.from_numpy(X).attributes

        domain = Orange.data.Domain(attrs, metas=metas)
        d = Orange.data.Table.from_numpy(domain, X, None, M)
        self.Outputs.data.send(d)

        self.set_modified(False)

    def onDeleteWidget(self):
        super().onDeleteWidget()

    def _saveState(self):
        maxitems = 15

        def dataiter(model, role=Qt.UserRole):
            return (model.data(model.index(i, 0), role)
                    for i in range(model.rowCount()))

        def recent_paths(model):
            return [el.abspath for el in dataiter(model)
                    if isinstance(el, RecentPath)][:maxitems]

        self._recent = recent_paths(self.recent_model)
        self._recent_row_annotations = recent_paths(self.row_annotations_combo.model())
        self._recent_col_annotations = recent_paths(self.col_annotations_combo.model())
        self._last_path = self._current_path

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

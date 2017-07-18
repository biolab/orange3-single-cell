import os
import logging
from warnings import catch_warnings

import numpy as np
from AnyQt.QtWidgets import \
    QStyle, QComboBox, QMessageBox, QFileDialog, QGridLayout, QLabel, \
    QLineEdit, QSizePolicy as Policy
from AnyQt.QtCore import Qt, QTimer, QSize

from Orange.canvas.gui.utils import OSX_NSURL_toLocalFile
from Orange.data import StringVariable
from Orange.data.table import Table, get_sample_datasets_dir
from Orange.data.io import FileFormat, UrlReader
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting, ContextSetting, \
    PerfectDomainContextHandler, SettingProvider
from Orange.widgets.utils.domaineditor import DomainEditor
from Orange.widgets.utils.itemmodels import PyListModel
from Orange.widgets.utils.filedialogs import RecentPathsWComboMixin
from Orange.widgets.widget import Output

from Orange.widgets.utils.filedialogs import RecentPath

from orangecontrib.variants.reader import VariantData


log = logging.getLogger(__name__)


class OWVcfFile(widget.OWWidget, RecentPathsWComboMixin):
    name = "VCF File"
    id = "orangecontrib.variants.widgets.vcf"
    description = "Read data from a VCF file."
    icon = "icons/VCFFile.svg"
    priority = 10
    category = "Variants"
    keywords = ["data", "vcf", "file", "load", "read"]

    class Outputs:
        data = Output("Data", Table,
                      doc="Attribute-valued data set read from the input file.")

    want_main_area = False

    SEARCH_PATHS = [("sample-datasets", get_sample_datasets_dir())]
    SIZE_LIMIT = 1e7

    settingsHandler = PerfectDomainContextHandler(
        match_values=PerfectDomainContextHandler.MATCH_VALUES_ALL
    )

    # Overload RecentPathsWidgetMixin.recent_paths to set defaults
    recent_paths = Setting([
        RecentPath("", "sample-datasets", "small.vcf"),
    ])
    quality = Setting(1)
    cb_qual = Setting(True)
    frequency = Setting(1)
    cb_freq = Setting(True)

    class Warning(widget.OWWidget.Warning):
        file_too_big = widget.Msg("The file is too large to load automatically."
                                  " Press Reload to load.")

    class Error(widget.OWWidget.Error):
        file_not_found = widget.Msg("File not found.")

    def __init__(self):
        super().__init__()
        RecentPathsWComboMixin.__init__(self)
        self.domain = None
        self.variants = None
        self.table = None
        self.loaded_file = ""

        layout = QGridLayout()
        gui.widgetBox(self.controlArea, margin=0, orientation=layout)
        label = gui.widgetLabel(self, " File:  ")
        layout.addWidget(label, 0, 0, Qt.AlignVCenter)

        box = gui.hBox(None, addToLayout=False, margin=0)
        box.setSizePolicy(Policy.MinimumExpanding, Policy.Fixed)
        self.file_combo.setSizePolicy(Policy.MinimumExpanding, Policy.Fixed)
        self.file_combo.activated[int].connect(self.select_file)
        box.layout().addWidget(self.file_combo)
        layout.addWidget(box, 0, 1)

        file_button = gui.button(
            None, self, '...', callback=self.browse_file, autoDefault=False)
        file_button.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        file_button.setSizePolicy(Policy.Maximum, Policy.Fixed)
        layout.addWidget(file_button, 0, 2)

        reload_button = gui.button(
            None, self, "Reload", callback=self.load_data, autoDefault=False)
        reload_button.setIcon(self.style().standardIcon(
            QStyle.SP_BrowserReload))
        reload_button.setSizePolicy(Policy.Fixed, Policy.Fixed)
        layout.addWidget(reload_button, 0, 3)

        box = gui.vBox(self.controlArea, "Info")
        self.info = gui.widgetLabel(box, 'No data loaded.')
        self.warnings = gui.widgetLabel(box, '')

        def enable_apply():
            self.apply_button.setEnabled(True)

        box = gui.vBox(self.controlArea, "Filtering")
        _, qspin = gui.spin(
            box, self, 'quality', 0, 999, step=1,
            label='Quality threshold (QT)', callback=enable_apply,
            checked='cb_qual', checkCallback=enable_apply)
        qspin.setToolTip("Minimum quality to use reads.")
        _, fspin = gui.spin(
            box, self, 'frequency', 0, 999, step=1,
            label='Frequency threshold (FT)', callback=enable_apply,
            checked='cb_freq', checkCallback=enable_apply)
        fspin.setToolTip("Keep only variants with at least this many "
                         "occurrences of alternative alleles.")

        gui.rubber(self.controlArea)

        box = gui.hBox(self.controlArea)
        box.layout().addWidget(self.report_button)
        self.report_button.setFixedWidth(170)
        gui.rubber(box)

        self.apply_button = gui.button(
            box, self, "Apply", callback=self.apply)
        self.apply_button.setEnabled(False)
        self.apply_button.setFixedWidth(170)

        self.set_file_list()
        # Must not call open_file from within __init__. open_file
        # explicitly re-enters the event loop (by a progress bar)

        self.setAcceptDrops(True)

        last_path = self.last_path()
        if last_path and os.path.exists(last_path) and \
                os.path.getsize(last_path) > self.SIZE_LIMIT:
            self.Warning.file_too_big()
            return

        QTimer.singleShot(0, self.load_data)

    def sizeHint(self):
        return QSize(500, 200)

    def select_file(self, n):
        assert n < len(self.recent_paths)
        super().select_file(n)
        if self.recent_paths:
            self.load_data()
            self.set_file_list()

    def browse_file(self):
        start_file = self.last_path() or os.path.expanduser("~/")
        dialog_formats = "VCF files (*.vcf);;All files (*)"

        filename, _ = QFileDialog.getOpenFileName(
            self, 'Open Orange Data File', start_file, dialog_formats)
        if not filename:
            return
        self.add_path(filename)
        self.load_data()

    # Open a file, create data from it and send it over the data channel
    def load_data(self):
        # We need to catch any exception type since anything can happen in
        # file readers
        # pylint: disable=broad-except
        self.apply_button.setEnabled(False)
        self.clear_messages()
        self.set_file_list()
        if not self.last_path() or not os.path.exists(self.last_path()):
            if self.last_path():
                self.Error.file_not_found()
            self.Outputs.data.send(None)
            self.info.setText("No data.")
            return

        error = None

        if not error:
            with catch_warnings(record=True) as warnings:
                try:
                    variants = VariantData(self.last_path())
                except Exception as ex:
                    log.exception(ex)
                    error = ex
                self.warning(warnings[-1].message.args[0] if warnings else '')

        if error:
            self.variants = self.table = None
            self.Outputs.data.send(None)
            self.info.setText("An error occurred:\n{}".format(error))
            return

        self.loaded_file = self.last_path()
        self.variants = variants
        self.apply()  # sends data

    def update_info(self):
        pl = lambda x: '' if x == 1 else 's'
        text = ""
        if self.variants is not None:
            nsamples, nvariants = self.variants.gt.T.shape
            text += ("<p>Before filtering:<br/>" +
                     "&nbsp; {} sample{}, {} variant{}</p>").\
                format(nsamples, pl(nsamples), nvariants, pl(nvariants), )
        if self.table is not None:
            nsamples, nvariants = self.table.X.shape
            below = np.isnan(self.table.X).sum() / self.table.X.size * 100
            text += ("<p>After filtering:<br/>" +
                     "&nbsp; {} sample{}, {} variant{}<br/>" +
                     "&nbsp; {:.2f}% reads below QT</p>").\
                format(nsamples, pl(nsamples), nvariants, pl(nvariants), below)
        self.info.setText(text)

    def apply(self):
        if self.variants is None:
            self.table = None
        else:
            q = self.quality if self.cb_qual else None
            f = self.frequency if self.cb_freq else None
            self.table = self.variants.get_data(q, f)

        self.update_info()
        self.Outputs.data.send(self.table)
        self.apply_button.setEnabled(False)

    def get_widget_name_extension(self):
        _, name = os.path.split(self.loaded_file)
        return os.path.splitext(name)[0]

    def send_report(self):
        if self.table is None:
            self.report_paragraph("VCF File", "No file.")
            return
        home = os.path.expanduser("~")
        if self.loaded_file.startswith(home):
            # os.path.join does not like ~
            name = "~" + os.path.sep + \
                   self.loaded_file[len(home):].lstrip("/").lstrip("\\")
        else:
            name = self.loaded_file
        self.report_items("VCF File", [("File name", name),])
        parameters = [("Quality", self.quality, self.cb_qual),
                      ("Frequency", self.frequency, self.cb_freq)]
        self.report_items(
            "Filtering parameters",
            [(name, value) for name, value, enabled in parameters if enabled])
        self.report_data("Data", self.table)

    def dragEnterEvent(self, event):
        """Accept drops of valid file urls"""
        urls = event.mimeData().urls()
        if urls:
            try:
                FileFormat.get_reader(OSX_NSURL_toLocalFile(urls[0]) or
                                      urls[0].toLocalFile())
                event.acceptProposedAction()
            except IOError:
                pass

    def dropEvent(self, event):
        """Handle file drops"""
        urls = event.mimeData().urls()
        if urls:
            self.add_path(OSX_NSURL_toLocalFile(urls[0]) or
                          urls[0].toLocalFile())  # add first file
            self.load_data()

if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication
    a = QApplication(sys.argv)
    ow = OWVcfFile()
    ow.show()
    a.exec_()
    ow.saveSettings()

import numbers
import sys
from collections import namedtuple
from enum import IntEnum
import numpy as np

from AnyQt.QtCore import Qt
from AnyQt.QtGui import QStandardItemModel, QStandardItem
from AnyQt.QtWidgets import (
    QStyledItemDelegate, QTreeView, QHeaderView, QApplication
)

from Orange.data import Table
from Orange.widgets import gui
from Orange.widgets.settings import (
    ContextSetting, Setting, PerfectDomainContextHandler
)
from Orange.widgets.widget import Input, Output, Msg, OWWidget

from orangecontrib.single_cell.preprocess.scbnorm import (
    ScBatchScorer, SCBatchNormalizer, LINKS
)

VariableRole = next(gui.OrangeUserRole)


class IntegralDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, numbers.Number):
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter


class RealDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        data = index.data(Qt.DisplayRole)
        if isinstance(data, numbers.Number):
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter
            option.text = "?" if np.isnan(data) else "{:.3f}".format(data)


class LinkMethod(IntEnum):
    IDENTITY_LINK, LOG_LINK = range(2)

    @staticmethod
    def items():
        return sorted(list(LINKS.keys()))


class OWBatchNorm(OWWidget):
    name = "Batch Effect Removal"
    description = "Batch effect normalization on Single Cell data set."
    icon = "icons/BatchEffectRemoval.svg"
    priority = 230

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)

    class Error(OWWidget.Error):
        general_error = Msg({})
        discrete_attributes = Msg("Data with discrete attributes "
                                  "can not be processed.")

    class Warning(OWWidget.Warning):
        missing_values = Msg("Missing values have been replaced with 0.")
        negative_values = Msg("Unable to use current settings due "
                              "to negative values in data.")

    resizing_enabled = False
    want_main_area = False

    settingsHandler = PerfectDomainContextHandler()
    batch_vars = ContextSetting([])
    link_method = Setting(LinkMethod.IDENTITY_LINK)
    skip_zeros = Setting(False)
    auto_commit = Setting(True)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.data = None

        # Info
        infobox = gui.widgetBox(self.controlArea, "Info")
        self.info_label = gui.widgetLabel(infobox, "No data on input.")

        # Link method
        method_box = gui.widgetBox(self.controlArea, "Method")
        gui.comboBox(
            method_box, self, "link_method", items=LinkMethod.items(),
            callback=self.__link_method_changed
        )
        gui.separator(method_box)
        self.skip_zeros_check = gui.checkBox(
            method_box, self, "skip_zeros", "Skip zero expressions",
            enabled=self.link_method != LinkMethod.LOG_LINK,
            callback=lambda: self.commit()
        )

        # Batch Variable Selection
        header_shema = (
            ("selected", ""),
            ("variable", "Variable"),
            ("count", "#"),
            ("score", "Score")
        )
        header_labels = labels = [label for _, label in header_shema]
        header = namedtuple("header", [tag for tag, _ in header_shema])
        self.Header = header(*[index for index, _ in enumerate(labels)])

        batch_box = gui.widgetBox(self.controlArea, "Batch Variable Selection")
        self.view = QTreeView()
        self.model = QStandardItemModel()
        self.model.itemChanged.connect(self.__selected_batch_vars_changed)
        self.model.setHorizontalHeaderLabels(header_labels)
        batch_box.layout().addWidget(self.view)
        self._setup_view()

        gui.auto_commit(self.controlArea, self, "auto_commit",
                        "Apply", "Apply Automatically")

    def __link_method_changed(self):
        enable = self.link_method != LinkMethod.LOG_LINK
        self.skip_zeros_check.setEnabled(enable)
        if not enable:
            self.skip_zeros_check.setChecked(True)
        self.commit()

    def __selected_batch_vars_changed(self, item):
        if item.checkState():
            self.batch_vars.append(item.data(VariableRole))
        else:
            self.batch_vars.remove(item.data(VariableRole))
        self.commit()

    def _setup_view(self):
        self.view.setModel(self.model)
        self.view.setSelectionMode(QTreeView.NoSelection)
        self.view.setSortingEnabled(True)
        self.view.setRootIsDecorated(False)
        self.view.setItemDelegateForColumn(self.Header.count,
                                           IntegralDelegate(self))
        self.view.setItemDelegateForColumn(self.Header.score,
                                           RealDelegate(self))
        self.view.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.view.header().setStretchLastSection(False)
        self.view.header().setSectionResizeMode(self.Header.variable,
                                                QHeaderView.Stretch)
        self.view.setFocus()

    @Inputs.data
    def set_data(self, data):
        self.closeContext()
        self.clear()
        self.data = data
        self._setup_info_label()
        self._check_data()
        self.openContext(data)
        if self.data is not None:
            self.batch_vars = [data.domain[v.name] for v in self.batch_vars]
            self._setup_model()
        self.commit()

    def clear(self):
        self.batch_vars = []
        if self.model:
            n_rows = self.model.rowCount()
            self.model.removeRows(0, n_rows)

    def _setup_info_label(self):
        text = "No data on input."
        if self.data is not None:
            domain, attrs = self.data.domain, self.data.domain.attributes
            text = "{} cells, {} genes\n".format(len(self.data), len(attrs))
            text += "{} meta features".format(len(domain.metas)) \
                if len(domain.metas) else "(no meta features)"
        self.info_label.setText(text)

    def _check_data(self):
        self.clear_messages()
        if self.data and self.data.domain.has_discrete_attributes():
            self.data = None
            self.Error.discrete_attributes()
        if self.data and np.isnan(self.data.X).any():
            self.data.X = np.nan_to_num(self.data.X)
            self.Warning.missing_values()

    def _setup_model(self):
        estimator = ScBatchScorer()
        for var in self.data.domain.class_vars + self.data.domain.metas:
            if not var.is_primitive():
                continue
            try:
                score = float(estimator.score_data(self.data, var))
            except Exception:
                score = np.nan
            self.model.appendRow([
                self.__selected_item(var),
                self.__variable_item(var),
                self.__count_item(var),
                self.__score_item(score)
            ])

    def __selected_item(self, var):
        item = QStandardItem()
        item.setData(var, VariableRole)
        item.setCheckable(True)
        select = var in self.batch_vars
        item.setCheckState(Qt.Checked if select else Qt.Unchecked)
        item.setEditable(False)
        return item

    def __variable_item(self, var):
        item = QStandardItem()
        item.setData(var.name, Qt.DisplayRole)
        item.setData(gui.attributeIconDict[var], Qt.DecorationRole)
        item.setEditable(False)
        return item

    def __count_item(self, var):
        item = QStandardItem()
        if var.is_discrete:
            item.setData(len(var.values), Qt.DisplayRole)
        item.setEditable(False)
        return item

    def __score_item(self, score):
        item = QStandardItem()
        item.setData(score, Qt.DisplayRole)
        item.setEditable(False)
        return item

    def commit(self):
        data = None
        self.Error.general_error.clear()
        self.Warning.negative_values.clear()
        if self.data is not None:
            if (self.data.X < 0).any() and self.skip_zeros:
                self.Warning.negative_values()
                data = self.data
            else:
                try:
                    data = SCBatchNormalizer(
                        LinkMethod.items()[self.link_method],
                        self.skip_zeros, self.batch_vars
                    )(self.data)
                except Exception as e:
                    self.Error.general_error(str(e))
                    data = None
        self.Outputs.data.send(data)

    def send_report(self):
        method = LinkMethod.items()[self.link_method]
        if self.skip_zeros:
            method += " (Skip zero expressions)"
        variables = ", ".join([v.name for v in self.batch_vars]) \
            if self.batch_vars else "None"
        self.report_items("", [("Method", method),
                               ("Batch variable selection", variables)])


def main(args=None):
    app = QApplication(args or [])
    w = OWBatchNorm()
    w.set_data(Table("iris"))
    w.show()
    w.raise_()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()


if __name__ == "__main__":
    sys.exit(main())

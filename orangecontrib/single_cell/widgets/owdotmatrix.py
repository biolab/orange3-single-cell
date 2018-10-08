import numpy as np
from AnyQt.QtCore import Qt
from Orange.data import (DiscreteVariable, Table, Domain)
from Orange.data.filter import Values, FilterDiscrete
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.annotated_data import ANNOTATED_DATA_SIGNAL_NAME, create_annotated_table
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.sql import check_sql_input

from orangecontrib.single_cell.preprocess.clusteranalysis import ClusterAnalysis
from orangecontrib.single_cell.widgets.contingency_table import ContingencyTable


class OWDotMatrix(widget.OWWidget):
    name = "Dot Matrix"
    description = "Perform cluster analysis."
    icon = "icons/DotMatrix.svg"
    priority = 410

    class Inputs:
        data = Input("Data", Table, default=True)

    class Outputs:
        selected_data = Output("Selected Data", Table, default=True)
        annotated_data = Output(ANNOTATED_DATA_SIGNAL_NAME, Table)
        contingency = Output("Contingency Table", Table)

    GENE_MAXIMUM = 100
    CELL_SIZES = (14, 22, 30)
    AGGREGATE_F = [
        lambda x: np.mean(x, axis=0),
        lambda x: np.median(x, axis=0),
        lambda x: np.min(x, axis=0),
        lambda x: np.max(x, axis=0),
        lambda x: np.mean(x > 0, axis=0),
    ]
    AGGREGATE_NAME = [
        "Mean expression",
        "Median expression",
        "Min expression",
        "Max expression",
        "Fraction expressing"
    ]

    settingsHandler = DomainContextHandler(metas_in_res=True)
    cluster_var = ContextSetting(None)
    aggregate_ix = ContextSetting(0)  # type: int
    biclustering = ContextSetting(True)
    transpose = ContextSetting(False)
    log_scale = ContextSetting(False)
    cell_size_ix = ContextSetting(2)  # type: int
    selection = ContextSetting(set())
    auto_apply = Setting(True)

    want_main_area = True

    def __init__(self):
        super().__init__()

        self.data = None  # type: Table
        self.matrix = None
        self.clusters = None
        self.cluster_order = None
        self.genes = None
        self.gene_order = None
        self.rows = None
        self.columns = None
        self.feature_model = DomainModel(valid_types=DiscreteVariable)

        box = gui.vBox(self.controlArea, "Info")
        self.infobox = gui.widgetLabel(box, self._get_info_string())

        box = gui.vBox(self.controlArea, "Cluster Variable")
        gui.comboBox(box, self, "cluster_var", sendSelectedValue=True,
                     model=self.feature_model, callback=self._calculate_table_values)

        box = gui.vBox(self.controlArea, "Aggregation")
        gui.comboBox(box, self, "aggregate_ix", sendSelectedValue=False,
                     items=self.AGGREGATE_NAME, callback=self._calculate_table_values)

        box = gui.vBox(self.controlArea, "Options")
        gui.checkBox(box, self, "biclustering", "Biclustering of cells and genes",
                     callback=self._calculate_table_values)
        gui.checkBox(box, self, "transpose", "Transpose",
                     callback=self._refresh_table)
        gui.checkBox(box, self, "log_scale", "Log scale",
                     callback=self._refresh_table)

        box = gui.vBox(self.controlArea, "Plot Size")
        gui.radioButtons(box, self, "cell_size_ix", btnLabels=("S", "M", "L"),
                         callback=lambda: self.tableview.set_cell_size(self.CELL_SIZES[self.cell_size_ix]),
                         orientation=Qt.Horizontal)

        gui.rubber(self.controlArea)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.tableview = ContingencyTable(self)
        self.mainArea.layout().addWidget(self.tableview)

    def _get_info_string(self):
        formatstr = "{} genes, {} cells\n{} clusters"
        if self.data:
            return formatstr.format(len(self.data.domain.attributes),
                                    len(self.data),
                                    len(self.clusters))
        else:
            return formatstr.format(*([0] * 3))

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        if self.feature_model:
            self.closeContext()

        self.data = data
        self.matrix = None
        self.feature_model.set_domain(None)
        self.cluster_var = None
        self.clusters = None
        self.cluster_order = None
        self.genes = None
        self.gene_order = None
        self.rows = None
        self.columns = None

        if self.data:
            self.feature_model.set_domain(self.data.domain)
            if self.feature_model:
                self.openContext(self.data)
                if self.cluster_var is None:
                    self.cluster_var = self.feature_model[0]
                self._calculate_table_values()
            else:
                self.tableview.clear()
        else:
            self.tableview.clear()

    @staticmethod
    def _group_by(table: Table, var: DiscreteVariable):
        column = table.get_column_view(var)[0]
        return (table[column == value] for value in np.unique(column))

    def _calculate_table_values(self):
        if self.data is None:
            self.Warning.clear()
        else:
            self.clusters = [self.cluster_var.values[int(ix)]
                             for ix in np.unique(self.data.get_column_view(self.cluster_var)[0])]
            self.genes = [var.name for var in self.data.domain.attributes]
            self.infobox.setText(self._get_info_string())

            if len(self.genes) > 100:
                self.warning("Too many genes on input, first {} genes displayed.".format(self.GENE_MAXIMUM))
            else:
                self.Warning.clear()

            self.matrix = np.stack((self.AGGREGATE_F[self.aggregate_ix](cluster.X[:self.GENE_MAXIMUM])
                                    for cluster in self._group_by(self.data, self.cluster_var)),
                                   axis=0)

            if self.biclustering:
                self.cluster_order, self.gene_order = ClusterAnalysis.biclustering(self.matrix,
                                                                                   ClusterAnalysis.neighbor_distance)
            else:
                self.cluster_order, self.gene_order = np.arange(len(self.clusters)), np.arange(len(self.genes))
            self.matrix = self.matrix[self.cluster_order][:,self.gene_order]

            self._refresh_table()
            self._invalidate()

    def _refresh_table(self):
        if self.matrix is not None:
            if not self.transpose:
                self.rows, self.columns = self.clusters, self.genes
                row_order, column_order = self.cluster_order, self.gene_order
            else:
                self.rows, self.columns = self.genes, self.clusters
                row_order, column_order = self.gene_order, self.cluster_order
            self.tableview.set_headers(np.array(self.rows)[row_order], np.array(self.columns)[column_order], circles=True,
                                       cell_size=self.CELL_SIZES[self.cell_size_ix], bold_headers=False)
            if self.matrix.size > 0:
                matrix = self.matrix
                if self.log_scale:
                    matrix = np.log(matrix + 1 - matrix.min())
                matrix = matrix / matrix.max()
                if self.transpose:
                    matrix = matrix.T

                def tooltip(i,j):
                    if not self.transpose:
                        cluster, gene, value = self.clusters[i], self.genes[j], self.matrix[i,j]
                    else:
                        cluster, gene, value = self.clusters[j], self.genes[i], self.matrix[j,i]
                    return "Cluster: {}\nGene: {}\n{}: {:.1f}".format(
                        cluster, gene, self.AGGREGATE_NAME[self.aggregate_ix], value)

                self.tableview.update_table(matrix, tooltip=tooltip)

    def commit(self):
        if len(self.selection):
            cluster_ids = set()
            gene_ids = set()
            for (ir, ic) in self.selection:
                if not self.transpose:
                    cluster_ids.add(ir)
                    gene_ids.add(ic)
                else:
                    cluster_ids.add(ic)
                    gene_ids.add(ir)

            new_domain = Domain([self.data.domain[self.genes[i]] for i in gene_ids],
                                self.data.domain.class_vars,
                                self.data.domain.metas)
            selected_data = Values([FilterDiscrete(self.cluster_var, [self.clusters[i]])
                                    for i in cluster_ids],
                                   conjunction=False)(self.data)
            selected_data = selected_data.transform(new_domain)
            annotated_data = create_annotated_table(self.data.transform(new_domain),
                                                    np.where(np.in1d(self.data.ids, selected_data.ids, True)))
        else:
            selected_data = None
            annotated_data = create_annotated_table(self.data, [])
        if self.matrix is not None:
            table = ClusterAnalysis.contingency_table(self.matrix,
                                                      DiscreteVariable(self.cluster_var.name,
                                                                       np.array(self.clusters)),
                                                      np.array(self.genes)[self.gene_order],
                                                      self.cluster_order[...,np.newaxis])
        else:
            table = None
        self.Outputs.selected_data.send(selected_data)
        self.Outputs.annotated_data.send(annotated_data)
        self.Outputs.contingency.send(table)

    def _invalidate(self):
        self.selection = self.tableview.get_selection()
        self.commit()

    def send_report(self):
        rows = None
        columns = None
        if self.data is not None:
            rows = self.cluster_var
            if rows in self.data.domain:
                rows = self.data.domain[rows]
            columns = self.columns
            if columns in self.data.domain:
                columns = self.data.domain[columns]
        self.report_items((
            ("Rows", rows),
            ("Columns", columns),
        ))


def test():
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])

    w = OWDotMatrix()
    data = Table("iris")
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    app.exec_()


if __name__ == "__main__":
    test()

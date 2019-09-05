import numpy as np
from AnyQt.QtCore import Qt
from Orange.clustering import hierarchical
from Orange.clustering.hierarchical import leaves
from Orange.data import (DiscreteVariable, Table, Domain, ContinuousVariable)
from Orange.data.filter import Values, FilterDiscrete
from Orange.distance import Euclidean
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.annotated_data import ANNOTATED_DATA_SIGNAL_NAME, create_annotated_table
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.sql import check_sql_input
from Orange.preprocess import SklImpute
from Orange.widgets.widget import OWWidget, Msg

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

    class Error(OWWidget.Error):
        no_discrete_variable = Msg("No discrete variables in data.")

    class Warning(OWWidget.Warning):
        to_many_attributes = Msg("Too many genes on input, first {} genes displayed.")

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

    settingsHandler = DomainContextHandler()
    cluster_var = ContextSetting(None)
    aggregate_ix = ContextSetting(0)  # type: int
    biclustering = ContextSetting(True)
    transpose = ContextSetting(False)
    log_scale = ContextSetting(False)
    normalize = ContextSetting(True)
    cell_size_ix = ContextSetting(2)  # type: int
    selection_indices = ContextSetting(set())
    auto_apply = Setting(True)

    want_main_area = True

    def __init__(self):
        super().__init__()

        self.feature_model = DomainModel(valid_types=DiscreteVariable)
        self._init_vars()
        self._set_info_string()

        box = gui.vBox(self.controlArea, "Cluster Variable")
        gui.comboBox(box, self, "cluster_var", sendSelectedValue=True,
                     model=self.feature_model, callback=self._aggregate_data)

        box = gui.vBox(self.controlArea, "Aggregation")
        gui.comboBox(box, self, "aggregate_ix", sendSelectedValue=False,
                     items=self.AGGREGATE_NAME, callback=self._aggregate_data)

        box = gui.vBox(self.controlArea, "Options")
        gui.checkBox(box, self, "biclustering", "Order cells and genes",
                     callback=self._calculate_table_values)
        gui.checkBox(box, self, "transpose", "Transpose",
                     callback=self._calculate_table_values)
        gui.checkBox(box, self, "log_scale", "Log scale",
                     callback=self._calculate_table_values)
        gui.checkBox(box, self, "normalize", "Normalize data",
                     callback=self._calculate_table_values)

        box = gui.vBox(self.controlArea, "Plot Size")
        gui.radioButtons(box, self, "cell_size_ix", btnLabels=("S", "M", "L"),
                         callback=lambda: self.tableview.set_cell_size(self.CELL_SIZES[self.cell_size_ix]),
                         orientation=Qt.Horizontal)

        gui.rubber(self.controlArea)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.tableview = ContingencyTable(self)
        self.mainArea.layout().addWidget(self.tableview)

    def _init_vars(self):
        self.data = None  # type: Table
        self.matrix = None
        self.clusters = None
        self.clusters_unordered = None
        self.aggregated_data = None
        self.cluster_var = None
        self.selected_names = {}

    def _set_info_string(self):
        formatstr = "{} genes\n{} cells\n{} clusters"
        if self.data:
            self.info.set_input_summary(
                str(len(self.data)),
                formatstr.format(len(self.data.domain.attributes), len(self.data), len(self.clusters_unordered)))
        else:
            self.info.set_input_summary(self.info.NoInput)

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        if self.feature_model:
            self.closeContext()

        self._init_vars()
        self.Error.no_discrete_variable.clear()
        self.Warning.clear()
        self.data = data
        if self.data:
            self.feature_model.set_domain(self.data.domain)
            if self.feature_model:
                self.openContext(self.data)
                if self.cluster_var is None:
                    self.cluster_var = self.feature_model[0]
                self._aggregate_data()
            else:
                self.tableview.clear()
                self.Error.no_discrete_variable()
                self.data = None
                self._set_info_string()
                self.commit()
        else:
            self.tableview.clear()
            self._set_info_string()
            self.commit()

    @staticmethod
    def _group_by(table: Table, var: DiscreteVariable):
        column = table.get_column_view(var)[0]
        return (table[column == value] for value in np.unique(column))

    @staticmethod
    def _transpose(matrix: np.ndarray, clusters, genes):
        return matrix.T, np.array([str(g) for g in genes]), [ContinuousVariable(c) for c in clusters]

    @staticmethod
    def _normalize(matrix: np.ndarray):
        matrix = (matrix - np.mean(matrix, axis=0, keepdims=True)) / (
                np.std(matrix, axis=0, keepdims=True) + 1e-10)
        matrix[matrix < -3] = -3
        matrix[matrix > 3] = 3
        return matrix

    @staticmethod
    def _norm_min_max(matrix: np.ndarray):
        matrix = matrix - matrix.min()
        return matrix / (matrix.max() + 1e-12)

    def _aggregate_data(self):
        self.Warning.clear()
        if self.data is None:
            return

        self.clusters_unordered = np.array(
            [self.cluster_var.values[int(ix)]
             for ix in np.unique(self.data.get_column_view(self.cluster_var)[0])])
        self._set_info_string()

        if len(self.data.domain.attributes) > self.GENE_MAXIMUM:
            self.Warning.to_many_attributes(self.GENE_MAXIMUM)

        self.aggregated_data = np.stack([self.AGGREGATE_F[self.aggregate_ix](cluster.X[:, :self.GENE_MAXIMUM])
                                         for cluster in self._group_by(self.data, self.cluster_var)],
                                        axis=0)
        self._calculate_table_values()

    def _calculate_table_values(self):
        genes = self.data.domain.attributes[:self.GENE_MAXIMUM]
        matrix = self.aggregated_data
        clusters = self.clusters_unordered
        if self.transpose:
            matrix, clusters, genes = self._transpose(matrix, clusters, genes)

        # create data table since imputation of nan values is required
        matrix = Table(Domain(genes), matrix)
        matrix_before_norm = matrix.copy()  # for tooltip
        matrix = SklImpute()(matrix)

        if self.log_scale:
            matrix.X = np.log(matrix.X + 1)
        if self.normalize:
            matrix.X = self._normalize(matrix.X)

        # values must be in range [0, 1] for visualisation
        matrix.X = self._norm_min_max(matrix.X)

        if self.biclustering:
            cluster_order, gene_order = self.cluster_data(matrix)
        else:
            cluster_order, gene_order = np.arange(matrix.X.shape[0]), np.arange(matrix.X.shape[1])

        # reorder
        self.matrix = matrix[cluster_order][:, gene_order]
        self.matrix_before_norm = matrix_before_norm[cluster_order][:, gene_order]
        self.clusters = clusters[cluster_order]

        self._refresh_table()
        self._update_selection()
        self._invalidate()

    def cluster_data(self, matrix):
        with self.progressBar():
            # cluster rows
            if len(matrix) > 1:
                rows_distances = Euclidean(matrix)
                cluster = hierarchical.dist_matrix_clustering(rows_distances)
                row_order = hierarchical.optimal_leaf_ordering(
                    cluster, rows_distances, progress_callback=self.progressBarSet)
                row_order = np.array([x.value.index for x in leaves(row_order)])
            else:
                row_order = np.array([0])

            # cluster columns
            if matrix.X.shape[1] > 1:
                columns_distances = Euclidean(matrix, axis=0)
                cluster = hierarchical.dist_matrix_clustering(columns_distances)
                columns_order = hierarchical.optimal_leaf_ordering(
                    cluster, columns_distances,
                    progress_callback=self.progressBarSet)
                columns_order = np.array([x.value.index for x in leaves(columns_order)])
            else:
                columns_order = np.array([0])
        return row_order, columns_order

    def _refresh_table(self):
        if self.matrix is None:
            return

        columns = np.array([str(x) for x in self.matrix.domain.attributes])
        rows = self.clusters
        # row_order, column_order = self.gene_order, self.cluster_order
        self.tableview.set_headers(rows, columns, circles=True,
                                   cell_size=self.CELL_SIZES[self.cell_size_ix], bold_headers=False)
        if self.matrix.X.size > 0:
            matrix = self.matrix.X

            def tooltip(i, j):
                cluster, gene, value = rows[i], columns[j], self.matrix_before_norm[i, j]
                return "Cluster: {}\nGene: {}\n{}: {:.1f}".format(
                    cluster, gene, self.AGGREGATE_NAME[self.aggregate_ix], value)

            self.tableview.update_table(matrix, tooltip=tooltip)

    def commit(self):
        if self.data is None:
            self.Outputs.selected_data.send(None)
            self.Outputs.annotated_data.send(None)
            self.Outputs.contingency.send(None)
            return

        if len(self.selection_indices):
            cluster_ids = set()
            gene_ids = set()
            for (ir, ic) in self.selection_indices:
                if not self.transpose:
                    cluster_ids.add(ir)
                    gene_ids.add(ic)
                else:
                    cluster_ids.add(ic)
                    gene_ids.add(ir)

            columns = self.clusters if self.transpose else [str(x) for x in self.matrix.domain.attributes]
            rows = self.clusters if not self.transpose else [str(x) for x in self.matrix.domain.attributes]
            new_domain = Domain([self.data.domain[columns[i]] for i in gene_ids],
                                self.data.domain.class_vars,
                                self.data.domain.metas)
            selected_data = Values([FilterDiscrete(self.cluster_var, [rows[i]])
                                    for i in cluster_ids],
                                   conjunction=False)(self.data)
            selected_data = selected_data.transform(new_domain)
            annotated_data = create_annotated_table(self.data,
                                                    np.where(np.in1d(self.data.ids, selected_data.ids, True)))
        else:
            selected_data = None
            annotated_data = create_annotated_table(self.data, [])

        clusters_values = list(set(self.clusters))
        table = ClusterAnalysis.contingency_table(
            self.matrix,
            DiscreteVariable("Gene" if self.transpose else self.cluster_var.name, clusters_values),
            [str(x) for x in self.matrix.domain.attributes],
            [[clusters_values.index(c)] for c in self.clusters]
        )

        self.Outputs.selected_data.send(selected_data)
        self.Outputs.annotated_data.send(annotated_data)
        self.Outputs.contingency.send(table)

    def _update_selection(self):
        """
        This function updates widget selection in case when any item is selected.
        It updates selection when order has changed
        """
        rows = self.clusters.tolist()
        columns = [str(x) for x in self.matrix.domain.attributes]
        if self.transpose:
            new_selection = {(rows.index(g), columns.index(c)) for g, c in self.selected_names}
        else:
            new_selection = {(rows.index(c), columns.index(g)) for g, c in self.selected_names}
        self.tableview.set_selection(new_selection)

    def _invalidate(self):
        self.save_selection_names()
        self.commit()

    def save_selection_names(self):
        """
        With this method we save the names of selected genes-clusters pairs, since options changes
        the columns, rows orders and we want to keep the selection.
        """
        self.selection_indices = self.tableview.get_selection()
        genes = self.clusters if self.transpose else [str(x) for x in self.matrix.domain.attributes]
        clusters = self.clusters if not self.transpose else [str(x) for x in self.matrix.domain.attributes]

        if self.transpose:
            self.selected_names = {(genes[g], clusters[c]) for g, c in self.selection_indices}
        else:
            self.selected_names = {(genes[g], clusters[c]) for c, g in self.selection_indices}


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

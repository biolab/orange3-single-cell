from AnyQt.QtGui import QStandardItemModel
from Orange.data import (DiscreteVariable, Table)
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.annotated_data import ANNOTATED_DATA_SIGNAL_NAME
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.sql import check_sql_input

from orangecontrib.single_cell.widgets.contingency_table import ContingencyTable
from orangecontrib.single_cell.widgets.clusteranalysis import ClusterAnalysis


class OWClusterAnalysis(widget.OWWidget):
    name = "Cluster Analysis"
    description = "Perform cluster analysis."
    priority = 2011

    inputs = [("Data", Table, "set_data", widget.Default)]
    outputs = [("Selected Data", Table),
               (ANNOTATED_DATA_SIGNAL_NAME, Table),
               ("Contingency Table", Table)]

    settingsHandler = DomainContextHandler(metas_in_res=True)
    rows = ContextSetting(None)
    columns = ContextSetting(None)
    clustering_var = ContextSetting(None)
    selection = ContextSetting(set())
    auto_apply = Setting(True)

    want_main_area = True

    def __init__(self):
        super().__init__()

        self.data = None
        self.feature_model = DomainModel(valid_types=DiscreteVariable)
        self.table = None

        box = gui.vBox(self.controlArea, "Info")
        self.infobox = gui.widgetLabel(box, self._get_info_string(None))

        box = gui.vBox(self.controlArea, "Rows")
        gui.comboBox(box, self, "clustering_var", sendSelectedValue=True,
                     model=self.feature_model, callback=self._run_cluster_analysis)

        gui.rubber(self.controlArea)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.tablemodel = QStandardItemModel(self)
        self.tableview = ContingencyTable(self, self.tablemodel)
        self.mainArea.layout().addWidget(self.tableview)

    def _get_info_string(self, cluster_variable):
        formatstr = "Cells: {0}\nGenes: {1}\nClusters: {2}"
        if self.data:
            return formatstr.format(len(self.data),
                                    len(self.data.domain.attributes),
                                    len(self.data.domain[cluster_variable].values))
        else:
            return formatstr.format(*["No input data"]*3)

    @check_sql_input
    def set_data(self, data):
        if self.feature_model:
            self.closeContext()
        self.data = data
        self.feature_model.set_domain(None)
        self.rows = None
        self.columns = None
        if self.data:
            self.feature_model.set_domain(self.data.domain)
            if self.feature_model:
                self.clustering_var = self.feature_model[0]
                self._run_cluster_analysis()
            else:
                self.tablemodel.clear()
        else:
            self.tablemodel.clear()

    def _run_cluster_analysis(self):
        self.infobox.setText(self._get_info_string(self.clustering_var.name))
        CA = ClusterAnalysis(self.data, n_enriched=2, genes=[1, 2, 3, 4], clustering_var=self.clustering_var.name)
        CA.percentage_expressing()
        self.table = CA.sort_percentage_expressing()
        # Referencing the Cluster variable directly doesn't preserve the order of clusters.
        clusters = [self.clustering_var.values[ix] for ix in self.table.get_column_view("Cluster")[0]]
        self.rows = DiscreteVariable("Cluster", clusters, ordered=True)
        self.columns = DiscreteVariable("Gene", [var.name for var in self.table.domain.variables], ordered=True)
        self.tableview.set_variables(self.rows, self.columns)
        self.tableview.update_table(self.table.X, formatstr="{:.2f}")
        self._invalidate()

    def handleNewSignals(self):
        self._invalidate()

    def commit(self):
        # TODO: Implement outputs for Selected Data and Data here.
        self.send("Contingency Table", self.table)
        pass

    def _invalidate(self):
        self.selection = self.tableview.get_selection()
        self.commit()

    def send_report(self):
        rows = None
        columns = None
        if self.data is not None:
            rows = self.rows
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

    w = OWClusterAnalysis()
    data = Table("testdata.tab")
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    app.exec_()


if __name__ == "__main__":
    test()

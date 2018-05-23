import concurrent.futures
from functools import partial

import numpy as np
from AnyQt.QtCore import Qt, Slot, QThread
from AnyQt.QtWidgets import QGridLayout
from Orange.data import (DiscreteVariable, Table, Domain)
from Orange.data.filter import Values, FilterDiscrete
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.annotated_data import ANNOTATED_DATA_SIGNAL_NAME, create_annotated_table
from Orange.widgets.utils.concurrent import ThreadExecutor, FutureWatcher
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.sql import check_sql_input
from orangecontrib.bioinformatics.widgets.utils.data import GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE

from orangecontrib.single_cell.preprocess.clusteranalysis import ClusterAnalysis
from orangecontrib.single_cell.widgets.contingency_table import ContingencyTable


class Task:
    future = None
    watcher = None
    cancelled = False

    def __init__(self, type):
        self.type = type

    def cancel(self):
        self.cancelled = True
        self.future.cancel()
        concurrent.futures.wait([self.future])


class OWClusterAnalysis(widget.OWWidget):
    name = "Cluster Analysis"
    description = "Perform cluster analysis."
    icon = "icons/ClusterAnalysis.svg"
    priority = 2011

    inputs = [("Data", Table, "set_data", widget.Default),
              ("Genes", Table, "set_genes")]
    outputs = [("Selected Data", Table),
               (ANNOTATED_DATA_SIGNAL_NAME, Table),
               ("Contingency Table", Table)]

    N_GENES_PER_CLUSTER_MAX = 10
    N_MOST_ENRICHED_MAX = 50

    settingsHandler = DomainContextHandler(metas_in_res=True)
    cluster_var = ContextSetting(None)
    selection = ContextSetting(set())
    gene_selection = ContextSetting(0)
    differential_expression = ContextSetting(0)
    _diff_exprs = ("high", "low", "either")
    n_genes_per_cluster = ContextSetting(3)
    n_most_enriched = ContextSetting(20)
    biclustering = ContextSetting(True)
    auto_apply = Setting(True)

    want_main_area = True

    def __init__(self):
        super().__init__()

        self.ca = None
        self.clusters = None
        self.data = None
        self.feature_model = DomainModel(valid_types=DiscreteVariable)
        self.gene_list = None
        self.model = None
        self.pvalues = None

        self._executor = ThreadExecutor()
        self._gene_selection_history = (self.gene_selection, self.gene_selection)
        self._task = None

        box = gui.vBox(self.controlArea, "Info")
        self.infobox = gui.widgetLabel(box, self._get_info_string())

        box = gui.vBox(self.controlArea, "Cluster Variable")
        gui.comboBox(box, self, "cluster_var", sendSelectedValue=True,
                     model=self.feature_model, callback=self._run_cluster_analysis)

        layout = QGridLayout()
        self.gene_selection_radio_group = gui.radioButtonsInBox(
            self.controlArea, self, "gene_selection", orientation=layout,
            box="Gene Selection", callback=self._gene_selection_changed)

        def conditional_set_gene_selection(id):
            def f():
                if self.gene_selection == id:
                    return self._set_gene_selection()

            return f

        layout.addWidget(gui.appendRadioButton(self.gene_selection_radio_group, "", addToLayout=False), 1, 1)
        cb = gui.hBox(None, margin=0)
        gui.widgetLabel(cb, "Top")
        self.n_genes_per_cluster_spin = gui.spin(
            cb, self, "n_genes_per_cluster", minv=1, maxv=self.N_GENES_PER_CLUSTER_MAX,
            controlWidth=60, alignment=Qt.AlignRight, callback=conditional_set_gene_selection(0))
        gui.widgetLabel(cb, "genes per cluster")
        gui.rubber(cb)
        layout.addWidget(cb, 1, 2, Qt.AlignLeft)

        layout.addWidget(gui.appendRadioButton(self.gene_selection_radio_group, "", addToLayout=False), 2, 1)
        mb = gui.hBox(None, margin=0)
        gui.widgetLabel(mb, "Top")
        self.n_most_enriched_spin = gui.spin(
            mb, self, "n_most_enriched", minv=1, maxv=self.N_MOST_ENRICHED_MAX,
            controlWidth=60, alignment=Qt.AlignRight, callback=conditional_set_gene_selection(1))
        gui.widgetLabel(mb, "highest enrichments")
        gui.rubber(mb)
        layout.addWidget(mb, 2, 2, Qt.AlignLeft)

        layout.addWidget(gui.appendRadioButton(self.gene_selection_radio_group, "", addToLayout=False, disabled=True),
                         3, 1)
        sb = gui.hBox(None, margin=0)
        gui.widgetLabel(sb, "User-provided list of genes")
        gui.rubber(sb)
        layout.addWidget(sb, 3, 2)

        layout = QGridLayout()
        self.differential_expression_radio_group = gui.radioButtonsInBox(
            self.controlArea, self, "differential_expression", orientation=layout,
            box="Differential Expression", callback=self._set_gene_selection)

        layout.addWidget(gui.appendRadioButton(self.differential_expression_radio_group,
                                               "Overexpressed in cluster", addToLayout=False), 1, 1)
        layout.addWidget(gui.appendRadioButton(self.differential_expression_radio_group,
                                               "Underexpressed in cluster", addToLayout=False), 2, 1)
        layout.addWidget(gui.appendRadioButton(self.differential_expression_radio_group,
                                               "Either", addToLayout=False), 3, 1)

        box = gui.vBox(self.controlArea, "Sorting")
        gui.checkBox(box, self, "biclustering", "Biclustering of analysis results", callback=self._set_gene_selection)

        gui.rubber(self.controlArea)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.tableview = ContingencyTable(self)
        self.mainArea.layout().addWidget(self.tableview)

    def _get_current_gene_selection(self):
        return self._gene_selection_history[0]

    def _get_previous_gene_selection(self):
        return self._gene_selection_history[1]

    def _progress_gene_selection_history(self, new_gene_selection):
        self._gene_selection_history = (new_gene_selection, self._gene_selection_history[0])

    def _get_info_string(self):
        formatstr = "Cells: {0}\nGenes: {1}\nClusters: {2}"
        if self.data:
            return formatstr.format(len(self.data),
                                    len(self.data.domain.attributes),
                                    len(self.cluster_var.values))
        else:
            return formatstr.format(*["No input data"] * 3)

    @check_sql_input
    def set_data(self, data):
        if self.feature_model:
            self.closeContext()
        self.data = data
        self.feature_model.set_domain(None)
        self.ca = None
        self.cluster_var = None
        self.columns = None
        self.clusters = None
        self.gene_list = None
        self.model = None
        self.pvalues = None
        self.n_genes_per_cluster_spin.setMaximum(self.N_GENES_PER_CLUSTER_MAX)
        self.n_most_enriched_spin.setMaximum(self.N_MOST_ENRICHED_MAX)
        if self.data:
            self.feature_model.set_domain(self.data.domain)
            if self.feature_model:
                self.openContext(self.data)
                if self.cluster_var is None:
                    self.cluster_var = self.feature_model[0]
                self._run_cluster_analysis()
            else:
                self.tableview.clear()
        else:
            self.tableview.clear()

    def set_genes(self, data):
        self.Error.clear()
        gene_list_radio = self.gene_selection_radio_group.group.buttons()[2]

        if (data is None
                or GENE_AS_ATTRIBUTE_NAME not in data.attributes
                or not data.attributes[GENE_AS_ATTRIBUTE_NAME] and GENE_ID_COLUMN not in data.attributes
                or data.attributes[GENE_AS_ATTRIBUTE_NAME] and GENE_ID_ATTRIBUTE not in data.attributes):
            if data is not None:
                self.error("Gene annotations missing in the input data. Use Gene Name Matching widget.")
            self.gene_list = None
            gene_list_radio.setDisabled(True)
            if self.gene_selection == 2:
                self.gene_selection_radio_group.group.buttons()[self._get_previous_gene_selection()].click()
        else:
            if data.attributes[GENE_AS_ATTRIBUTE_NAME]:
                gene_id_attribute = data.attributes.get(GENE_ID_ATTRIBUTE, None)

                self.gene_list = tuple(str(var.attributes[gene_id_attribute]) for var in data.domain.attributes
                                       if gene_id_attribute in var.attributes
                                       and var.attributes[gene_id_attribute] != "?")
            else:
                gene_id_column = data.attributes.get(GENE_ID_COLUMN, None)
                self.gene_list = tuple(str(v) for v in data.get_column_view(gene_id_column)[0]
                                       if v not in ("", "?"))
            gene_list_radio.setDisabled(False)
            if self.gene_selection == 2:
                self._set_gene_selection()
            else:
                gene_list_radio.click()

    def _run_cluster_analysis(self):
        self.infobox.setText(self._get_info_string())
        gene_count = len(self.data.domain.attributes)
        cluster_count = len(self.cluster_var.values)
        self.n_genes_per_cluster_spin.setMaximum(min(self.N_GENES_PER_CLUSTER_MAX, gene_count // cluster_count))
        self.n_most_enriched_spin.setMaximum(min(self.N_MOST_ENRICHED_MAX, gene_count))
        # TODO: what happens if error occurs? If CA fails, widget should properly handle it.
        self._start_task_init(partial(ClusterAnalysis, self.data, self.cluster_var.name))

    def _start_task_init(self, f):
        if self._task is not None:
            self.cancel()
        assert self._task is None

        self._task = Task("init")

        def callback(finished):
            if self._task.cancelled:
                raise KeyboardInterrupt()
            self.progressBarSet(finished * 50)

        f = partial(f, callback=callback)

        self.progressBarInit()
        self._task.future = self._executor.submit(f)
        self._task.watcher = FutureWatcher(self._task.future)
        self._task.watcher.done.connect(self._init_task_finished)

    def _start_task_gene_selection(self, f):
        if self._task is not None:
            self.cancel()
        assert self._task is None

        self._task = Task("gene_selection")

        def callback(finished):
            if self._task.cancelled:
                raise KeyboardInterrupt()
            self.progressBarSet(50 + finished * 50)

        f = partial(f, callback=callback)

        self.progressBarInit()
        self.progressBarSet(50)
        self._task.future = self._executor.submit(f)
        self._task.watcher = FutureWatcher(self._task.future)
        self._task.watcher.done.connect(self._gene_selection_task_finished)

    @Slot(concurrent.futures.Future)
    def _init_task_finished(self, f):
        assert self.thread() is QThread.currentThread()
        assert self._task is not None
        assert self._task.future is f
        assert f.done()

        self._task = None
        self.progressBarFinished()

        self.ca = f.result()
        self._set_gene_selection()

    @Slot(concurrent.futures.Future)
    def _gene_selection_task_finished(self, f):
        assert self.thread() is QThread.currentThread()
        assert self._task is not None
        assert self._task.future is f
        assert f.done()

        self._task = None
        self.progressBarFinished()

        self.clusters, genes, self.model, self.pvalues = f.result()
        genes = [str(gene) for gene in genes]
        self.columns = DiscreteVariable("Gene", genes, ordered=True)
        self.tableview.set_headers(self.clusters, self.columns.values, circles=True, bold_headers=False)

        def tooltip(i, j):
            return ("<b>cluster</b>: {}<br /><b>gene</b>: {}<br /><b>fraction expressing</b>: {:.2f}<br />\
                                <b>p-value</b>: {:.2e}".format(
                self.clusters[i],
                self.columns.values[j],
                self.model[i, j],
                self.pvalues[i, j])
            )

        self.tableview.update_table(self.model, tooltip=tooltip)
        self._invalidate()

    def cancel(self):
        """
        Cancel the current task (if any).
        """
        if self._task is not None:
            self._task.cancel()
            assert self._task.future.done()
            # disconnect the `_task_finished` slot
            if self._task.type == "init":
                self._task.watcher.done.disconnect(self._init_task_finished)
            else:
                self._task.watcher.done.disconnect(self._gene_selection_task_finished)
            self._task = None

    def onDeleteWidget(self):
        self.cancel()
        super().onDeleteWidget()

    def _gene_selection_changed(self):
        if self.gene_selection != self._get_current_gene_selection():
            self._progress_gene_selection_history(self.gene_selection)
            self.differential_expression_radio_group.setDisabled(self.gene_selection == 2)
            self._set_gene_selection()

    def _set_gene_selection(self):
        self.Warning.clear()
        if self.ca is not None and (self._task is None or self._task.type != "init"):
            if self.gene_selection == 0:
                f = partial(self.ca.enriched_genes_per_cluster, self.n_genes_per_cluster)
            elif self.gene_selection == 1:
                f = partial(self.ca.enriched_genes_data, self.n_most_enriched)
            else:
                if self.data is not None and GENE_ID_ATTRIBUTE not in self.data.attributes:
                    self.error("Gene annotations missing in the input data. Use Gene Name Matching widget.")
                    if self.gene_selection == 2:
                        self.gene_selection_radio_group.group.buttons()[self._get_previous_gene_selection()].click()
                    return
                relevant_genes = tuple(self.ca.intersection(self.gene_list))
                if len(relevant_genes) > self.N_MOST_ENRICHED_MAX:
                    self.warning("Only first {} reference genes shown.".format(self.N_MOST_ENRICHED_MAX))
                f = partial(self.ca.enriched_genes, relevant_genes[:self.N_MOST_ENRICHED_MAX])
            f = partial(f, enrichment=self._diff_exprs[self.differential_expression], biclustering=self.biclustering)
            self._start_task_gene_selection(f)
        else:
            self._invalidate()

    def handleNewSignals(self):
        self._invalidate()

    def commit(self):
        if len(self.selection):
            cluster_ids = set()
            column_ids = set()
            for (ir, ic) in self.selection:
                cluster_ids.add(ir)
                column_ids.add(ic)
            new_domain = Domain([self.data.domain[self.columns.values[col]] for col in column_ids],
                                self.data.domain.class_vars,
                                self.data.domain.metas)
            selected_data = Values([FilterDiscrete(self.cluster_var, [self.clusters[ir]])
                                    for ir in cluster_ids],
                                   conjunction=False)(self.data)
            selected_data = selected_data.transform(new_domain)
            annotated_data = create_annotated_table(self.data.transform(new_domain),
                                                    np.where(np.in1d(self.data.ids, selected_data.ids, True)))
            table = self.ca.create_contingency_table()
        else:
            selected_data = None
            annotated_data = create_annotated_table(self.data, [])
            table = None
        self.send("Selected Data", selected_data)
        self.send(ANNOTATED_DATA_SIGNAL_NAME, annotated_data)
        self.send("Contingency Table", table)

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

    w = OWClusterAnalysis()
    data = Table("../../../../testdata.tab")
    data.X = data.X > 0
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    app.exec_()


if __name__ == "__main__":
    test()

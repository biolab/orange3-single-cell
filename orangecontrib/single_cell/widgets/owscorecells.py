import numpy as np

from Orange.data import ContinuousVariable, StringVariable, Domain, Table
from Orange.widgets import widget, gui
from Orange.widgets.settings import (Setting, ContextSetting,
                                     DomainContextHandler)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg


class OWScoreCells(widget.OWWidget):
    name = "Score Cells"
    description = "Add a cell score based on the given set of genes"
    icon = "icons/ScoreCells.svg"
    priority = 180

    settingsHandler = DomainContextHandler(metas_in_res=True)
    gene = ContextSetting(None)
    auto_apply = Setting(True)

    want_main_area = False

    class Error(OWWidget.Error):
        no_genes = Msg("No matching genes in data")

    class Warning(OWWidget.Warning):
        some_genes = Msg("{} (of {}) genes not found in data")

    class Inputs:
        data = Input("Data", Table)
        genes = Input("Genes", Table)

    class Outputs:
        data = Output("Data", Table)

    def __init__(self):
        super().__init__()

        self.data = None
        self.genes = None
        self.feature_model = DomainModel(valid_types=StringVariable)

        box = gui.vBox(self.controlArea, "Gene name")
        gui.comboBox(box, self, 'gene', sendSelectedValue=True,
                     model=self.feature_model, callback=self._invalidate)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        self.data = data

    @Inputs.genes
    @check_sql_input
    def set_genes(self, genes):
        self.closeContext()
        self.genes = genes
        self.feature_model.set_domain(None)
        self.gene = None
        if self.genes:
            self.feature_model.set_domain(self.genes.domain)
            if self.feature_model:
                self.gene = self.feature_model[0]
                self.openContext(genes)

    def handleNewSignals(self):
        self._invalidate()

    def commit(self):
        table = None
        self.Error.no_genes.clear()
        self.Warning.some_genes.clear()
        if self.data and self.genes and self.gene:
            available_genes = set(f.name for f in self.data.domain.variables)
            gene_list_all = [str(ins[self.gene]) for ins in self.genes]
            gene_list = [g for g in gene_list_all if g in available_genes]
            if not gene_list:
                self.Error.no_genes()
                score = np.zeros(len(self.data))
            else:
                if len(gene_list) < len(gene_list_all):
                    self.Warning.some_genes(len(gene_list_all) - len(gene_list),
                                            len(gene_list_all))
                values = self.data[:, gene_list].X
                score = np.nanmax(values, axis=1)
            d = self.data.domain
            score_var = ContinuousVariable('Score')
            dom = Domain(d.attributes, d.class_vars,
                         d.metas + (score_var,))
            table = self.data.transform(dom)
            col, sparse = table.get_column_view(score_var)
            col[:] = score
        self.Outputs.data.send(table)

    def _invalidate(self):
        self.commit()

    def send_report(self):
        gene = None
        if self.genes is not None:
            gene = self.gene
            if gene in self.genes.domain:
                gene = self.genes.domain[gene]
        self.report_items((
            ("Gene", gene),
        ))


def test():
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])

    w = OWScoreCells()
    data = Table("iris")
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    app.exec_()


if __name__ == "__main__":
    test()

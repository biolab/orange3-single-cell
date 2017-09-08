from Orange.data import ContinuousVariable, StringVariable, Domain, Table
from Orange.widgets import widget, gui
from Orange.widgets.settings import (Setting, ContextSetting,
                                     DomainContextHandler)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.sql import check_sql_input


class OWScoreCells(widget.OWWidget):
    name = "Score Cells"
    description = "Add a cell score based on the given set of genes"
    icon = "icons/ScoreCells.svg"
    priority = 112

    settingsHandler = DomainContextHandler(metas_in_res=True)
    gene = ContextSetting(None)
    auto_apply = Setting(True)

    want_main_area = False

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
        if self.feature_model:
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
        if self.data and self.genes and self.gene:
            gene_list = [str(ins[self.gene]) for ins in self.genes]
            values = self.data[:, gene_list].X
            score = values.max(axis=1)
            d = self.data.domain
            score_var = ContinuousVariable('Gene score')
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

import numpy as np

from AnyQt.QtCore import QSize

from Orange.data import ContinuousVariable, StringVariable, Domain, Table
from Orange.widgets import widget, gui
from Orange.widgets.settings import (Setting, ContextSetting,
                                     DomainContextHandler)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg

from orangecontrib.bioinformatics.ncbi.gene.config import MAP_GENE_ID


class OWScoreCells(widget.OWWidget):
    name = "Score Cells"
    description = "Add a cell score based on the given set of genes"
    icon = "icons/ScoreCells.svg"
    priority = 180

    settingsHandler = DomainContextHandler()
    gene = ContextSetting(None)
    auto_apply = Setting(True)

    want_main_area = False

    UserAdviceMessages = [
        widget.Message("Gene variables must have 'gene_id' attribute! "
                       "See example data in datasets widget.",
                       'gene_id')]

    class Error(OWWidget.Error):
        no_gene_ids = Msg("gene_id not found")
        column_name = Msg("Set column name")

    class Warning(OWWidget.Warning):
        no_genes = Msg("No matching genes in data")
        imperfect_match = Msg("Some input genes can't be matched")

    class Inputs:
        data = Input("Data", Table)
        genes = Input("Genes", Table)

    class Outputs:
        data = Output("Data", Table)

    def __init__(self):
        super().__init__()

        self.data = None
        self.data_genes = set()

        self.genes = None
        self.marker_genes = set()
        self.feature_model = DomainModel(valid_types=StringVariable)

        self.score_variable_name = ''

        box = gui.vBox(self.controlArea, "Gene name: ")
        gui.comboBox(box, self, 'gene', sendSelectedValue=True,
                     model=self.feature_model, callback=self._invalidate)

        box = gui.vBox(self.controlArea, "Score output: ")
        self.line_edit = gui.lineEdit(box, self, 'score_variable_name', callback=self._invalidate)
        self.line_edit.setPlaceholderText('Column name ...')

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.sizeHint()

    def sizeHint(self):
        return QSize(320, 240)

    def handleNewSignals(self):
        self._invalidate()

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        self.data = data
        if self.data:
            for variable in self.data.domain.attributes:
                self.data_genes.add(variable.attributes.get(MAP_GENE_ID, None))

        self.data_genes = set(filter(None, self.data_genes))
        if not self.data_genes:
            self.Error.no_gene_ids()

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
                self.marker_genes = set(self.genes.get_column_view(self.gene)[0])

    def __score_cells(self):
        score = np.zeros(len(self.data))
        matched = [gene_id for gene_id in self.data_genes.intersection(self.marker_genes)]

        if not matched:
            self.Warning.no_genes()
        else:
            values = self.data[:, matched].X
            score = np.nanmax(values, axis=1)

        return score

    def commit(self):
        if self.data is None:
            self.Outputs.data.send(None)
            return

        if not self.score_variable_name:
            self.Error.column_name()
            return

        score_var = ContinuousVariable(self.score_variable_name)

        score = self.__score_cells()
        domain = Domain(self.data.domain.attributes,
                        self.data.domain.class_vars,
                        self.data.domain.metas + (score_var,))

        table = self.data.transform(domain)
        col, sparse = table.get_column_view(score_var)
        col[:] = score

        self.Outputs.data.send(table)
        self.clear_messages()

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

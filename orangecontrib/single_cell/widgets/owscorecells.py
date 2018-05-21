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

from orangecontrib.bioinformatics.widgets.utils.data import (
    GENE_ID_ATTRIBUTE, GENE_ID_COLUMN, GENE_AS_ATTRIBUTE_NAME, TAX_ID
)


class OWScoreCells(widget.OWWidget):
    name = "Score Cells"
    description = "Add a cell score based on the given set of genes"
    icon = "icons/ScoreCells.svg"
    priority = 180

    settingsHandler = DomainContextHandler()
    gene = ContextSetting(None)
    auto_apply = Setting(True)

    want_main_area = False

    class Warning(OWWidget.Warning):
        no_genes = Msg("No matching genes in data")
        some_genes = Msg("{} (of {}) genes not found in data")
        no_gene_id_attribute = Msg('Unable to locate gene id attribute')
        no_gene_id_column = Msg('Unable to locate gene id column')

    class Inputs:
        data = Input("Data", Table)
        genes = Input("Genes", Table)

    class Outputs:
        data = Output("Data", Table)

    def __init__(self):
        super().__init__()

        self.data = None
        self.map_var_to_attribute = {}
        self.data_genes = set()

        self.genes = None
        self.marker_genes = set()
        self.feature_model = DomainModel(valid_types=StringVariable)

        self.score_variable_name = 'Score'

        box = gui.vBox(self.controlArea, "Gene name: ")
        gui.comboBox(box, self, 'gene', sendSelectedValue=True,
                     model=self.feature_model, callback=self._invalidate)

        box = gui.vBox(self.controlArea, "Column name: ")
        self.line_edit = gui.lineEdit(box, self, 'score_variable_name', callback=self._invalidate)
        self.line_edit.setPlaceholderText('Column name ...')

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

        self.sizeHint()

    def sizeHint(self):
        return QSize(320, 240)

    def handleNewSignals(self):
        self.Outputs.data.send(None)
        self._invalidate()

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        self.Warning.no_gene_id_attribute.clear()
        self.data = data

        if self.data:
            id_attribute = self.data.attributes.get(GENE_ID_ATTRIBUTE, None)

            if id_attribute is None:
                self._invalidate()
                self.Warning.no_gene_id_attribute()
                return

            for variable in self.data.domain.attributes:
                self.map_var_to_attribute[str(variable.name)] = str(variable.attributes.get(id_attribute, None))

    @Inputs.genes
    @check_sql_input
    def set_genes(self, genes):
        self.Warning.no_gene_id_column.clear()
        self.closeContext()
        self.genes = genes
        self.feature_model.set_domain(None)
        self.gene = None

        if self.genes:
            id_column = self.genes.attributes.get(GENE_ID_COLUMN, None)

            if id_column is None:
                self._invalidate()
                self.Warning.no_gene_id_column()
                return

            self.feature_model.set_domain(self.genes.domain)
            self.gene = self.genes.domain[id_column]

            self.openContext(genes)

    def __score_cells(self):
        score = np.zeros(len(self.data))
        matched = dict([(key, value) for key, value in self.map_var_to_attribute.items()
                        if value in self.marker_genes])

        self.Warning.no_genes.clear()
        self.Warning.some_genes.clear()

        if not matched:
            self.Warning.no_genes()
        else:
            if len(matched.keys()) < len(self.marker_genes):
                self.Warning.some_genes(len(self.marker_genes) - len(matched.keys()),
                                        len(self.marker_genes))

            values = self.data[:, matched.keys()].X
            score = np.nanmax(values, axis=1)

        return score

    def commit(self):

        if self.data is None:
            self.Outputs.data.send(None)
            return

        if self.genes and self.gene:
            self.marker_genes = set(self.genes.get_column_view(self.gene)[0])
            self.marker_genes = list(filter(None, self.marker_genes))
            score = self.__score_cells()

            score_var = ContinuousVariable(self.score_variable_name)
            domain = Domain(self.data.domain.attributes,
                            self.data.domain.class_vars,
                            self.data.domain.metas + (score_var,))

            table = self.data.transform(domain)
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

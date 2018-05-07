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

from orangecontrib.bioinformatics.ncbi.gene import NCBI_ID


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
        imperfect_match = Msg("Some input genes can't be matched")

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
        self._invalidate()

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        self.data = data
        if self.data:

            for variable in self.data.domain.attributes:
                self.map_var_to_attribute[str(variable.name)] = str(variable.attributes.get(NCBI_ID, None))

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

    def __score_cells(self):
        score = np.zeros(len(self.data))
        matched = dict([(key, value) for key, value in self.map_var_to_attribute.items()
                        if key in self.marker_genes or value in self.marker_genes])
        print(self.marker_genes)
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
        self.clear_messages()

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

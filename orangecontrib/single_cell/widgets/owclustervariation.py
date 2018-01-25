import numpy as np

from Orange.data import (ContinuousVariable, DiscreteVariable, StringVariable,
                         Domain, Table)
from Orange.widgets import widget, gui
from Orange.widgets.settings import (Setting, ContextSetting,
                                     DomainContextHandler)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.sql import check_sql_input


class OWDifferentialExpression(widget.OWWidget):
    name = "Cluster Variation"
    description = "Compare expressions in groups of samples"
    icon = "icons/ClusterVariation.svg"
    priority = 190

    settingsHandler = DomainContextHandler()
    groupby = ContextSetting(None)
    auto_apply = Setting(True)

    want_main_area = False

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        differential_expression = Output("Differential Expressions", Table)

    def __init__(self):
        super().__init__()

        self.data = None
        self.feature_model = DomainModel(valid_types=DiscreteVariable)

        box = gui.vBox(self.controlArea, "Group by")
        gui.comboBox(box, self, 'groupby', sendSelectedValue=True,
                     model=self.feature_model, callback=self._invalidate)

        self.apply_button = gui.auto_commit(
            self.controlArea, self, "auto_apply", "&Apply", box=False)

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        if self.feature_model:
            self.closeContext()
        self.data = data
        self.feature_model.set_domain(None)
        self.groupby = None
        if self.data:
            self.feature_model.set_domain(self.data.domain)
            if self.feature_model:
                self.groupby = self.feature_model[0]
                self.openContext(data)

    def handleNewSignals(self):
        self._invalidate()

    def commit(self):
        table = None
        if self.data and self.groupby:
            table = differential_expression(self.data, self.groupby)
        self.Outputs.differential_expression.send(table)

    def _invalidate(self):
        self.commit()

    def send_report(self):
        groupby = None
        method = 'Log2 Fold change'
        if self.data is not None:
            groupby = self.groupby
            if groupby in self.data.domain:
                groupby = self.data.domain[groupby]
        self.report_items((
            ("Group by", groupby),
            ("Method", method),
        ))


def differential_expression(data, groupby):
    group_col, sparse = data.get_column_view(groupby)
    groups = np.unique(group_col)
    X = np.zeros((data.X.shape[1], len(groups)))
    for gi, g in enumerate(groups):
        ind = group_col == g
        in_group = np.nanmean(data.X[ind], axis=0)
        out_group = np.nanmean(data.X[~ind], axis=0)
        fc = in_group / (out_group + 0.0001)
        logfc = np.log2(abs(fc) + 0.0001)
        X[:, gi] = logfc

    metavar = StringVariable('name')
    metas = [[f.name] for f in data.domain.attributes]
    domain = Domain([ContinuousVariable(groupby.repr_val(g)) for g in groups], metas=[metavar])
    return Table(domain, X, metas=metas)


def test():
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])

    w = OWDifferentialExpression()
    data = Table("iris")
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    app.exec_()


if __name__ == "__main__":
    test()

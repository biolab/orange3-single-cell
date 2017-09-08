import numpy as np

from Orange.data import Table, Domain, DiscreteVariable, ContinuousVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel


class OWClusterStatistics(widget.OWWidget):
    name = 'Cluster Statistics'
    description = ('Compute variant calling statistics for clusters')
    icon = 'icons/ClusterStatistics.svg'
    priority = 110

    inputs = [("Data", Table, 'set_data')]
    outputs = [("Data", Table)]

    want_main_area = False
    resizing_enabled = False

    settingsHandler = settings.DomainContextHandler(metas_in_res=True)
    selected_attr = settings.ContextSetting("")
    autocommit = settings.Setting(True)

    class Warning(widget.OWWidget.Warning):
        no_discrete_attributes = widget.Msg("No cluster variables in input data.")

    def __init__(self):
        self.data = None
        self.info = gui.label(self.controlArea, self,
                              "No data on input", box="Info")

        attrs_model = self.attrs_model = DomainModel(
            order=(DomainModel.CLASSES, DomainModel.METAS),
            valid_types=DiscreteVariable)
        combo_attrs = self.combo_attrs = gui.comboBox(
            self.controlArea, self, 'selected_attr',
            box="Cluster Column",
            callback=self.on_changed,
            sendSelectedValue=True)
        combo_attrs.setModel(attrs_model)

        gui.auto_commit(self.controlArea, self, 'autocommit', '&Apply')

    def set_data(self, data):
        self.closeContext()
        self.data = data
        self.Warning.clear()

        if self.data is None:
            self.attrs_model.clear()
            self.commit()
            self.info.setText("No data on input")
            return

        self.info.setText("%d samples, %d mutations" %
                          (len(data), len(data.domain.attributes)))

        self.attrs_model.set_domain(data.domain)
        if len(self.attrs_model) != 0:
            self.selected_attr = str(self.attrs_model[0])
        else:
            self.Warning.no_discrete_attributes()
            self.send("Data", None)
            return

        self.openContext(self.data.domain)
        self.on_changed()

    def on_changed(self):
        self.commit()

    def commit(self):
        if self.data is None or len(self.attrs_model) == 0 or \
                        self.selected_attr not in self.data.domain:
            self.send("Data", None)
            return

        cluster_var = self.data.domain[self.selected_attr]
        metas = [
            cluster_var,
            ContinuousVariable("size", number_of_decimals=0)
        ]
        atts = [ContinuousVariable("#%s" % a.name)
                for a in self.data.domain.attributes]
        for new_a, old_a in zip(atts, self.data.domain.attributes):
            new_a.attributes = old_a.attributes
        domain = Domain(atts, None, metas=metas)

        clusters = self.data.get_column_view(self.selected_attr)[0]
        X = self.data.X
        Z = np.array([np.nansum(X[clusters == v], axis=0)
                      for v in range(len(cluster_var.values))])

        cluster_counts = np.bincount(
            clusters.astype(int), minlength=len(cluster_var.values))

        M = np.column_stack((np.arange(len(cluster_var.values)),
                             cluster_counts))
        count_data = Table(domain, Z, metas=M)
        self.send("Data", count_data)


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication

    app = QApplication([])
    ow = OWClusterStatistics()

    table = Table("/Users/blaz/Desktop/sparse.xlsx")
    ow.set_data(table)

    ow.show()
    app.exec()
    ow.saveSettings()

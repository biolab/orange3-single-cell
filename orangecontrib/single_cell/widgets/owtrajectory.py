from itertools import chain
from xml.sax.saxutils import escape

import numpy as np
from scipy.stats import linregress
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import r2_score

from AnyQt.QtCore import Qt, QTimer, QPointF, QLineF
from AnyQt.QtGui import QColor, QPainterPath
from AnyQt.QtWidgets import QApplication, QFormLayout, QGraphicsPathItem, QGraphicsLineItem

import pyqtgraph as pg

import Orange
from Orange.data import Table, Domain, ContinuousVariable, DiscreteVariable
from Orange.data.sql.table import SqlTable, AUTO_DL_LIMIT
from Orange.preprocess.score import ReliefF, RReliefF

from Orange.widgets import gui, report, widget
from Orange.widgets.settings import (
    DomainContextHandler, Setting, ContextSetting, SettingProvider
)
from Orange.widgets.utils import get_variable_values_sorted
from Orange.widgets.utils.annotated_data import (
    create_annotated_table, create_groups_table, ANNOTATED_DATA_SIGNAL_NAME
)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.visualize.owscatterplotgraph import (
    OWScatterPlotGraph, OWScatterPlotBase, OWProjectionWidget
)
from Orange.widgets.visualize.utils.plotutils import AnchorItem
from Orange.widgets.visualize.utils import VizRankDialogAttrPair
from Orange.widgets.widget import AttributeList, Msg, Input, Output

from orangecontrib.single_cell.preprocess.pseudotimemst import Pseudotimemst

MAX_COMPONENTS = 20
MAX_CLUSTERS = 20
DIM_REDUCTIONS = [
    ("PCA", Orange.projection.PCA),
    #("UMAP", Orange.projection.UMAP),
]

class OWPseudotimeGraph(OWScatterPlotBase):
    jitter_size = Setting(0)
    jitter_continuous = Setting(False)

    def __init__(self, scatter_widget, parent):
        super().__init__(scatter_widget, parent)
        self.reg_line_item = None
        self.mst_coord = None
        self.proj_lines = None

    def get_mst(self):
        self.mst_coord = self.master.get_mst_data()
        return self.mst_coord

    def update_mst(self):
        pen = pg.mkPen(QColor(Qt.black), width=1, cosmetic=True)
        self._mst = self.get_mst()
        for a, b in self._mst:
            self._mst_line = QGraphicsLineItem()
            self._mst_line.setLine(a[0],a[1],b[0],b[1])
            self._mst_line.setPen(pen)
            self.plot_widget.addItem(self._mst_line)

    def get_proj_lines(self):
        dest, [x, y] = self.master.get_proj_lines_data()
        res = list()
        for i in range(len(dest)):
            res.append([dest[i], [x[i], y[i]]])
        return res

    def update_proj_lines(self):
        pen = pg.mkPen(QColor(Qt.black), width=.5, cosmetic=True)
        proj_lines = self.get_proj_lines()

        for a, b in proj_lines:
            self._mst_line = QGraphicsLineItem()
            self._mst_line.setLine(a[0], a[1], b[0], b[1])
            self._mst_line.setPen(pen)
            self.plot_widget.addItem(self._mst_line)

    def reset_graph(self):
        super().reset_graph()
        self.update_mst()
        self.update_proj_lines()

    def clear(self):
        super().clear()
        self.reg_line_item = None
        self.mst_coord = None
        self.proj_lines = None

    def set_axis_labels(self, axis, labels):
        axis = self.plot_widget.getAxis(axis)
        if labels:
            axis.setTicks([list(enumerate(labels))])
        else:
            axis.setTicks(None)

    def set_axis_title(self, axis, title):
        self.plot_widget.setLabel(axis=axis, text=title)

    def update_coordinates(self):
        super().update_coordinates()



class OWTrajectory(OWProjectionWidget):
    name = "Trajectory"
    description = "Placeholder"
    icon = "icons/Trajectory.svg"
    priority = 1000

    class Inputs:
        data = Input("Data", Orange.data.Table, default=True)
        data_subset = Input("Data Subset", Orange.data.Table)

    class Outputs:
        pseudotimes = Output("Pseudotimes for Cells", Orange.data.Table, default=True)

    settingsHandler = DomainContextHandler()

    auto_send_selection = Setting(True)
    auto_sample = Setting(True)
    tooltip_shows_all = Setting(True)

    n_components = Setting(2)
    n_clusters = Setting(10)
    dim_reduction = Setting("PCA")


    graph = SettingProvider(OWPseudotimeGraph)
    graph_name = "graph.plot_widget.plotItem"

    class Warning(widget.OWWidget.Warning):
        placeholder = Msg("Placeholder warning")

    class Information(widget.OWWidget.Information):
        missing_size = Msg(
            "Points with undefined '{}' are shown in smaller size")
        missing_shape = Msg(
            "Points with undefined '{}' are shown as crossed circles")


    def __init__(self):
        super().__init__()

        box = gui.vBox(self.mainArea, True, margin=0)
        self.graph = OWPseudotimeGraph(self, box)
        box.layout().addWidget(self.graph.plot_widget)

        self.subset_data = None
        self.subset_indices = None
        self.sql_data = None
        self.xy = None
        #self.data = None  # neaded?
        self.__timer = QTimer(self, interval=1200)
        self.__timer.timeout.connect(self.add_data)

        common_options = dict(
            labelWidth=50, orientation=Qt.Horizontal, sendSelectedValue=True,
            valueType=str, contentsLength=14
        )

        box = gui.vBox(self.controlArea, "Controls")
        dmod = DomainModel
        self.dom_model = DomainModel(dmod.MIXED, valid_types=dmod.PRIMITIVE)

        form = QFormLayout()

        gui.comboBox(
            box, self, "dim_reduction",
            callback=self.attr_changed,
            editable=False,
            items=[item for item, _ in DIM_REDUCTIONS],
            sendSelectedValue=True
        )
        gui.spin(
            box, self, "n_components", 2, MAX_COMPONENTS,
            callback=self.attr_changed,
            keyboardTracking=False,
        )

        gui.spin(
            box, self, "n_clusters", 2, MAX_CLUSTERS,
            callback=self.attr_changed,
            keyboardTracking=False,
        )
        form.addRow("Dimensionality reduction:", self.controls.dim_reduction)
        form.addRow("Num. of components:", self.controls.n_components)
        form.addRow("Num. of clusters:", self.controls.n_clusters)
        box.layout().addLayout(form)

        g = self.graph.gui

        g.point_properties_box(self.controlArea)

        box = g.create_gridbox(self.controlArea, True)
        g.add_widgets([
            g.PointSize,
            g.AlphaValue,
            #g.JitterSizeSlider,
            ], box
        )
        #g.add_widget(g.JitterNumericValues, box)

        box_plot_prop = g.plot_properties_box(self.controlArea)
        g.add_widgets([
            g.ShowGridLines,
            g.ToolTipShowsAll],
            box_plot_prop)

        self.controlArea.layout().addStretch(100)
        self.graph.box_zoom_select(self.controlArea)
        gui.auto_commit(self.controlArea, self, "auto_send_selection",
                        "Send Selection", "Send Automatically")

    @Inputs.data
    def set_data(self, data):
        #self.Information.sampled_sql.clear()

        if isinstance(data, SqlTable):
            if data.approx_len() < 5000:
                data = Table(data)
            else:
                self.Information.sampled_sql()
                self.sql_data = data
                data_sample = data.sample_time(0.8, no_cache=True)
                data_sample.download_data(2000, partial=True)
                data = Table(data_sample)
                self.sampling.setVisible(True)
                if self.auto_sample:
                    self.__timer.start()

        if data is not None and (len(data) == 0 or len(data.domain) == 0):
            data = None
        if self.data and data and self.data.checksum() == data.checksum():
            return
        self.closeContext()
        same_domain = (self.data and data and
                       data.domain.checksum() == self.data.domain.checksum())
        self.data = data

        if not same_domain:
            self.init_attr_values()
        self.openContext(self.data)
        self.attr_changed()


    def fit_transform(self):
        if self.data is None:
            return
        p = Pseudotimemst(n_components=self.n_components, n_clusters=self.n_clusters,
                          projection=[ext for name, ext in DIM_REDUCTIONS if name == self.dim_reduction][0])
        p.fit_transform(self.data)
        self.p = p
        self.xy = [np.array([row[0] for row in p.transformed_data]), np.array([row[1] for row in p.transformed_data])]

    def attr_changed(self):
        self.fit_transform()
        self.graph.reset_graph()
        self.commit()

        for axis, var in (("bottom", "PCA component1"), ("left", "PCA component 2")):
            self.graph.set_axis_title(axis, None)
            self.graph.set_axis_labels(axis, None)


    def get_mst_data(self):
        return self.p.mst_coordinates

    def get_proj_lines_data(self):
        return self.p.proj, self.xy

    def get_coordinates_data(self):
        return self.xy[0], self.xy[1]

    @Inputs.data_subset
    def set_subset_data(self, subset_data):
        self.warning()
        pass

    def add_data(self, time=0.4):
        if self.data and len(self.data) > 2000:
            return self.__timer.stop()
        data_sample = self.sql_data.sample_time(time, no_cache=True)
        if data_sample:
            data_sample.download_data(2000, partial=True)
            data = Table(data_sample)
            self.data = Table.concatenate((self.data, data), axis=0)
            self.handleNewSignals()


    def handle_new_signals(self):
        if self.data is None:
            return
        self.attr_changed()
        self.graph.set_
        self.cb_class_density.setEnabled(self.can_draw_density())
        self.commit()


    def send_pseudotimes(self):
        pt = self.p.pseudotime
        table = Table.from_numpy(
            Domain(
                [ContinuousVariable.make("pseudotime")],
                self.data.domain.class_vars,
                self.data.domain.metas
            ),
            np.array([pt]).T,
            self.data.Y,
            self.data.metas
        )
        self.Outputs.pseudotimes.send(table)

    def commit(self):
        self.send_pseudotimes()


if __name__ == '__main__':
    a = QApplication([])

    ow = OWTrajectory()

    tab = Orange.data.Table("brown-selected")
    ow.set_data(tab)
    ow.show()
    ow.raise_()

    rval = a.exec()

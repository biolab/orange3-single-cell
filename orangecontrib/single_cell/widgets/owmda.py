from builtins import property

from AnyQt.QtWidgets import QFormLayout
from AnyQt.QtGui import QColor
from AnyQt.QtCore import Qt

import numpy as np
import pyqtgraph as pg
from pyqtgraph.graphicsItems import LegendItem
from pyqtgraph import functions as fn
from scipy.stats import multivariate_normal as mvn

from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.data import Table, Domain, ContinuousVariable, DiscreteVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.widget import Input, Output
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.utils.annotated_data import add_columns

from orangecontrib.single_cell.preprocess.alignment import SeuratAlignmentModel

# Maximum number of components and metagenes
MAX_COMPONENTS = 50
MAX_COMPONENTS_DEFAULT = 50
MAX_GENES = 100

SCORINGS = {
    "Pearson correlation": "pearson",
    "Spearman correlation": "spearman",
    "Biweights midcorrelation": "bicor",
}
SCORINGS_NAMES = (
    "Pearson correlation",
    "Spearman correlation",
    "Biweights midcorrelation",
)


class MyItemSample(LegendItem.ItemSample):
    def paint(self, p, *args):
        opts = self.item.opts
        opts['pen'].setWidth(2)
        p.setPen(fn.mkPen(opts['pen']))
        p.drawLine(0, 10, 9, 10)


def smooth_correlations(M, offset=2):
    """
    Smoothed correlations; Convolve with a small Gaussian bump.
    :param M: model.shared_correlations
    :param offset: First few values to leave unchanged
    :return:
    """
    g = mvn.pdf(x=np.linspace(-offset, offset, 2 * offset + 1))
    Z = np.zeros(M.shape)
    for i, m in enumerate(M):
        Z[i, :offset] = m[:offset]
        Z[i, offset:] = np.convolve(m, g, mode="same")[offset:]
    return Z


def interpolate_nans(A):
    ok = ~np.isnan(A)
    xp = ok.ravel().nonzero()[0]
    fp = A[ok]
    x = np.isnan(A).ravel().nonzero()[0]
    A[np.isnan(A)] = np.interp(x, xp, fp)
    return A



class OWAlignDatasets(widget.OWWidget):
    name = "Align Datasets"
    description = "Alignment of multiple datasets with a scree-diagram correlation visualization."
    icon = "icons/AlignDatasets.svg"
    priority = 3050

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        transformed_data = Output("Transformed Data", Table)
        genes_components = Output("Genes per n. Components", Table)

    settingsHandler = DomainContextHandler()
    auto_update = ContextSetting(True)
    auto_commit = Setting(True)
    axis_labels = ContextSetting(10)
    source_id = ContextSetting(None)
    ncomponents = ContextSetting(20)
    ngenes = ContextSetting(30)
    scoring = ContextSetting(SCORINGS_NAMES[0])
    quantile_normalization = ContextSetting(False)
    quantile_normalization_perc = ContextSetting(2.5)
    dynamic_time_warping = ContextSetting(False)

    graph_name = "plot.plotItem"

    class Error(widget.OWWidget.Error):
        no_features = widget.Msg("At least 1 feature is required")
        no_instances = widget.Msg("At least 2 data instances are required")
        sparse_data = widget.Msg("Sparse data is not supported")

    def __init__(self):
        super().__init__()
        self.data = None
        self.source_id = None
        self._mas = None
        self._Ws = None
        self._transformed = None
        self._components = None
        self._use_genes = None
        self._shared_correlations = None
        self._transformed_table = None
        self._line = False
        self._feature_model = DomainModel(valid_types=DiscreteVariable)
        self._feature_model.set_domain(None)
        self._init_mas()
        self._legend = None
        form = QFormLayout(
            labelAlignment=Qt.AlignLeft,
            formAlignment=Qt.AlignLeft,
            fieldGrowthPolicy=QFormLayout.AllNonFixedFieldsGrow,
            verticalSpacing=10
        )
        # Data source indicator
        box = gui.vBox(self.controlArea, "Data source indicator")

        gui.comboBox(
            box, self, "source_id", sendSelectedValue=True,
            callback=self._update_combo_source_id,
            model=self._feature_model,
        )

        # Canonical correlation analysis
        box = gui.vBox(self.controlArea, "Canonical correlation analysis")
        gui.spin(
            box, self, "ncomponents", 1, MAX_COMPONENTS,
            callback=self._update_selection_component_spin,
            keyboardTracking=False,
            label="Num. of components"
        )

        # Shared genes
        box = gui.vBox(self.controlArea, "Shared genes")
        gui.spin(
            box, self, "ngenes", 1, MAX_GENES,
            callback=self._update_ngenes_spin,
            keyboardTracking=False,
        )
        form.addRow(
            "Num. of genes",
            self.controls.ngenes
        )

        gui.comboBox(
            box, self, "scoring",
            callback=self._update_scoring_combo,
            items=SCORINGS_NAMES, sendSelectedValue=True,
            editable=False,
        )
        form.addRow(
            "Scoring:",
            self.controls.scoring
        )

        box.layout().addLayout(form)

        # Post-processing
        box = gui.vBox(self.controlArea, "Post-processing")
        gui.doubleSpin(
            box, self, "quantile_normalization_perc", minv=0, maxv=49, step=5e-1,
            callback=self._update_quantile_normalization,
            checkCallback=self._update_quantile_normalization,
            controlWidth=80, alignment=Qt.AlignRight,
            label="Quantile normalization", checked="quantile_normalization",
        )
        self.controls.quantile_normalization_perc.setSuffix("%")
        gui.checkBox(
            box, self, "dynamic_time_warping",
            callback=self._update_dynamic_time_warping,
            label="Dynamic time warping"
        )

        self.controlArea.layout().addStretch()

        gui.auto_commit(self.controlArea, self, "auto_commit", "Apply",
                        callback=self._invalidate_selection(),
                        checkbox_label="Apply automatically")

        self.plot = pg.PlotWidget(background="w")

        axis = self.plot.getAxis("bottom")
        axis.setLabel("Correlation components")
        axis = self.plot.getAxis("left")
        axis.setLabel("Correlation strength")
        self.plot_horlabels = []
        self.plot_horlines = []

        self.plot.getViewBox().setMenuEnabled(False)
        self.plot.getViewBox().setMouseEnabled(False, False)
        self.plot.showGrid(True, True, alpha=0.5)
        self.plot.setRange(xRange=(0.0, 1.0), yRange=(0.0, 1.0))

        self.mainArea.layout().addWidget(self.plot)

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        # if self._feature_model:
        #    self.closeContext()
        self.clear_messages()
        self.clear()
        self.information()
        self.clear_outputs()
        self._feature_model.set_domain(None)

        if data:
            self._feature_model.set_domain(data.domain)
            if self._feature_model:
                # self.openContext(data)
                if self.source_id is None or self.source_id == '':
                    self.source_id = self._feature_model[0]
            if len(data.domain.attributes) == 0:
                self.Error.no_features()
                return
            if len(data) == 0:
                self.Error.no_instances()
                return
            self.data = data

            global MAX_COMPONENTS
            if len(data.domain.attributes) < MAX_COMPONENTS_DEFAULT:
                MAX_COMPONENTS = len(data.domain.attributes)
            else:
                MAX_COMPONENTS = MAX_COMPONENTS_DEFAULT

        self.fit()

    def fit(self):
        if self.data is None:
            return

        self._init_mas()
        X = self.data.X
        y = self.data.get_column_view(self.source_id)[0]

        self._Ws = self._mas.fit(X, y)
        shared_correlations = self._mas.shared_correlations
        self._shared_correlations = np.array([interpolate_nans(x) for x in shared_correlations])
        self._use_genes = self._mas.use_genes

        self._setup_plot()
        if self.auto_commit:
            self.commit()

    def clear(self):
        self.data = None
        self.source_id = None
        self._mas = None
        self._Ws = None
        self._transformed = None
        self._transformed_table = None
        self._components = None
        self._use_genes = None
        self._shared_correlations = None
        self._line = False
        # self._legend = None
        self._feature_model.set_domain(None)
        self.plot_horlabels = []
        self.plot_horlines = []
        self.plot.clear()

    def clear_outputs(self):
        self.Outputs.transformed_data.send(None)
        self.Outputs.genes_components.send(None)

    def _init_mas(self):
        self._mas = SeuratAlignmentModel(
            n_components=MAX_COMPONENTS,
            n_metagenes=self.ngenes,
            gene_scoring=SCORINGS[self.scoring],
        )

    def get_model(self):
        if self.data is None:
            return

        self.fit()
        self._setup_plot()
        self.commit()

    def _setup_plot(self):
        self.plot.clear()
        if self._mas is None:
            return

        shared_correlations = self._shared_correlations
        p = MAX_COMPONENTS

        # Colors chosen based on: http://colorbrewer2.org/?type=qualitative&scheme=Set1&n=9
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

        if self._legend is not None:
            self._legend.scene().removeItem(self._legend)
        self._legend = self.plot.addLegend(offset=(-1, 1))
        # correlation lines
        smoothed_correlations = smooth_correlations(shared_correlations, offset=2)
        plotitem = dict()
        for i, corr in enumerate(smoothed_correlations):
            plotitem[i] = self.plot.plot(np.arange(p), corr, pen=pg.mkPen(QColor(colors[i]), width=2),
                                            antialias=True) # name=self.source_id.values[i]
        #self.plot.plotItem.legend.addItem(3, "maximum value")

        for i in range(len(plotitem)):
            self._legend.addItem(MyItemSample(pg.ScatterPlotItem(pen=colors[i])), self.source_id.values[i])

        # vertical movable line
        cutpos = self.ncomponents - 1
        self._line = pg.InfiniteLine(
            angle=90, pos=cutpos, movable=True, bounds=(0, p - 1))
        self._line.setCursor(Qt.SizeHorCursor)
        self._line.setPen(pg.mkPen(QColor(Qt.black), width=2))
        self._line.sigPositionChanged.connect(self._on_cut_changed)
        self.plot.addItem(self._line)

        # horizontal lines
        self.plot_horlines = tuple(
            pg.PlotCurveItem(pen=pg.mkPen(QColor(colors[i]), style=Qt.DashLine)) for i in
            range(len(shared_correlations))
        )
        self.plot_horlabels = tuple(
            pg.TextItem(color=QColor('k'), anchor=(0, 1)) for _ in range(len(shared_correlations))
        )

        for item in self.plot_horlabels + self.plot_horlines:
            self.plot.addItem(item)
        self._set_horline_pos()

        #self.plot.setRange(xRange=(0.0, p - 1), yRange=(0.0, 1.0))
        self.plot.setXRange(0.0, p - 1, padding=0)
        self.plot.setYRange(0.0, 1.0, padding=0)
        self._update_axis()

    def _set_horline_pos(self):
        cutidx = self.ncomponents - 1
        for line, label, curve in zip(self.plot_horlines, self.plot_horlabels,
                                      self._shared_correlations):
            y = curve[cutidx]
            line.setData([-1, cutidx], 2 * [y])
            label.setPos(cutidx, y)
            label.setPlainText("{:.3f}".format(y))

    def _on_cut_changed(self, line):
        # cut changed by means of a cut line over the scree plot.
        value = int(round(line.value()))
        components = value + 1

        if not (self.ncomponents == 0 and
                components == len(self._components)):
            self.ncomponents = components

        self._line.setValue(value)
        self._set_horline_pos()
        self.commit()

    def _update_selection_component_spin(self):
        # cut changed by "ncomponents" spin.
        if self._mas is None:
            self._invalidate_selection()

        if np.floor(self._line.value()) + 1 != self.ncomponents:
            self._line.setValue(self.ncomponents - 1)

        self.commit()

    def _invalidate_selection(self):
        if self.data is not None:
            self._transformed = None
            self.commit()

    def _update_ngenes_spin(self):
        self.fit()
        self._invalidate_selection()

    def _update_scoring_combo(self):
        self.fit()
        self._invalidate_selection()

    def _update_dynamic_time_warping(self):
        self._invalidate_selection()

    def _update_quantile_normalization(self):
        self._invalidate_selection()

    def _update_combo_source_id(self):
        self.fit()
        self._invalidate_selection()

    def _update_axis(self):
        p = MAX_COMPONENTS
        axis = self.plot.getAxis("bottom")
        d = max((p - 1) // (self.axis_labels - 1), 1)
        axis.setTicks([[(i, str(i + 1)) for i in range(0, p, d)]])

    def commit(self):
        transformed_table = meta_genes = genes_componenets = None
        if self._mas is not None:
            if self._transformed is None:
                # Compute the full transform (MAX_COMPONENTS components) only once.
                X = self.data.X
                y = self.data.get_column_view(self.source_id)[0]
                self._transformed = self._mas.transform(X, y, normalize=self.quantile_normalization,
                                                        quantile=self.quantile_normalization_perc,
                                                        dtw=self.dynamic_time_warping)

                attributes = tuple(ContinuousVariable.make("CCA{}".format(x + 1)) for x in
                                   range(MAX_COMPONENTS))
                dom = Domain(
                    attributes,
                    self.data.domain.class_vars,
                    self.data.domain.metas
                )

                # Meta-genes
                meta_genes = self.data.transform(dom)
                genes_components = np.zeros((self.data.X.shape[1], MAX_COMPONENTS))
                for key, genes in self._mas.use_genes.items():
                    for gene in genes:
                        genes_components[gene - 1, key] = genes.index(gene) + 1
                genes_components[genes_components == 0] = np.NaN
                meta_genes.X = genes_components
                self.meta_genes = Table.from_numpy(Domain(attributes), genes_components)

                # Transformed data
                transformed = self._transformed
                new_domain = add_columns(self.data.domain, attributes=attributes)
                transformed_table_temp = self.data.transform(new_domain)
                transformed_table_temp.X[:, -MAX_COMPONENTS:] = transformed
                self.transformed_table = Table.from_table(dom, transformed_table_temp)

            ncomponents_attributes = tuple(ContinuousVariable.make("CCA{}".format(x + 1)) for x in
                      range(self.ncomponents))
            ncomponents_domain = Domain(
                    ncomponents_attributes,
                    self.data.domain.class_vars,
                    self.data.domain.metas
                )

            meta_genes = self.meta_genes.transform(Domain(ncomponents_attributes))
            transformed_table = self.transformed_table.transform(ncomponents_domain)

        self.Outputs.transformed_data.send(transformed_table)
        self.Outputs.genes_components.send(meta_genes)

    def send_report(self):
        if self.data is None:
            return
        self.report_items((
            ("Source ID", self.source_id),
            ("Selected num. of components", self.ncomponents),
            ("Selected num. of genes", self.ngenes),
            ("Scoring", self.scoring),
            ("Quantile normalization", True if self.quantile_normalization else "False"),
            ("Quantile normalization percentage",
             self.quantile_normalization_perc if self.quantile_normalization else False),
            ("Dynamic time warping", True if self.dynamic_time_warping else "False")
        ))
        self.report_plot()

    """
    @classmethod
    def migrate_settings(cls, settings, version):
        if "variance_covered" in settings:
            # Due to the error in gh-1896 the variance_covered was persisted
            # as a NaN value, causing a TypeError in the widgets `__init__`.
            vc = settings["variance_covered"]
            if isinstance(vc, numbers.Real):
                if np.isfinite(vc):
                    vc = int(vc)
                else:
                    vc = 100
                settings["variance_covered"] = vc
        if settings.get("ncomponents", 0) > MAX_COMPONENTS:
            settings["ncomponents"] = MAX_COMPONENTS
    """


def main():
    import gc
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    w = OWAlignDatasets()
    in_file = "../tutorials/Showcase-SampleAlignment-data/data_kang2018.tab.gz"
    data = Table(in_file)
    w.set_data(data)
    w.show()
    w.raise_()
    rval = w.exec()
    w.deleteLater()
    del w
    app.processEvents()
    gc.collect()
    return rval


if __name__ == "__main__":
    main()

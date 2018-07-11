from AnyQt.QtWidgets import QFormLayout
from AnyQt.QtGui import QColor
from AnyQt.QtCore import Qt

import numpy as np
import pyqtgraph as pg

from Orange.data import Table, Domain, ContinuousVariable, DiscreteVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.widget import Input, Output

from orangecontrib.single_cell.preprocess.alignment import SeuratAlignmentModel

# Maximum number of COMPONENTS AND METAGENES
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


class OWMDA(widget.OWWidget):
    name = "MDA"
    description = "Multi-data alignment with a scree-diagram."
    icon = "icons/MDA.svg"
    priority = 3050

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        transformed_data = Output("Transformed data", Table)
        genes = Output("Meta-genes", np.ndarray)
        genes_components = Output("Genes per n. components", Table)

    # settingsHandler = settings.DomainContextHandler()
    auto_update = settings.Setting(True)
    auto_commit = settings.Setting(True)
    axis_labels = settings.Setting(10)
    source_id = settings.Setting('')
    ncomponents = settings.Setting(20)
    ngenes = settings.Setting(30)
    scoring = settings.Setting(SCORINGS_NAMES[0])
    quantile_normalization = settings.Setting(False)
    quantile_normalization_perc = settings.Setting(2.5)
    dynamic_time_warping = settings.Setting(False)

    graph_name = "plot.plotItem"

    class Warning(widget.OWWidget.Warning):
        trivial_components = widget.Msg(
            "All components of the PCA are trivial (explain 0 variance). "
            "Input data is constant (or near constant).")

    class Error(widget.OWWidget.Error):
        no_features = widget.Msg("At least 1 feature is required")
        no_instances = widget.Msg("At least 2 data instances are required")
        sparse_data = widget.Msg("Sparse data is not supported")

    def __init__(self):
        super().__init__()
        self.data = None
        self._mas = None
        self._Ws = None
        self._transformed = None
        self._components = None
        self._use_genes = None
        self._shared_correlations = None
        self._transformed_table = None
        self._line = False
        self._feature_model = DomainModel(valid_types=DiscreteVariable)
        self._init_mas()

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
    def set_data(self, data):
        self.clear_messages()
        self.clear()
        self.information()
        self.clear_outputs()
        self._feature_model.set_domain(None)

        if data:
            self._feature_model.set_domain(data.domain)
            if self._feature_model:
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

    def fit(self):
        if self.data is None:
            return

        self._init_mas()
        X = self.data.X
        y = self.data.get_column_view(self.source_id)[0]

        self._Ws = self._mas.fit(X, y)
        self._shared_correlations = self._mas.shared_correlations
        self._use_genes = self._mas.use_genes

        self._setup_plot()

    def clear(self):
        self._mas = None
        self._Ws = None
        self._transformed = None
        self._transformed_table = None
        self._components = None
        self._use_genes = None
        self._shared_correlations = None
        self._line = False
        self.plot_horlabels = []
        self.plot_horlines = []
        self.plot.clear()

    def clear_outputs(self):
        self.Outputs.transformed_data.send(None)
        self.Outputs.genes.send(None)
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
        self._transformed = None
        self.commit()

    def _setup_plot(self):
        self.plot.clear()
        if self._mas is None:
            return

        shared_correlations = self._shared_correlations
        p = MAX_COMPONENTS

        # Colors chosen based on: http://colorbrewer2.org/?type=qualitative&scheme=Set1&n=9
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

        # correlation lines
        for i, corr in enumerate(shared_correlations):
            # Smoothening - not working, needs fix for NaN values in shared_correlations
            # f = interp1d(np.arange(p), corr, bounds_error=False, fill_value='extrapolate', kind='cubic')
            # xnew = np.linspace(0, p-1, num=3*p, endpoint=True)
            self.plot.plot(np.arange(p), corr,
                           pen=pg.mkPen(QColor(colors[i]), width=2),
                           antialias=True,
                           name="Correlation")

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
            pg.TextItem(color=QColor('k'), anchor=(1, 1)) for i in range(len(shared_correlations))
        )

        for item in self.plot_horlabels + self.plot_horlines:
            self.plot.addItem(item)
        self._set_horline_pos()

        self.plot.setRange(xRange=(0.0, p - 1), yRange=(0.0, 1.0))
        self._update_axis()
        return

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
            self._transformed_table = None
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
        transformed_table = genes = genes_componenets = None
        if self._mas is not None:
            if self._transformed_table is None:
                # Compute the full transform (MAX_COMPONENTS components) only once.
                X = self.data.X
                y = self.data.get_column_view(self.source_id)[0]
                self._transformed = self._mas.transform(X, y, normalize=self.quantile_normalization,
                                                        quantile=self.quantile_normalization_perc,
                                                        dtw=self.dynamic_time_warping)
            transformed = self._transformed[:, :self.ncomponents]

            # Meta genes for selected n components
            genes = np.array(self.data.domain.attributes)[self._use_genes[self.ncomponents - 1]]

            # Meta genes for all up to MAX n components
            genes_componenets = Table.from_numpy(
                Domain([DiscreteVariable.make("{}".format(x + 1)) for x in range(MAX_COMPONENTS)]),
                np.array([x for _, x in self._mas.use_genes.items()]).T
                )
            genes_componenets.domain.attributes = Domain(
                [DiscreteVariable.make("{}".format(x + 1)) for x in range(MAX_COMPONENTS)])

            # Transformed data
            transformed_table = Table.from_numpy(
                Domain(
                    [ContinuousVariable.make("CC{}".format(x)) for x in
                     range(self.ncomponents)],
                    self.data.domain.class_vars,
                    self.data.domain.metas
                ),
                transformed,
                self.data.Y,
                self.data.metas,
                self.data.W
            )
            transformed_table.domain.attributes = [ContinuousVariable.make("CC{}".format(x)) for x in
                                                   range(self.ncomponents)]
            self._transformed_table = transformed_table

        self.Outputs.transformed_data.send(transformed_table)
        self.Outputs.genes.send(genes)
        self.Outputs.genes_components.send(genes_componenets)

    def send_report(self):
        if self.data is None:
            return
        self.report_items((
            ("Source ID", self.source_id),
            ("Selected num. of components", self.ncomponents),
            ("Selected num. of genes", self.ngenes),
            ("Scoring", self.scoring),
            ("Quantile normalization", self.quantile_normalization),
            ("Quantile normalization percentage", self.quantile_normalization_perc),
            ("Dynamic time warping", self.dynamic_time_warping)
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
    w = OWMDA()
    in_file = "../tutorials/Showcase-SampleAlignment-data/data_kang2018.tab.gz"
    data = Table(in_file)
    data.Y[-2000:] = 2
    w.set_data(data)
    w.show()
    w.get_model()
    w.raise_()
    rval = w.exec()
    w.deleteLater()
    del w
    app.processEvents()
    gc.collect()
    return rval


if __name__ == "__main__":
    main()

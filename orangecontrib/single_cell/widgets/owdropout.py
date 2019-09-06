from collections import namedtuple
import numpy as np

from AnyQt.QtCore import Qt, QSize, QRectF, QPointF, pyqtSignal as Signal
from AnyQt.QtGui import QColor, QMouseEvent
from AnyQt.QtWidgets import QHBoxLayout, QVBoxLayout, QToolTip, \
    QGraphicsSceneMouseEvent

import pyqtgraph as pg

from Orange.data import Table
from Orange.widgets import gui, report
from Orange.widgets.settings import Setting, ContextSetting, \
    DomainContextHandler
from Orange.widgets.visualize.utils.plotutils import MouseEventDelegate
from Orange.widgets.widget import OWWidget, Input, Output, Msg

from orangecontrib.single_cell.preprocess.scpreprocess import \
    DropoutGeneSelection

ENTREZ_ID = "Entrez ID"
DropoutResults = namedtuple("DropoutResults",
                            ["zero_rate", "mean_expr", "decay",
                             "x_offset", "y_offset", "threshold"])


class States:
    WAITING, ON_CURVE, HOLDING_CURVE, MOVING_CURVE = range(4)


class DropoutGraph(pg.PlotWidget):
    pg.setConfigOption("foreground", "k")

    curve_moved = Signal(float, float)
    CURVE_PEN = pg.mkPen(color=QColor(Qt.darkCyan), width=4)
    MOVING_CURVE_PEN = pg.mkPen(color=QColor(179, 215, 255), width=4)

    def __init__(self, parent):
        super().__init__(parent, background="w")
        self.__dots = None  # type: pg.ScatterPlotItem
        self.__markers = None  # type: pg.ScatterPlotItem
        self.__curve = None  # type: pg.PlotCurveItem
        self.__decay = None  # type: float
        self.__x_offset = None  # type: float
        self.__y_offset = None  # type: float
        self._initial_x = None  # type: float
        self._initial_y = None  # type: float
        self._state = States.WAITING
        self._delegate = MouseEventDelegate(self.help_event, self.cursor_event)
        self.setMouseEnabled(False, False)
        self.hideButtons()
        self.plotItem.setContentsMargins(0, 20, 20, 0)
        self.setLabel("bottom", "Mean log2 nonzero expression")
        self.setLabel("left", "Frequency of zero expression")
        self.scene().installEventFilter(self._delegate)

    def set_data(self, results: DropoutResults, data: Table, genes: Table):
        self.__plot_dots(results.mean_expr, results.zero_rate, data)
        self.update_markers(data, genes)
        self.update_curve(results)
        self.__set_range(results.threshold, results.mean_expr)

    def update_markers(self, data: Table, genes: Table):
        self.removeItem(self.__markers)
        if data is None or self.__dots is None or genes is None:
            return
        x, y = self.__dots.getData()
        self.__plot_markers(x, y, data, genes)

    def update_curve(self, results: DropoutResults):
        self.__set_curve_params(results)
        self.__plot_curve(results)

    def __plot_dots(self, x: np.ndarray, y: np.ndarray, data: Table):
        data = list(zip(data.domain.attributes,
                        (data.X > 0).sum(axis=0),
                        np.full_like(y, len(data))))
        self.__dots = pg.ScatterPlotItem(x=x, y=y, size=5, data=data)
        self.addItem(self.__dots)

    def __plot_markers(self, x: np.ndarray, y: np.ndarray,
                       data: Table, markers: Table):
        col = markers.get_column_view(ENTREZ_ID)[0]
        mask = [str(a.attributes.get(ENTREZ_ID, None)) in col
                for a in data.domain.attributes]
        self.__markers = pg.ScatterPlotItem(
            x=x[mask], y=y[mask], size=7,
            brush=pg.mkBrush(color=QColor(Qt.magenta)))
        self.addItem(self.__markers)

    def __set_curve_params(self, results: DropoutResults):
        self.__decay = results.decay
        self.__x_offset = results.x_offset
        self.__y_offset = results.y_offset

    def __plot_curve(self, results: DropoutResults):
        self.removeItem(self.__curve)
        xmin, xmax = self.__get_xlim(results.threshold, results.mean_expr)
        x = np.arange(xmin, xmax + 0.01, 0.01)
        y = np.exp(-results.decay * (x - results.x_offset)) + results.y_offset
        pen = self.MOVING_CURVE_PEN if self._state == States.MOVING_CURVE \
            else self.CURVE_PEN
        self.__curve = pg.PlotCurveItem(
            x=x, y=y, fillLevel=1, pen=pen,
            brush=pg.mkBrush(color=QColor(0, 250, 0, 50)),
            antialias=True)
        self.addItem(self.__curve)

    def __set_range(self, threshold: float, x: np.ndarray):
        xmin, xmax = self.__get_xlim(threshold, x)
        rect = QRectF(xmin, 0, xmin + xmax, 1)
        self.setRange(rect, padding=0)

    def help_event(self, ev):
        if self.__dots is None:
            return False
        dot = self._dotAt(self.__dots.mapFromScene(ev.scenePos()))
        if dot is not None and dot.data() is not None:
            var, p, n = dot.data()
            q = round(p / n * 100, 1)
            text = f"{var.name}\nExpressed in {int(p)}/{int(n)} cells ({q}%)"
            QToolTip.showText(ev.screenPos(), text, widget=self)
            return True
        else:
            return False

    def _dotAt(self, pos: QPointF):
        pw, ph = self.__dots.pixelWidth(), self.__dots.pixelHeight()
        for s in self.__dots.points():
            sx, sy = s.pos().x(), s.pos().y()
            s2x = s2y = s.size()
            if self.__dots.opts['pxMode']:
                s2x *= pw
                s2y *= ph
            if sx + s2x > pos.x() > sx - s2x and sy + s2y > pos.y() > sy - s2y:
                return s
        return None

    def cursor_event(self, ev: QGraphicsSceneMouseEvent):
        if self.__curve is None:
            return False
        if self._state == States.HOLDING_CURVE:
            return False

        pos = self.__curve.mapFromScene(ev.scenePos())
        if self._state == States.MOVING_CURVE:
            self.curve_moved.emit(pos.x() - self._initial_x,
                                  pos.y() - self._initial_y)
            self._initial_x = pos.x()
            self._initial_y = pos.y()
            return False

        if self._on_curve(pos.x(), pos.y()):
            self._initial_x = pos.x()
            self._initial_y = pos.y()
            if self._state == States.WAITING:
                self._state = States.ON_CURVE
                self.getViewBox().setCursor(Qt.SizeAllCursor)
        else:
            self.getViewBox().setCursor(Qt.ArrowCursor)
            if self._state == States.ON_CURVE:
                self._state = States.WAITING
        return False

    def _on_curve(self, x: float, y: float):
        if self.__curve is None:
            return False
        return abs(DropoutGeneSelection.y(
            x, self.__decay, self.__x_offset, self.__y_offset) - y) < 0.01

    def mousePressEvent(self, ev: QMouseEvent):
        super().mousePressEvent(ev)
        if self._state == States.ON_CURVE:
            self._state = States.HOLDING_CURVE
            self.__curve.setPen(self.MOVING_CURVE_PEN)

    def mouseMoveEvent(self, ev: QMouseEvent):
        super().mouseMoveEvent(ev)
        if self._state == States.HOLDING_CURVE:
            self._state = States.MOVING_CURVE

    def mouseReleaseEvent(self, ev: QMouseEvent):
        super().mouseReleaseEvent(ev)
        if self._state in (States.HOLDING_CURVE, States.MOVING_CURVE):
            self._state = States.WAITING
            self.__curve.setPen(self.CURVE_PEN)

    def clear_all(self):
        self.clear()
        self.__dots = None
        self.__curve = None
        self.__decay = None
        self.__x_offset = None
        self.__y_offset = None

    @staticmethod
    def __get_xlim(threshold: float, x: np.ndarray):
        if not len(x):
            return 0, 0
        xmin = 0 if threshold == 0 else np.log2(threshold)
        return xmin, np.ceil(np.nanmax(x))


class FilterType:
    ByNumber, ByEquation = range(2)


class OWDropout(OWWidget):
    name = "Dropout Gene Selection"
    description = "Dropout-based gene selection"
    icon = 'icons/Dropout.svg'
    priority = 205

    class Inputs:
        data = Input("Data", Table, default=True)
        genes = Input("Genes", Table)

    class Outputs:
        data = Output("Data", Table)

    class Warning(OWWidget.Warning):
        less_selected = Msg("Cannot select more than {} genes.")
        missing_entrez_id = Msg("'Entred ID' is missing in Genes table.")

    settingsHandler = DomainContextHandler()
    filter_type = Setting(FilterType.ByNumber)
    n_genes = ContextSetting(1000)
    x_offset = ContextSetting(5)
    y_offset = ContextSetting(0.02)
    decay = ContextSetting(1)
    auto_commit = Setting(True)

    graph_name = "graph.plotItem"

    def __init__(self):
        super().__init__()
        self.data = None  # type: Table
        self.genes = None  # type: Table
        self.zero_rate = None  # type: np.ndarray
        self.mean_expr = None  # type: np.ndarray
        self.selected = None  # type: np.ndarray
        self.setup_gui()

    def setup_gui(self):
        self._add_graph()
        self._add_controls()

    def _add_graph(self):
        box = gui.vBox(self.mainArea, True, margin=0)
        self.graph = DropoutGraph(self)
        self.graph.curve_moved.connect(self.__manual_move)
        box.layout().addWidget(self.graph)

    def _add_controls(self):
        info_box = gui.widgetBox(self.controlArea, "Info")
        self.info_label = gui.widgetLabel(info_box, "")

        filter_box = gui.radioButtons(
            self.controlArea, self, "filter_type", box="Filter",
            orientation=QVBoxLayout(), callback=self.__filter_type_changed)

        genes_layout = QHBoxLayout()
        formula_layout = QVBoxLayout()
        genes_layout.addWidget(gui.appendRadioButton(
            filter_box, "Number of genes:", addToLayout=False))
        genes_layout.addWidget(gui.spin(
            filter_box, self, "n_genes", 0, 100000, addToLayout=False,
            callback=self.__param_changed))
        formula_layout.addWidget(gui.appendRadioButton(
            filter_box, "Apply exp(-a(x-b))+c", addToLayout=False))
        filter_box.layout().addLayout(genes_layout)
        filter_box.layout().addLayout(formula_layout)

        gui.separator(filter_box, height=1)
        coef_box = gui.hBox(filter_box, False)
        gui.separator(coef_box, width=15)
        common = dict(orientation=Qt.Horizontal,
                      callback=self.__param_changed,
                      alignment=Qt.AlignRight, controlWidth=60)
        gui.doubleSpin(
            coef_box, self, "decay", 0.0, 10.0, 0.01, label="a:", **common)
        gui.doubleSpin(
            coef_box, self, "x_offset", 0.0, 10.0, 0.01, label="b:", **common)
        gui.doubleSpin(
            coef_box, self, "y_offset", 0.0, 1.0, 0.01, label="c:", **common)

        gui.rubber(self.controlArea)
        gui.auto_commit(self.controlArea, self, "auto_commit",
                        "Send Selection", "Send Automatically")

        self.setup_info_label()
        self.enable_controls()

    def __filter_type_changed(self):
        self.enable_controls()
        self.__param_changed()

    def __param_changed(self):
        self.update_selection()
        self.setup_info_label()
        self.commit()

    def __manual_move(self, delta_x, delta_y):
        self.x_offset += delta_x
        self.y_offset += delta_y
        if self.x_offset < 0:
            self.x_offset = 0
        if self.y_offset < 0:
            self.y_offset = 0
        self.filter_type = FilterType.ByEquation
        self.__filter_type_changed()

    @property
    def filter_by_nr_of_genes(self):
        return self.filter_type == FilterType.ByNumber

    @Inputs.data
    def set_data(self, data):
        self.closeContext()
        self.data = data
        self.openContext(data)
        self.select_genes()
        self.setup_info_label()
        self.unconditional_commit()

    @Inputs.genes
    def set_genes(self, genes):
        self.genes = genes
        self.check_genes()
        self.graph.update_markers(self.data, self.genes)

    def check_genes(self):
        self.Warning.missing_entrez_id.clear()
        if self.genes:
            if ENTREZ_ID not in self.genes.domain:
                self.Warning.missing_entrez_id()
                self.genes = None

    def select_genes(self):
        self.selected = None
        self.graph.clear_all()
        self.Warning.less_selected.clear()
        if not self.data:
            return

        selector = self.__get_selector()
        self.zero_rate, self.mean_expr = selector.detection(self.data.X)
        results = self.__select(selector)
        self.graph.set_data(DropoutResults(*results), self.data, self.genes)

    def __get_selector(self):
        kwargs = {"decay": self.decay, "y_offset": self.y_offset}
        if self.filter_by_nr_of_genes:
            kwargs["n_genes"] = self.n_genes
        else:
            kwargs["x_offset"] = self.x_offset
        return DropoutGeneSelection(**kwargs)

    def __select(self, selector):
        self.selected = selector.select_genes(self.zero_rate, self.mean_expr)
        n_selected = sum(self.selected)
        if n_selected < self.n_genes and self.filter_by_nr_of_genes:
            self.Warning.less_selected(n_selected)
        self.n_genes = n_selected
        self.decay = selector.decay
        self.x_offset = selector.x_offset
        self.y_offset = selector.y_offset
        return (self.zero_rate, self.mean_expr, selector.decay,
                selector.x_offset, selector.y_offset, selector.threshold)

    def update_selection(self):
        self.Warning.less_selected.clear()
        if not self.data:
            return
        selector = self.__get_selector()
        self.graph.update_curve(DropoutResults(*self.__select(selector)))

    def setup_info_label(self):
        text = "No data on input."
        if self.selected is not None:
            k = sum(self.selected)
            n, m = len(self.data), len(self.data.domain.attributes)
            ks = "s" if k != 1 else ""
            ns, ms = "s" if n != 1 else "", "s" if m != 1 else ""
            text = f"Data with {n} cell{ns} and {m} gene{ms}" \
                   f"\n{k} gene{ks} in selection"
        self.info_label.setText(text)

    def enable_controls(self):
        self.controls.n_genes.setEnabled(self.filter_by_nr_of_genes)
        self.controls.decay.setEnabled(not self.filter_by_nr_of_genes)
        self.controls.x_offset.setEnabled(not self.filter_by_nr_of_genes)
        self.controls.y_offset.setEnabled(not self.filter_by_nr_of_genes)

    def commit(self):
        data = None
        if self.selected is not None:
            data = DropoutGeneSelection.filter_columns(self.data,
                                                       self.selected)
        self.Outputs.data.send(data)

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(800, 400))

    def send_report(self):
        if self.selected is None:
            return
        self.report_plot()
        self.report_caption(report.render_items_vert((
            ("Number of genes", self.n_genes),
            ("a", round(self.decay, 2)),
            ("b", round(self.x_offset, 2)),
            ("a", round(self.y_offset, 2)))))


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    table = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")
    WidgetPreview(OWDropout).run(set_data=table)

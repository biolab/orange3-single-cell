import sys
from types import SimpleNamespace
from typing import Optional, Sequence

import numpy as np
from scipy import stats

from AnyQt.QtCore import Qt, QSize, QPointF, QRectF, QLineF, QTimer
from AnyQt.QtCore import pyqtSignal as Signal
from AnyQt.QtGui import (
    QPolygonF, QPainterPath, QPalette, QPen, QBrush, QColor, QKeySequence
)
from AnyQt.QtWidgets import (
    QLabel, QSpinBox, QGroupBox, QAction, QGraphicsPathItem, QGraphicsRectItem,
    QFormLayout, QApplication, QButtonGroup, QRadioButton
)

import pyqtgraph as pg

import Orange.data
import Orange.widgets.utils.plot.owpalette
from Orange.widgets import widget, gui, settings

#: Filter type
Cells, Genes, Data = 0, 1, 2

#: Filter descriptions for various roles in UI
#: (short name, name, description)
FilterInfo = {
    Cells: ("Cells", "Cell Filter",
            "Filter cells (rows) by number of positive measurements"),
    Genes: ("Genes", "Gene Filter",
            "Filter genes (columns) by number of positive measurements"),
    Data: ("Data", "Data filter",
           "Filter out (zero) small measurements")
}


class OWFilter(widget.OWWidget):
    name = "Filter"
    description = "Filter cells/genes"

    class Inputs:
        data = widget.Input("Data", Orange.data.Table)

    class Outputs:
        data = widget.Output("Data", Orange.data.Table)

    class Warning(widget.OWWidget.Warning):
        invalid_range = widget.Msg(
            "Negative values in input data.\n"
            "This filter only makes sense for non-negative measurements"
            "where 0 indicates a lack (of) and/or a neutral reading."
        )
        sampling_in_effect = widget.Msg(
            "Too many data points to display.\n"
            "Sampling {} of {} data points."
        )

    Cells, Genes, Data = Cells, Genes, Data

    #: Augment the violin plot with a dot plot (strip plot) of the counts
    display_dotplot = settings.Setting(True)  # type: bool

    #: The selected filter type
    selected_filter_type = settings.Setting(Cells)  # type: int

    #: Is min/max count range selection enabled
    range_filter_enabled = settings.Setting(True)  # type: bool

    #: The lower and upper selection limits stored as absolute counts
    #: (does not transfer well between datasets)
    limit_lower = settings.Setting(0, schema_only=True)            # type: int
    limit_upper = settings.Setting(2 ** 31 - 1, schema_only=True)  # type: int

    auto_commit = settings.Setting(True)   # type: bool

    def __init__(self):
        super().__init__()
        self.data = None      # type: Optional[Orange.data.Table]
        self._counts = None   # type: Optional[np.ndarray]

        box = gui.widgetBox(self.controlArea, "Info")
        self._info = QLabel(box, wordWrap=True)
        self._info.setText("No data in input\n")

        box.layout().addWidget(self._info)

        box = gui.widgetBox(self.controlArea, "Filter Type")
        rbg = QButtonGroup(box, exclusive=True)
        for id_ in [Cells, Genes, Data]:
            name, _, tip = FilterInfo[id_]
            b = QRadioButton(
                name, toolTip=tip, checked=id_ == self.selected_filter_type
            )
            box.layout().addWidget(b)
            rbg.addButton(b, id_)
        rbg.buttonClicked[int].connect(self.set_filter_type)

        box = gui.widgetBox(self.controlArea, "View")
        self._showpoints = gui.checkBox(
            box, self, "display_dotplot", "Show data points",
            callback=self._update_dotplot
        )
        form = QFormLayout(
            labelAlignment=Qt.AlignLeft,
            formAlignment=Qt.AlignLeft,
            fieldGrowthPolicy=QFormLayout.AllNonFixedFieldsGrow
        )
        self._filter_box = box = gui.widgetBox(
            self.controlArea, "Filter", orientation=form
        )  # type: QGroupBox
        box.setCheckable(True)
        box.setChecked(self.range_filter_enabled)
        box.toggled.connect(self._set_filter_enabled)

        self.mincountspin = gui.spin(
            box, self, "limit_lower", 0, 2 ** 31 - 1, keyboardTracking=False,
            callback=self._limitchanged, addToLayout=False
        )  # type: QSpinBox
        # TODO: Selectively disable upper (lower?) bound (auto adjust to data)
        self.maxcountspin = gui.spin(
            box, self, "limit_upper", 0, 2 ** 31 - 1, keyboardTracking=False,
            callback=self._limitchanged, addToLayout=False
        )  # type: QSpinBox

        form.addRow("Min:", self.mincountspin)
        form.addRow("Max:", self.maxcountspin)

        self.controlArea.layout().addStretch(10)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")

        self._view = pg.GraphicsView()
        self._view.enableMouse(False)
        self._view.setAntialiasing(True)
        self._plot = plot = ViolinPlot()
        self._plot.setDataPointsVisible(self.display_dotplot)
        self._plot.setSelectionEnabled(self.range_filter_enabled)
        self._plot.selectionEdited.connect(self._limitchanged_plot)
        self._view.setCentralWidget(self._plot)
        self._plot.setTitle("Detected genes")

        left = self._plot.getAxis("left")  # type: pg.AxisItem
        left.setLabel("Detected genes")
        bottom = self._plot.getAxis("bottom")  # type: pg.AxisItem
        bottom.hide()
        plot.setMouseEnabled(False, False)
        plot.hideButtons()
        self.mainArea.layout().addWidget(self._view)

        # Coalescing commit timer
        self._committimer = QTimer(self, singleShot=True)
        self._committimer.timeout.connect(self.commit)

        self.addAction(
            QAction("Select All", self, shortcut=QKeySequence.SelectAll,
                    triggered=self._select_all)
        )

    def sizeHint(self):
        sh = super().sizeHint()  # type: QSize
        return sh.expandedTo(QSize(800, 600))

    def set_filter_type(self, type_):
        if self.selected_filter_type != type_:
            assert type_ in (Cells, Genes, Data), str(type_)
            self.selected_filter_type = type_
            if self.data is not None:
                self._setup(self.data, type_)

    def filter_type(self):
        return self.selected_filter_type

    @Inputs.data
    def set_data(self, data):
        # type: (Optional[Orange.data.Table]) -> None
        self.clear()
        self.data = data
        if data is not None:
            if np.any(data.X < 0):
                self.Warning.invalid_range()
            self._setup(data, self.filter_type())

        self.unconditional_commit()

    def clear(self):
        self._plot.clear()
        self.data = None
        self._counts = None
        self._update_info()
        self.Warning.clear()

    def _update_info(self):
        text = []
        if self.data is None:
            text += ["No data on input.\n"]
        else:
            N, M = len(self.data), len(self.data.domain.attributes)
            text = []
            text += [
                "Data with {N} cell{Np} and {M} gene{Mp}"
                .format(N=N, Np="s" if N != 1 else "",
                        M=M, Mp="s" if N != 1 else "")
            ]
            if self.range_filter_enabled and \
                    self.filter_type() in [Cells, Genes]:
                mask = ((self.limit_lower <= self._counts) &
                        (self._counts <= self.limit_upper))
                n = np.count_nonzero(mask)
                subject = "cell" if self.filter_type() == Cells else "gene"
                if n == 0:
                    text += ["All {}s filtered out".format(subject)]
                else:
                    text += [
                        "{} {subject}{s} in selection"
                        .format(n, subject=subject, s="s" if n != 1 else "")
                    ]
            else:
                text += [""]
        self._info.setText("\n".join(text))

    def _set_filter_enabled(self, value):
        self.range_filter_enabled = value
        self._filter_box.setChecked(value)
        self._plot.setSelectionEnabled(value)
        self._update_info()
        self._schedule_commit()

    def _select_all(self):
        self.limit_lower = 0
        self.limit_upper = 2 ** 31 - 1
        self._limitchanged()

    def _setup(self, data, filter_type):
        self._plot.clear()
        self._counts = None
        title = None
        sample_range = None

        if filter_type in [Cells, Genes]:
            if filter_type == Cells:
                axis = 1
                title = "Cell Filter"
                axis_label = "Detected Genes"
            else:
                axis = 0
                title = "Gene Filter"
                axis_label = "Detected Cells"

            mask = (data.X != 0) & (np.isfinite(data.X))
            counts = np.count_nonzero(mask, axis=axis)
            x = counts
            self._counts = counts
            self.Warning.sampling_in_effect.clear()
        elif filter_type == Data:
            x = data.X.ravel()
            x = x[np.isfinite(x)]
            self._counts = x
            MAX_DISPLAY_SIZE = 20000
            if x.size > MAX_DISPLAY_SIZE:
                self.Warning.sampling_in_effect(MAX_DISPLAY_SIZE, x.size)
                # tails to preserve exactly
                tails = 1
                x = np.sort(x)
                x1, x2, x3 = x[:tails], x[tails:x.size - tails], x[x.size-tails:]
                assert x1.size + x2.size + x3.size == x.size
                x2 = np.random.RandomState(0x667).choice(
                    x2, size=MAX_DISPLAY_SIZE - 2 * tails, replace=False,
                )
                x = np.r_[x1, x2, x3]
            else:
                self.Warning.sampling_in_effect.clear()
            title = "Data Filter"
            axis_label = "Gene Expression"
        else:
            assert False

        if x.size:
            xmin, xmax = np.min(x), np.max(x)
            self.limit_lower = np.clip(self.limit_lower, xmin, xmax)
            self.limit_upper = np.clip(self.limit_upper, xmin, xmax)

        if x.size > 0:
            # TODO: Need correction for lower bounded distribution (counts)
            # Use reflection around 0, but gaussian_kde does not provide
            # sufficient flexibility w.r.t bandwidth selection.
            self._plot.setData(x, 1000)
            self._plot.setBoundary(self.limit_lower, self.limit_upper)

        ax = self._plot.getAxis("left")  # type: pg.AxisItem
        ax.setLabel(axis_label)
        self._plot.setTitle(title)
        self._update_info()

    def _update_dotplot(self):
        self._plot.setDataPointsVisible(self.display_dotplot)

    def _limitchanged(self):
        # Low/high limit changed via the spin boxes
        if self._counts is not None and self._counts.size:
            xmin = np.min(self._counts)
            xmax = np.max(self._counts)
            self._plot.setBoundary(
                np.clip(self.limit_lower, xmin, xmax),
                np.clip(self.limit_upper, xmin, xmax)
            )
            # TODO: Only when the actual selection/filter mask changes
            self._schedule_commit()
            self._update_info()

    def _limitchanged_plot(self):
        # Low/high limit changed via the plot
        if self._counts is not None:
            self.limit_lower, self.limit_upper = self._plot.boundary()
            # TODO: Only when the actual selection/filter mask changes
            self._schedule_commit()
            self._update_info()

    def _schedule_commit(self):
        self._committimer.start()

    def commit(self):
        self._committimer.stop()
        data = self.data
        if data is not None and self.range_filter_enabled:
            if self.filter_type() in [Cells, Genes]:
                counts = self._counts
                cmax = self.limit_upper
                cmin = self.limit_lower
                mask = (cmin <= counts) & (counts <= cmax)
                if self.filter_type() == Cells:
                    assert counts.size == len(data)
                    data = data[mask]
                else:
                    assert counts.size == len(data.domain.attributes)
                    atts = [v for v, m in zip(data.domain.attributes, mask)
                            if m]
                    data = data.from_table(
                        Orange.data.Domain(
                            atts, data.domain.class_vars, data.domain.metas
                        ),
                        data
                    )
                if len(data) == 0 or \
                        len(data.domain) + len(data.domain.metas) == 0:
                    data = None
            elif self.filter_type() == Data:
                dmin, dmax = self.limit_lower, self.limit_upper
                data = data.copy()
                assert data.X.base is None
                mask = data.X < dmin
                mask |= data.X > dmax
                data.X[mask] = 0.0
            else:
                assert False

        self.Outputs.data.send(data)

    def onDeleteWidget(self):
        self.clear()
        self._plot.close()
        super().onDeleteWidget()


class ViolinPlot(pg.PlotItem):
    """
    A violin plot item with interactive data boundary selection.
    """
    #: Emitted when the selection boundary has changed
    selectionChanged = Signal()
    #: Emitted when the selection boundary has been edited by the user
    #: (by dragging the boundary lines)
    selectionEdited = Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.__data = None
        #: min/max cutoff line positions
        self.__min = 0
        self.__max = 0
        self.__dataPointsVisible = True
        self.__selectionEnabled = True
        self._plotitems = None

    def setData(self, data, nsamples, sample_range=None, color=Qt.magenta):
        assert np.all(np.isfinite(data))

        if data.size > 0:
            xmin, xmax = np.min(data), np.max(data)
        else:
            xmin = xmax = 0.0

        if sample_range is None:
            xrange = xmax - xmin
            sample_min = xmin - xrange * 0.025
            sample_max = xmax + xrange * 0.025
        else:
            sample_min, sample_max = sample_range

        sample = np.linspace(sample_min, sample_max, nsamples)
        if data.size < 2:
            est = np.full(sample.size, 1. / sample.size, )
        else:
            try:
                density = stats.gaussian_kde(data)
                est = density.evaluate(sample)
            except np.linalg.LinAlgError:
                est = np.zeros(sample.size)

        item = QGraphicsPathItem(violin_shape(sample, est))
        color = QColor(color)
        color.setAlphaF(0.5)
        item.setBrush(QBrush(color))
        pen = QPen(self.palette().color(QPalette.Shadow))
        pen.setCosmetic(True)
        item.setPen(pen)
        est_max = np.max(est)

        x = np.random.RandomState(0xD06F00D).uniform(
            -est_max, est_max, size=data.size
        )
        dots = pg.ScatterPlotItem(
            x=x, y=data, size=3,
        )
        dots.setVisible(self.__dataPointsVisible)
        cmax = SelectionLine(
            angle=0, pos=xmax, movable=True, bounds=(xmin, xmax)
        )
        cmin = SelectionLine(
            angle=0, pos=xmin, movable=True, bounds=(xmin, xmax)
        )
        cmax.setCursor(Qt.SizeVerCursor)
        cmin.setCursor(Qt.SizeVerCursor)

        selection_item = QGraphicsRectItem(
            QRectF(-est_max, xmin, est_max * 2, xmax - xmin)
        )
        selection_item.setPen(QPen(Qt.NoPen))
        selection_item.setBrush(QColor(0, 250, 0, 50))

        def update_selection_rect():
            rect = selection_item.rect()  # type: QRectF
            rect.setTop(cmax.value())
            rect.setBottom(cmin.value())
            selection_item.setRect(rect.normalized())

        cmax.sigPositionChanged.connect(update_selection_rect)
        cmin.sigPositionChanged.connect(update_selection_rect)

        def setupper(line):
            upper = line.value()
            lower = min(self.__min, upper)
            if (lower, upper) != (self.__min, self.__max):
                self.__min = lower
                self.__max = upper
                cmin.setValue(self.__min)
                self.selectionEdited.emit()
                self.selectionChanged.emit()

        def setlower(line):
            lower = line.value()
            upper = max(self.__max, lower)
            if (lower, upper) != (self.__min, self.__max):
                self.__min = lower
                self.__max = upper
                cmax.setValue(self.__max)
                self.selectionChanged.emit()
                self.selectionEdited.emit()

        cmax.sigPositionChanged.connect(setupper)
        cmin.sigPositionChanged.connect(setlower)

        cmax.setVisible(self.__selectionEnabled)
        cmin.setVisible(self.__selectionEnabled)
        selection_item.setVisible(self.__selectionEnabled)

        self.addItem(dots)
        self.addItem(item)
        self.addItem(cmax)
        self.addItem(cmin)
        self.addItem(selection_item)

        self.setRange(
            QRectF(-est_max, np.min(sample), est_max * 2, np.ptp(sample))
        )
        self._plotitems = SimpleNamespace(
            pointsitem=dots,
            densityitem=item,
            cmax=cmax,
            cmin=cmin,
            selection_item=selection_item
        )
        self.__min = xmin
        self.__max = xmax

    def setDataPointsVisible(self, visible):
        self.__dataPointsVisible = visible
        if self._plotitems is not None:
            self._plotitems.pointsitem.setVisible(visible)

    def setSelectionEnabled(self, enabled):
        self.__selectionEnabled = enabled
        if self._plotitems is not None:
            items = [self._plotitems.cmax,
                     self._plotitems.cmin,
                     self._plotitems.selection_item]

            for item in items:
                item.setVisible(enabled)

    def setBoundary(self, low, high):
        """
        Set the lower and upper selection boundary value.
        """
        changed = False
        if self.__min != low:
            self.__min = low
            changed = True
        if self.__max != high:
            self.__max = high
            changed = True

        if changed:
            if self._plotitems:
                b = self.blockSignals(True)
                try:
                    self._plotitems.cmin.setValue(low)
                    self._plotitems.cmax.setValue(high)
                finally:
                    self.blockSignals(b)

            self.selectionChanged.emit()

    def boundary(self):
        """
        Return the current lower and upper selection boundary values.
        """
        return self.__min, self.__max

    def clear(self):
        super().clear()
        self._plotitems = None


def violin_shape(x, p):
    # type: (Sequence[float], Sequence[float]) -> QPainterPath
    points = [QPointF(pi, xi) for xi, pi in zip(x, p)]
    points += [QPointF(-pi, xi) for xi, pi in reversed(list(zip(x, p)))]
    poly = QPolygonF(points)
    path = QPainterPath()
    path.addPolygon(poly)
    return path


class SelectionLine(pg.InfiniteLine):
    def paint(self, painter, option, widget=None):
        brect = self.boundingRect()
        c = brect.center()
        line = QLineF(brect.left(), c.y(), brect.right(), c.y())
        t = painter.transform()
        line = t.map(line)
        painter.save()
        painter.resetTransform()
        painter.setPen(self.currentPen)
        painter.drawLine(line)
        painter.restore()


def main(argv=None):
    app = QApplication(list(argv or sys.argv))
    argv = app.arguments()
    w = OWFilter()
    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = "brown-selected"  # bad example
    data = Orange.data.Table(filename)
    w.set_data(data)
    w.show()
    w.raise_()
    app.exec()
    w.onDeleteWidget()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

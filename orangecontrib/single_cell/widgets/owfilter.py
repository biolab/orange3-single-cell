import sys
import enum
import math
import numbers

from contextlib import contextmanager
from types import SimpleNamespace

import typing
from typing import Optional, Sequence, Tuple, Dict, Callable, Union, Iterable

import numpy as np
from scipy import stats

from AnyQt.QtCore import Qt, QSize, QPointF, QRectF, QLineF, QTimer
from AnyQt.QtCore import pyqtSignal as Signal, pyqtSlot as Slot
from AnyQt.QtGui import (
    QPainter, QPolygonF, QPainterPath, QPalette, QPen, QBrush, QColor,
    QKeySequence
)
from AnyQt.QtWidgets import (
    QLabel, QDoubleSpinBox, QGroupBox, QHBoxLayout, QAction,
    QGraphicsPathItem, QGraphicsRectItem, QGraphicsItem, QFormLayout,
    QApplication, QButtonGroup, QRadioButton, QCheckBox, QStackedWidget
)

import pyqtgraph as pg

import Orange.data
import Orange.widgets.utils.plot.owpalette
from Orange.widgets import widget, gui, settings

#: Filter type
Cells, Genes, Data = 0, 1, 2

#: Filter quality control measure (apply to Cell/Genes type only)
DetectionCount = 0  # number of genes/features with non-zero expression level
TotalCounts = 1  # total counts by cell/gene


#: Filter descriptions for various roles in UI
#: (short name, name, description)
FilterInfo = {
    Cells: ("Cells", "Cell Filter",
            "Filter cells (rows) by total counts (library size) or number "
            "of expressed genes."),
    Genes: ("Genes", "Gene Filter",
            "Filter genes (columns) by total counts mapped to a gene or "
            "number of cells in which the gene is expressed in."),
    Data: ("Data", "Data filter",
           "Filter out (zero) small measurements")
}

# Quality control measure descriptions for UI
MeasureInfo = {
    TotalCounts: ("Total counts",
                  "Sum of all counts across cell/gene"),
    DetectionCount: ("Detection count",
                     "Number of cells/genes with non-zero expression")
}


# Plot scales
class Scale(enum.Enum):
    Linear = "Linear"
    Log1p = "Log1p"


def log1p(x, base=10.):
    x = np.log1p(x)
    x /= np.log(base)
    return x


def expm1(x, base=10.):
    x = np.asarray(x)
    x *= np.log(base)
    return np.expm1(x)


class ScatterPlotItem(pg.ScatterPlotItem):
    def paint(self, painter, *args):
        if self.opts["antialias"]:
            painter.setRenderHint(QPainter.Antialiasing, True)
        if self.opts["pxMode"]:
            painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        super().paint(painter, *args)


if typing.TYPE_CHECKING:
    ArrayLike = Union[np.ndarray, np.generic, numbers.Number, Iterable]
    Transform = Callable[[ArrayLike], np.ndarray]


# filter state data class
class _FilterData:
    #: The array used for filtering/ploting. Is some marginal statistic
    #: of the input.
    x = ...   # type: np.ndarray
    #: The transformed array x (log1p).
    #: NOTE: xt.size need not match x.size (non finite values can be omitted)
    xt = ...  # type: np.ndarray
    #: The transformation function mapping x to xt
    transform = ...      # type: Transform
    #: The inverse of `transform`,
    transform_inv = ...  # type: Transform
    #: min/max bounds of x
    xmin = xmax = ...    # type: float
    #: min/max bounds of xt
    xtmin = xtmax = ...  # type: float


class OWFilter(widget.OWWidget):
    name = "Filter"
    icon = 'icons/Filter.svg'
    description = "Filter cells/genes"
    priority = 210

    class Inputs:
        data = widget.Input("Data", Orange.data.Table)

    class Outputs:
        data = widget.Output("Data", Orange.data.Table)

    class Warning(widget.OWWidget.Warning):
        sampling_in_effect = widget.Msg(
            "Too many data points to display.\n"
            "Sampling {} of {} data points."
        )

    class Error(widget.OWWidget.Error):
        invalid_range = widget.Msg(
            "Negative values in input data.\n"
            "This filter is only defined for non-negative values."
        )
        invalid_domain = widget.Msg(
            "Invalid domain\n"
            "Domain contains non numeric columns."
        )

    #: Filter mode.
    #: Filter out rows/columns or 'zap' data values in range.
    Cells, Genes, Data = Cells, Genes, Data

    settings_version = 3

    #: The selected filter mode
    selected_filter_type = settings.Setting(Cells)  # type: int

    #: Selected filter statistics / QC measure indexed by filter_type
    selected_filter_metric = settings.Setting(TotalCounts)  # type: int

    #: Augment the violin plot with a dot plot (strip plot) of the (non-zero)
    #: measurement counts in Cells/Genes mode or data matrix values in Data
    #: mode.
    display_dotplot = settings.Setting(True)  # type: bool

    #: Is min/max range selection enable
    limit_lower_enabled = settings.Setting(True)  # type: bool
    limit_upper_enabled = settings.Setting(True)  # type: bool

    #: The lower and upper selection limit for each filter type
    thresholds = settings.Setting({
        (Cells, DetectionCount): (0, 2 ** 31 - 1),
        (Cells, TotalCounts): (0, 2 ** 31 - 1),
        (Genes, DetectionCount): (0, 2 ** 31 - 1),
        (Genes, TotalCounts): (0, 2 ** 31 - 1),
        (Data, -1): (0.0, 2.0 ** 31 - 1)
    })  # type: Dict[Tuple[int, int], Tuple[float, float]]

    #: Plot scale: 'Linear' or 'Log1p'
    scale = settings.Setting(Scale.Linear.name)  # type: str

    auto_commit = settings.Setting(True)   # type: bool

    def __init__(self):
        super().__init__()
        self.data = None      # type: Optional[Orange.data.Table]
        self._state = None     # type: Optional[_FilterData]

        box = gui.widgetBox(self.controlArea, "Info")
        self._info = QLabel(box, wordWrap=True)
        self._info.setText("No data in input\n")

        box.layout().addWidget(self._info)

        box = gui.widgetBox(self.controlArea, "Filter Type", spacing=-1)
        rbg = QButtonGroup(box, exclusive=True)
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        for id_ in [Cells, Genes, Data]:
            name, _, tip = FilterInfo[id_]
            b = QRadioButton(
                name, toolTip=tip, checked=id_ == self.selected_filter_type
            )
            rbg.addButton(b, id_)
            layout.addWidget(b, stretch=10, alignment=Qt.AlignCenter)
        box.layout().addLayout(layout)

        rbg.buttonClicked[int].connect(self.set_filter_type)

        self.filter_metric_cb = gui.comboBox(
            box, self, "selected_filter_metric", callback=self._update_metric,
            enabled=self.selected_filter_type != Data
        )
        for id_ in [DetectionCount, TotalCounts]:
            text, ttip = MeasureInfo[id_]
            self.filter_metric_cb.addItem(text)
            idx = self.filter_metric_cb.count() - 1
            self.filter_metric_cb.setItemData(idx, ttip, Qt.ToolTipRole)
        self.filter_metric_cb.setCurrentIndex(self.selected_filter_metric)

        form = QFormLayout(
            labelAlignment=Qt.AlignLeft,
            formAlignment=Qt.AlignLeft,
            fieldGrowthPolicy=QFormLayout.AllNonFixedFieldsGrow
        )
        self._filter_box = box = gui.widgetBox(
            self.controlArea, "Filter", orientation=form
        )  # type: QGroupBox

        self.threshold_stacks = (
            QStackedWidget(enabled=self.limit_lower_enabled),
            QStackedWidget(enabled=self.limit_upper_enabled),
        )
        finfo = np.finfo(np.float64)
        for filter_ in [Cells, Genes, Data]:
            if filter_ in {Cells, Genes}:
                minimum = 0.0
                ndecimals = 1
                metric = self.selected_filter_metric
            else:
                minimum = finfo.min
                ndecimals = 3
                metric = -1
            spinlower = QDoubleSpinBox(
                self, minimum=minimum, maximum=finfo.max, decimals=ndecimals,
                keyboardTracking=False,
            )
            spinupper = QDoubleSpinBox(
                self, minimum=minimum, maximum=finfo.max, decimals=ndecimals,
                keyboardTracking=False,
            )

            lower, upper = self.thresholds.get((filter_, metric), (0, 0))

            spinlower.setValue(lower)
            spinupper.setValue(upper)

            self.threshold_stacks[0].addWidget(spinlower)
            self.threshold_stacks[1].addWidget(spinupper)

            spinlower.valueChanged.connect(self._limitchanged)
            spinupper.valueChanged.connect(self._limitchanged)

        self.threshold_stacks[0].setCurrentIndex(self.selected_filter_type)
        self.threshold_stacks[1].setCurrentIndex(self.selected_filter_type)

        self.limit_lower_enabled_cb = cb = QCheckBox(
            "Min", checked=self.limit_lower_enabled
        )
        cb.toggled.connect(self.set_lower_limit_enabled)
        cb.setAttribute(Qt.WA_LayoutUsesWidgetRect, True)
        form.addRow(cb, self.threshold_stacks[0])

        self.limit_upper_enabled_cb = cb = QCheckBox(
            "Max", checked=self.limit_upper_enabled
        )
        cb.toggled.connect(self.set_upper_limit_enabled)
        cb.setAttribute(Qt.WA_LayoutUsesWidgetRect, True)
        form.addRow(cb, self.threshold_stacks[1])

        box = gui.widgetBox(self.controlArea, "Plot Options")
        self._showpoints = gui.checkBox(
            box, self, "display_dotplot", "Show data points",
            callback=self._update_dotplot
        )
        self.log_scale_cb = QCheckBox(
            "Log scale", checked=self.scale == Scale.Log1p.name
        )
        self.log_scale_cb.toggled[bool].connect(
            lambda state:
                self.set_filter_scale(Scale.Log1p if state else Scale.Linear)
        )
        box.layout().addWidget(self.log_scale_cb)

        self.controlArea.layout().addStretch(10)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")

        self._view = pg.GraphicsView()
        self._view.enableMouse(False)
        self._view.setAntialiasing(True)
        self._plot = plot = ViolinPlot()
        self._plot.setDataPointsVisible(self.display_dotplot)
        self._plot.setSelectionMode(
            (ViolinPlot.Low if self.limit_lower_enabled else 0) |
            (ViolinPlot.High if self.limit_upper_enabled else 0)
        )
        self._plot.setRange(QRectF(-1., 0., 2., 1.))
        self._plot.selectionEdited.connect(self._limitchanged_plot)
        self._view.setCentralWidget(self._plot)

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
        self._setup_axes()

    def sizeHint(self):
        sh = super().sizeHint()  # type: QSize
        return sh.expandedTo(QSize(800, 600))

    def set_filter_type(self, type_):
        if self.selected_filter_type != type_:
            assert type_ in (Cells, Genes, Data), str(type_)
            self.selected_filter_type = type_
            self.threshold_stacks[0].setCurrentIndex(type_)
            self.threshold_stacks[1].setCurrentIndex(type_)
            self.filter_metric_cb.setEnabled(type_ != Data)
            self._setup_axes()
            if self.data is not None:
                self._setup(self.data, type_)
                self._schedule_commit()

    def filter_type(self):
        return self.selected_filter_type

    def set_filter_scale(self, scale):
        # type: (Scale) -> None
        if self.scale != scale:
            self.scale = scale.name
            self.log_scale_cb.setChecked(scale == Scale.Log1p)
            self._update_scale()

    def filter_scale(self):
        return Scale[self.scale]

    def _update_metric(self):
        self._update_scale()
        if self.data is not None:
            self._setup(self.data, self.selected_filter_type, )

    def set_upper_limit_enabled(self, enabled):
        if enabled != self.limit_upper_enabled:
            self.limit_upper_enabled = enabled
            self.threshold_stacks[1].setEnabled(enabled)
            self.limit_upper_enabled_cb.setChecked(enabled)
            self._update_filter()
            self._schedule_commit()

    def set_lower_limit_enabled(self, enabled):
        if enabled != self.limit_lower_enabled:
            self.limit_lower_enabled = enabled
            self.threshold_stacks[0].setEnabled(enabled)
            self.limit_lower_enabled_cb.setChecked(enabled)
            self._update_filter()
            self._schedule_commit()

    def _update_filter(self):
        mode = 0
        if self.limit_lower_enabled:
            mode |= ViolinPlot.Low
        if self.limit_upper_enabled:
            mode |= ViolinPlot.High
        self._plot.setSelectionMode(mode)
        self._update_info()
        self._schedule_commit()

    def _is_filter_enabled(self):
        return self.limit_lower_enabled or self.limit_upper_enabled

    @Inputs.data
    def set_data(self, data):
        # type: (Optional[Orange.data.Table]) -> None
        self.clear()

        if data is not None and \
                any(type(v) is not Orange.data.ContinuousVariable
                    for v in data.domain.attributes):
            self.Error.invalid_domain()
            data = None

        if data is not None and np.any(data.X < 0):
            self.Error.invalid_range()
            data = None

        self.data = data

        if data is not None:
            self._setup(data, self.filter_type())

        self.unconditional_commit()

    def clear(self):
        self.data = None
        self._state = None
        self._plot.clear()
        # reset the plot range
        self._plot.setRange(QRectF(-1., 0., 2., 1.))
        self._update_info()
        self.Warning.clear()
        self.Error.clear()

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
            if self._is_filter_enabled() and \
                    self.filter_type() in [Cells, Genes]:
                counts = self._state.x
                mask = np.ones(counts.shape, dtype=bool)
                if self.limit_lower_enabled:
                    mask &= self.limit_lower <= counts

                if self.limit_upper_enabled:
                    mask &= counts <= self.limit_upper

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

    def _select_all(self):
        self.limit_lower = 0
        self.limit_upper = 2 ** 31 - 1
        self._limitchanged()

    def _setup_axes(self):
        # Setup the plot axes and title
        filter_type = self.filter_type()
        info = FilterInfo[filter_type]
        _, title, _, *_ = info

        if filter_type in [Cells, Genes]:
            measure = self.selected_filter_metric
        else:
            measure = None

        if filter_type == Cells and measure == TotalCounts:
            axis_label = "Total counts (library size)"
        elif filter_type == Cells and measure == DetectionCount:
            axis_label = "Number of expressed genes"
        elif filter_type == Genes and measure == TotalCounts:
            axis_label = "Total counts"
        elif filter_type == Genes and measure == DetectionCount:
            # TODO: Too long
            axis_label = "Number of cells a gene is expressed in"
        elif filter_type == Data:
            axis_label = "Gene Expression"

        ax = self._plot.getAxis("left")
        if self.filter_scale() == Scale.Log1p:
            axis_label = "1 + '{}' <i>(in log scale)</i>".format(axis_label)
            ax.setLabel(axis_label)
            ax.setLogMode(True)
        else:
            ax.setLogMode(False)
            ax.setLabel(axis_label)
        # Reset the tick text area width
        ax.textWidth = 30
        ax.setWidth(None)

        self._plot.setTitle(title)

    def _setup(self, data, filter_type):
        self._plot.clear()
        self._state = None
        self._setup_axes()

        span = -1.0  # data span
        measure = self.selected_filter_metric if filter_type != Data else None
        state = _FilterData()

        if filter_type in [Cells, Genes]:
            if filter_type == Cells:
                axis = 1
            else:
                axis = 0
            if measure == TotalCounts:
                counts = np.nansum(data.X, axis=axis)
            else:
                mask = (data.X != 0) & (np.isfinite(data.X))
                counts = np.count_nonzero(mask, axis=axis)
            x = counts
            self.Warning.sampling_in_effect.clear()
        elif filter_type == Data:
            x = data.X.ravel()
            x = x[np.isfinite(x)]
            x = x[x != 0]
            MAX_DISPLAY_SIZE = 20000
            if x.size > MAX_DISPLAY_SIZE:
                self.Warning.sampling_in_effect(MAX_DISPLAY_SIZE, x.size)
                # tails to preserve exactly
                tails = 1
                assert x.flags.owndata
                x.partition(tails - 1)
                xrest = x[tails:]
                xrest.partition(xrest.size - tails)

                x1, x2, x3 = x[:tails], x[tails:x.size - tails], x[x.size-tails:]
                assert x1.size + x2.size + x3.size == x.size
                x2 = np.random.RandomState(0x667).choice(
                    x2, size=MAX_DISPLAY_SIZE - 2 * tails, replace=False,
                )
                x = np.r_[x1, x2, x3]
            else:
                self.Warning.sampling_in_effect.clear()
        else:
            assert False

        state.x = x

        scale = self.filter_scale()

        if scale == Scale.Log1p:
            scale_transform = log1p
            scale_transform_inv = expm1
        else:
            scale_transform = lambda x: x
            scale_transform_inv = scale_transform

        state.transform = scale_transform
        state.transform_inv = scale_transform_inv

        if x.size:
            xmin, xmax = np.min(x), np.max(x)
        else:
            xmin = xmax = 0., 1.

        state.xmin, state.xmax = xmin, xmax

        xs = scale_transform(x)
        xs = xs[np.isfinite(xs)]
        state.xt = xs

        if xs.size:
            xsmin, xsmax = np.min(xs), np.max(xs)
            # find effective xmin, xmax (valid in both original and transformed
            # space
            xmin_, xmax_ = scale_transform_inv([xsmin, xsmax])
            xmin, xmax = max(xmin, xmin_), min(xmax, xmax_)

            lower = np.clip(self.limit_lower, xmin, xmax)
            upper = np.clip(self.limit_upper, xmin, xmax)
        else:
            xmin, xmax = 0., 1.
            lower, upper = 0., 1.

        state.xtmin, state.xtmax = xsmin, xsmax

        spinlow = self.threshold_stacks[0].widget(filter_type)
        spinhigh = self.threshold_stacks[1].widget(filter_type)
        if filter_type == Data or measure == TotalCounts:
            span = xmax - xmin
            if span > 0:
                ndecimals = max(4 - int(np.floor(np.log10(span))), 1)
            else:
                ndecimals = 1
        else:
            ndecimals = 1

        # Round effective bounds (spin <=> plot cut lines)
        lower = round(lower, ndecimals)
        upper = round(upper, ndecimals)

        if xs.size > 0:
            # TODO: Need correction for lower bounded distribution (counts)
            # Use reflection around 0, but gaussian_kde does not provide
            # sufficient flexibility w.r.t bandwidth selection.

            self._plot.setData(xs, 1000)
            self._plot.setBoundary(*scale_transform([lower, upper]))

        spinlow.setDecimals(ndecimals)
        self.limit_lower = lower

        spinhigh.setDecimals(ndecimals)
        self.limit_upper = upper

        self._state = state
        self._update_info()

    def _update_dotplot(self):
        self._plot.setDataPointsVisible(self.display_dotplot)

    def current_filter_thresholds(self):
        if self.selected_filter_type in {Cells, Genes}:
            metric = self.selected_filter_metric
        else:
            metric = -1
        return self.thresholds[self.selected_filter_type, metric]

    def set_current_filter_thesholds(self, lower, upper):
        if self.selected_filter_type in {Cells, Genes}:
            metric = self.selected_filter_metric
        else:
            metric = -1
        self.thresholds[self.selected_filter_type, metric] = (lower, upper)

    def _update_scale(self):
        self._setup_axes()
        if self.data is not None:
            self._setup(self.data, self.filter_type())

    @property
    def limit_lower(self):
        return self.current_filter_thresholds()[0]

    @limit_lower.setter
    def limit_lower(self, value):
        _, upper = self.current_filter_thresholds()
        self.set_current_filter_thesholds(value, upper)
        stacklower, _ = self.threshold_stacks
        sb = stacklower.widget(self.selected_filter_type)
        # prevent changes due to spin box rounding
        sb.setValue(value)

    @property
    def limit_upper(self):
        return self.current_filter_thresholds()[1]

    @limit_upper.setter
    def limit_upper(self, value):
        lower, _ = self.current_filter_thresholds()
        self.set_current_filter_thesholds(lower, value)
        _, stackupper = self.threshold_stacks
        sb = stackupper.widget(self.selected_filter_type)
        sb.setValue(value)

    @Slot()
    def _limitchanged(self):
        # Low/high limit changed via the spin boxes
        stacklow, stackhigh = self.threshold_stacks
        filter_ = self.selected_filter_type

        lower = stacklow.widget(filter_).value()
        upper = stackhigh.widget(filter_).value()
        self.set_current_filter_thesholds(lower, upper)

        state = self._state
        if state is not None and state.x.size:
            xmin, xmax = state.xmin, state.xmax
            lower = np.clip(lower, xmin, xmax)
            upper = np.clip(upper, xmin, xmax)
            lower, upper = state.transform([lower, upper])
            self._plot.setBoundary(lower, upper)
            # TODO: Only when the actual selection/filter mask changes
            self._schedule_commit()
            self._update_info()

    def _limitchanged_plot(self):
        # Low/high limit changed via the plot
        if self._state is not None:
            state = self._state
            newlower_, newupper_ = self._plot.boundary()
            newlower, newupper = state.transform_inv([newlower_, newupper_])
            filter_ = self.selected_filter_type
            lower, upper = self.current_filter_thresholds()
            stacklow, stackhigh = self.threshold_stacks
            spin_lower = stacklow.widget(filter_)
            spin_upper = stackhigh.widget(filter_)
            # do rounding to match the spin box's precision
            if self.limit_lower_enabled:
                newlower = round(newlower, spin_lower.decimals())
            else:
                newlower = lower

            if self.limit_upper_enabled:
                newupper = round(newupper, spin_upper.decimals())
            else:
                newupper = upper

            if self.limit_lower_enabled and newlower != lower:
                self.limit_lower = newlower
            if self.limit_upper_enabled and newupper != upper:
                self.limit_upper = newupper

            newlower_, newupper_ = state.transform([newlower, newupper])
            self._plot.setBoundary(newlower_, newupper_)

            # TODO: Only when the actual selection/filter mask changes
            self._schedule_commit()
            self._update_info()

    def _schedule_commit(self):
        self._committimer.start()

    def commit(self):
        self._committimer.stop()
        data = self.data

        if data is not None and self._is_filter_enabled():
            if self.filter_type() in [Cells, Genes]:
                state = self._state
                assert state is not None
                counts = state.x
                cmax = self.limit_upper
                cmin = self.limit_lower
                mask = np.ones(counts.shape, dtype=bool)
                if self.limit_lower_enabled:
                    mask &= cmin <= counts
                if self.limit_upper_enabled:
                    mask &= counts <= cmax

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
                mask = None
                if self.limit_lower_enabled:
                    mask = data.X < dmin
                if self.limit_upper_enabled:
                    if mask is not None:
                        mask |= data.X > dmax
                    else:
                        mask = data.X < dmax
                data.X[mask] = 0.0
            else:
                assert False

        self.Outputs.data.send(data)

    def onDeleteWidget(self):
        self.clear()
        self._plot.close()
        self._view.scene().clear()
        super().onDeleteWidget()

    @classmethod
    def migrate_settings(cls, settings, version):
        if (version is None or version < 2) and \
                ("limit_lower" in settings and "limit_upper" in settings):
            # v2 changed limit_lower, limit_upper to per filter limits stored
            # in a single dict
            lower = settings.pop("limit_lower")
            upper = settings.pop("limit_upper")
            settings["thresholds"] = {
                (Cells, TotalCounts): (lower, upper),
                (Cells, DetectionCount): (lower, upper),
                (Genes, TotalCounts): (lower, upper),
                (Genes, DetectionCount): (lower, upper),
                (Data, -1): (lower, upper),
            }
        if version == 2:
            thresholds = settings["thresholds"]
            c = thresholds.pop(Cells)
            g = thresholds.pop(Genes)
            d = thresholds.pop(Data)
            thresholds = {
                (Cells, TotalCounts): c,
                (Cells, DetectionCount): c,
                (Genes, TotalCounts): g,
                (Genes, DetectionCount): g,
                (Data, -1): d,
            }
            settings["thresholds"] = thresholds



@contextmanager
def block_signals(qobj):
    b = qobj.blockSignals(True)
    try:
        yield
    finally:
        qobj.blockSignals(b)


class AxisItem(pg.AxisItem):
    def logTickStrings(self, values, scale, spacing):
        # reimplemented
        values = [10 ** v for v in values]
        return [render_exp(v, 1) for v in values]


class ViolinPlot(pg.PlotItem):
    """
    A violin plot item with interactive data boundary selection.
    """
    #: Emitted when the selection boundary has changed
    selectionChanged = Signal()
    #: Emitted when the selection boundary has been edited by the user
    #: (by dragging the boundary lines)
    selectionEdited = Signal()

    #: Selection Flags
    NoSelection, Low, High = 0, 1, 2

    def __init__(self, *args, enableMenu=False, axisItems=None, **kwargs):
        if axisItems is None:
            axisItems = {}
        for position in ("left", 'right', 'top', 'bottom'):
            axisItems.setdefault(position, AxisItem(position))

        super().__init__(*args, enableMenu=enableMenu, axisItems=axisItems,
                         **kwargs)
        self.__data = None
        #: min/max cutoff line positions
        self.__min = 0
        self.__max = 0
        self.__dataPointsVisible = True
        self.__selectionEnabled = True
        self.__selectionMode = ViolinPlot.High | ViolinPlot.Low
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
        dots = ScatterPlotItem(
            x=x, y=data, size=3,
        )
        dots.setVisible(self.__dataPointsVisible)
        pen = QPen(self.palette().color(QPalette.Shadow), 1)
        hoverPen = QPen(self.palette().color(QPalette.Highlight), 1.5)
        cmax = SelectionLine(
            angle=0, pos=xmax, movable=True, bounds=(sample_min, sample_max),
            pen=pen, hoverPen=hoverPen
        )
        cmin = SelectionLine(
            angle=0, pos=xmin, movable=True, bounds=(sample_min, sample_max),
            pen=pen, hoverPen=hoverPen
        )
        cmax.setCursor(Qt.SizeVerCursor)
        cmin.setCursor(Qt.SizeVerCursor)

        selection_item = QGraphicsRectItem(
            QRectF(-est_max, xmin, est_max * 2, xmax - xmin)
        )
        selection_item.setPen(QPen(Qt.NoPen))
        selection_item.setBrush(QColor(0, 250, 0, 50))

        def update_selection_rect():
            mode = self.__selectionMode
            p = selection_item.parentItem()  # type: Optional[QGraphicsItem]
            while p is not None and not isinstance(p, pg.ViewBox):
                p = p.parentItem()
            if p is not None:
                viewbox = p  # type: pg.ViewBox
            else:
                viewbox = None
            rect = selection_item.rect()  # type: QRectF
            if mode & ViolinPlot.High:
                rect.setTop(cmax.value())
            elif viewbox is not None:
                rect.setTop(viewbox.viewRect().bottom())
            else:
                rect.setTop(cmax.maxRange[1])

            if mode & ViolinPlot.Low:
                rect.setBottom(cmin.value())
            elif viewbox is not None:
                rect.setBottom(viewbox.viewRect().top())
            else:
                rect.setBottom(cmin.maxRange[0])

            selection_item.setRect(rect.normalized())

        cmax.sigPositionChanged.connect(update_selection_rect)
        cmin.sigPositionChanged.connect(update_selection_rect)
        cmax.visibleChanged.connect(update_selection_rect)
        cmin.visibleChanged.connect(update_selection_rect)

        def setupper(line):
            ebound = self.__effectiveBoundary()
            elower, eupper = ebound
            mode = self.__selectionMode
            if not mode & ViolinPlot.High:
                return
            upper = line.value()
            lower = min(elower, upper)
            if lower != elower and mode & ViolinPlot.Low:
                self.__min = lower
                cmin.setValue(lower)

            if upper != eupper:
                self.__max = upper

            if ebound != self.__effectiveBoundary():
                self.selectionEdited.emit()
                self.selectionChanged.emit()

        def setlower(line):
            ebound = self.__effectiveBoundary()
            elower, eupper = ebound
            mode = self.__selectionMode
            if not mode & ViolinPlot.Low:
                return
            lower = line.value()
            upper = max(eupper, lower)
            if upper != eupper and mode & ViolinPlot.High:
                self.__max = upper
                cmax.setValue(upper)

            if lower != elower:
                self.__min = lower

            if ebound != self.__effectiveBoundary():
                self.selectionEdited.emit()
                self.selectionChanged.emit()

        cmax.sigPositionChanged.connect(setupper)
        cmin.sigPositionChanged.connect(setlower)
        selmode = self.__selectionMode
        cmax.setVisible(selmode & ViolinPlot.High)
        cmin.setVisible(selmode & ViolinPlot.Low)
        selection_item.setVisible(selmode)

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

    def setSelectionMode(self, mode):
        oldlower, oldupper = self.__effectiveBoundary()
        oldmode = self.__selectionMode
        mode = mode & 0b11
        if self.__selectionMode == mode:
            return

        self.__selectionMode = mode

        if self._plotitems is None:
            return

        cmin = self._plotitems.cmin
        cmax = self._plotitems.cmax
        selitem = self._plotitems.selection_item

        cmin.setVisible(mode & ViolinPlot.Low)
        cmax.setVisible(mode & ViolinPlot.High)
        selitem.setVisible(bool(mode))

        lower, upper = self.__effectiveBoundary()
        # The recorded values are not bounded by each other on gui interactions
        # when one is disabled. Rectify this now.
        if (oldmode ^ mode) & ViolinPlot.Low and mode & ViolinPlot.High:
            # Lower activated and High enabled
            lower = min(lower, upper)
        if (oldmode ^ mode) & ViolinPlot.High and mode & ViolinPlot.Low:
            # High activated and Low enabled
            upper = max(lower, upper)

        with block_signals(self):
            if lower != oldlower and mode & ViolinPlot.Low:
                cmin.setValue(lower)
            if upper != oldupper and mode & ViolinPlot.High:
                cmax.setValue(upper)

        self.selectionChanged.emit()

    def setBoundary(self, low, high):
        """
        Set the lower and upper selection boundary value.
        """
        changed = 0
        mode = self.__selectionMode
        if self.__min != low:
            self.__min = low
            changed |= mode & ViolinPlot.Low
        if self.__max != high:
            self.__max = high
            changed |= mode & ViolinPlot.High

        if changed:
            if self._plotitems:
                with block_signals(self):
                    if changed & ViolinPlot.Low:
                        self._plotitems.cmin.setValue(low)
                    if changed & ViolinPlot.High:
                        self._plotitems.cmax.setValue(high)

            self.selectionChanged.emit()

    def boundary(self):
        """
        Return the current lower and upper selection boundary values.
        """
        return self.__min, self.__max

    def __effectiveBoundary(self):
        # effective boundary, masked by selection mode
        low, high = -np.inf, np.inf
        if self.__selectionMode & ViolinPlot.Low:
            low = self.__min
        if self.__selectionMode & ViolinPlot.High:
            high = self.__max
        return low, high

    def clear(self):
        super().clear()
        self._plotitems = None

    def mouseDragEvent(self, event):
        mode = self.__selectionMode
        if mode != ViolinPlot.NoSelection and event.buttons() & Qt.LeftButton:
            start = event.buttonDownScenePos(Qt.LeftButton)  # type: QPointF
            pos = event.scenePos()  # type: QPointF
            cmin, cmax = self._plotitems.cmin, self._plotitems.cmax
            assert cmin.parentItem() is cmax.parentItem()
            pos = self.mapToItem(cmin.parentItem(), pos)
            start = self.mapToItem(cmin.parentItem(), start)
            if mode & ViolinPlot.Low and mode & ViolinPlot.High:
                lower, upper = min(pos.y(), start.y()), max(pos.y(), start.y())
                cmin.setValue(lower)
                cmax.setValue(upper)
            elif mode & ViolinPlot.Low:
                lower = pos.y()
                cmin.setValue(lower)
            elif mode & ViolinPlot.High:
                upper = pos.y()
                cmax.setValue(upper)
            event.accept()


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


def render_exp(value, prec=2):
    # type: (float, int) -> str
    if not math.isfinite(value):
        return repr(value)
    exp = "{:.{prec}G}".format(value, prec=prec)
    try:
        frac, exp = exp.split("E", 1)
    except ValueError:
        return exp

    frac = float(frac)
    exp = int(exp)
    if exp == 0:
        return str(frac)
    elif frac == 1.0:
        return "10{exp}".format(exp=_superscript(str(exp)))
    else:
        return "{frac:g}\u00D710{exp}".format(
            frac=frac, exp=_superscript(str(exp))
        )


def _superscript(string):
    # type: (str) -> str
    table = str.maketrans(
        "0123456789+-",
        "\u2070\u00B9\u00B2\u00B3\u2074\u2075\u2076\u2077\u2078\u2079"
        "\u207A\u207B",
    )
    return string.translate(table)


def main(argv=None):  # pragma: no cover
    app = QApplication(list(argv or sys.argv))
    argv = app.arguments()
    w = OWFilter()
    if len(argv) > 1:
        filename = argv[1]
        data = Orange.data.Table(filename)
    else:
        X = np.random.exponential(size=(1000, 1050)) - 1
        X[X < 0] = 0
        data = Orange.data.Table.from_numpy(None, X)
    w.set_data(data)
    w.show()
    w.raise_()
    app.exec()
    w.saveSettings()
    w.onDeleteWidget()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv))

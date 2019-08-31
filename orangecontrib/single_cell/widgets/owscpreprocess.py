import sys
import pkg_resources
import warnings

import numpy as np

from AnyQt.QtCore import Signal
from AnyQt.QtWidgets import (
    QApplication, QVBoxLayout, QHBoxLayout, QFormLayout, QSpinBox, QComboBox,
    QButtonGroup, QLabel, QDoubleSpinBox, QGroupBox, QCheckBox, QRadioButton
)

from Orange.data import DiscreteVariable
from Orange.preprocess.preprocess import PreprocessorList
import Orange.widgets.data.owpreprocess
from Orange.widgets.data.owpreprocess import (
    PreprocessAction, Description, index_to_enum, enum_to_index
)
from Orange.widgets.data.utils.preprocess import (
    ParametersRole, DescriptionRole, Controller, BaseEditor
)
from Orange.widgets.settings import (
    DomainContextHandler, ContextSetting, Setting
)
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import Input, Output, Msg

from orangecontrib.single_cell.preprocess.scpreprocess import (
    LogarithmicScale, Binarize, Normalize, NormalizeSamples, Standardize,
    SelectMostVariableGenes, NormalizeGroups, DropoutGeneSelection,
    DropoutWarning
)


def icon_path(basename):
    return pkg_resources.resource_filename(__name__, "icons/" + basename)


class ScBaseEditor(BaseEditor):
    def __init__(self, parent=None, master=None, **kwargs):
        super().__init__(parent, **kwargs)


class LogarithmicScaleEditor(ScBaseEditor):
    DEFAULT_BASE = LogarithmicScale.BinaryLog

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.setLayout(QVBoxLayout())

        form = QFormLayout()
        self.base_cb = QComboBox()
        self.base_cb.addItems(["2 (Binary Logarithm)",
                               "e (Natural Logarithm)",
                               "10 (Common Logarithm)"])
        self.base_cb.currentIndexChanged.connect(self.changed)
        self.base_cb.activated.connect(self.edited)

        form.addRow("Logarithm Base:", self.base_cb)
        self.layout().addLayout(form)

    def setParameters(self, params):
        base = params.get("base", self.DEFAULT_BASE)
        self.base_cb.setCurrentIndex(
            enum_to_index(LogarithmicScale.Base, base))

    def parameters(self):
        return {"base": index_to_enum(LogarithmicScale.Base,
                                      self.base_cb.currentIndex())}

    @staticmethod
    def createinstance(params):
        base = params.get("base", LogarithmicScaleEditor.DEFAULT_BASE)
        return LogarithmicScale(base)

    def __repr__(self):
        return "Base: {}".format(self.base_cb.currentText())


class BinarizeEditor(ScBaseEditor):
    DEFAULT_CONDITION = Binarize.GreaterOrEqual
    DEFAULT_THRESHOLD = 1

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self._threshold = self.DEFAULT_THRESHOLD

        self.setLayout(QVBoxLayout())
        form = QFormLayout()
        self.cond_cb = QComboBox()
        self.cond_cb.addItems(["Greater or Equal", "Greater"])
        self.cond_cb.currentIndexChanged.connect(self.changed)
        self.cond_cb.activated.connect(self.edited)

        self.thr_spin = QDoubleSpinBox(
            minimum=0, singleStep=0.5, decimals=1, value=self._threshold
        )
        self.thr_spin.valueChanged[float].connect(self._set_threshold)
        self.thr_spin.editingFinished.connect(self.edited)

        form.addRow("Condition:", self.cond_cb)
        form.addRow("Threshold:", self.thr_spin)
        self.layout().addLayout(form)

    def _set_threshold(self, t):
        if self._threshold != t:
            self._threshold = t
            self.thr_spin.setValue(t)
            self.changed.emit()

    def setParameters(self, params):
        cond = params.get("condition", self.DEFAULT_CONDITION)
        self.cond_cb.setCurrentIndex(enum_to_index(Binarize.Condition, cond))
        self._set_threshold(params.get("threshold", self.DEFAULT_THRESHOLD))

    def parameters(self):
        cond = index_to_enum(Binarize.Condition, self.cond_cb.currentIndex())
        return {"condition": cond, "threshold": self._threshold}

    @staticmethod
    def createinstance(params):
        condition = params.get("condition", BinarizeEditor.DEFAULT_CONDITION)
        threshold = params.get("threshold", BinarizeEditor.DEFAULT_THRESHOLD)
        return Binarize(condition, threshold)

    def __repr__(self):
        return "Condition: {}, Threshold: {}".format(
            self.cond_cb.currentText(), self.thr_spin.value()
        )


class NormalizeEditor(ScBaseEditor):
    DEFAULT_GROUP_BY = False
    DEFAULT_GROUP_VAR = None
    DEFAULT_METHOD = Normalize.CPM

    def __init__(self, parent=None, master=None, **kwargs):
        super().__init__(parent, **kwargs)
        self._group_var = self.DEFAULT_GROUP_VAR
        self._master = master
        self._master.input_data_changed.connect(self._set_model)
        self.setLayout(QVBoxLayout())

        form = QFormLayout()
        cpm_b = QRadioButton("Counts per million", checked=True)
        med_b = QRadioButton("Median")
        self.group = QButtonGroup()
        self.group.buttonClicked.connect(self._on_button_clicked)
        for i, button in enumerate([cpm_b, med_b]):
            index = index_to_enum(Normalize.Method, i).value
            self.group.addButton(button, index - 1)
            form.addRow(button)

        self.group_by_check = QCheckBox("Cell Groups: ",
                                        enabled=self.DEFAULT_GROUP_BY)
        self.group_by_check.clicked.connect(self.edited)
        self.group_by_combo = QComboBox(enabled=self.DEFAULT_GROUP_BY)
        self.group_by_model = DomainModel(
            order=(DomainModel.METAS, DomainModel.CLASSES),
            valid_types=DiscreteVariable,
            alphabetical=True
        )
        self.group_by_combo.setModel(self.group_by_model)
        self.group_by_combo.currentIndexChanged.connect(self.changed)
        self.group_by_combo.activated.connect(self.edited)

        form.addRow(self.group_by_check, self.group_by_combo)
        self.layout().addLayout(form)

        self._set_model()

    def _set_model(self):
        data = self._master.data
        self.group_by_model.set_domain(data and data.domain)
        enable = bool(self.group_by_model)
        self.group_by_check.setChecked(False)
        self.group_by_check.setEnabled(enable)
        self.group_by_combo.setEnabled(enable)
        if self.group_by_model:
            self.group_by_combo.setCurrentIndex(0)
            if self._group_var and self._group_var in data.domain:
                index = self.group_by_model.indexOf(self._group_var)
                self.group_by_combo.setCurrentIndex(index)
        else:
            self.group_by_combo.setCurrentText(None)

    def _on_button_clicked(self):
        self.changed.emit()
        self.edited.emit()

    def setParameters(self, params):
        method = params.get("method", self.DEFAULT_METHOD)
        index = enum_to_index(Normalize.Method, method)
        self.group.buttons()[index].setChecked(True)
        self._group_var = params.get("group_var", self.DEFAULT_GROUP_VAR)
        group = bool(self._group_var and self.group_by_model)
        if group:
            index = self.group_by_model.indexOf(self._group_var)
            self.group_by_combo.setCurrentIndex(index)
        group_by = params.get("group_by", self.DEFAULT_GROUP_BY)
        self.group_by_check.setChecked(group_by and group)

    def parameters(self):
        index = self.group_by_combo.currentIndex()
        group_var = self.group_by_model[index] if index > -1 else None
        group_by = self.group_by_check.isChecked()
        method = index_to_enum(Normalize.Method, self.group.checkedId())
        return {"group_var": group_var, "group_by": group_by, "method": method}

    @staticmethod
    def createinstance(params):
        group_var = params.get("group_var")
        group_by = params.get("group_by", NormalizeEditor.DEFAULT_GROUP_BY)
        method = params.get("method", NormalizeEditor.DEFAULT_METHOD)
        return NormalizeGroups(group_var, method) \
            if group_by and group_var else NormalizeSamples(method)

    def __repr__(self):
        method = self.group.button(self.group.checkedId()).text()
        index = self.group_by_combo.currentIndex()
        group_var = self.group_by_model[index] if index > -1 else None
        group_by = self.group_by_check.isChecked()
        group_text = ", Grouped by: {}".format(group_var) if group_by else ""
        return "Method: {}".format(method) + group_text


class StandardizeEditor(ScBaseEditor):
    DEFAULT_LOWER_CLIP = False
    DEFAULT_UPPER_CLIP = False
    DEFAULT_LOWER_BOUND = -10
    DEFAULT_UPPER_BOUND = 10

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self._lower_bound = self.DEFAULT_LOWER_BOUND
        self._upper_bound = self.DEFAULT_UPPER_BOUND

        self.setLayout(QVBoxLayout())

        box = QGroupBox(title="Clipping", flat=True)
        form = QFormLayout()
        self.lower_check = QCheckBox("Lower Bound: ")
        self.lower_check.clicked.connect(self.edited)
        self.lower_spin = QSpinBox(
            minimum=-99, maximum=0, value=self._lower_bound
        )
        self.lower_spin.valueChanged[int].connect(self._set_lower_bound)
        self.lower_spin.editingFinished.connect(self.edited)

        self.upper_check = QCheckBox("Upper Bound: ")
        self.upper_check.clicked.connect(self.edited)
        self.upper_spin = QSpinBox(value=self._upper_bound)
        self.upper_spin.valueChanged[int].connect(self._set_upper_bound)
        self.upper_spin.editingFinished.connect(self.edited)

        form.addRow(self.lower_check, self.lower_spin)
        form.addRow(self.upper_check, self.upper_spin)
        box.setLayout(form)
        self.layout().addWidget(box)

    def _set_lower_bound(self, x):
        if self._lower_bound != x:
            self._lower_bound = x
            self.lower_spin.setValue(x)
            self.changed.emit()

    def _set_upper_bound(self, x):
        if self._upper_bound != x:
            self._upper_bound = x
            self.upper_spin.setValue(x)
            self.changed.emit()

    def setParameters(self, params):
        lower_clip = params.get("lower_clip", self.DEFAULT_LOWER_CLIP)
        self.lower_check.setChecked(lower_clip)
        self._set_lower_bound(params.get("lower", self.DEFAULT_LOWER_BOUND))
        upper_clip = params.get("upper_clip", self.DEFAULT_UPPER_CLIP)
        self.upper_check.setChecked(upper_clip)
        self._set_upper_bound(params.get("upper", self.DEFAULT_UPPER_BOUND))

    def parameters(self):
        return {"lower_clip": self.lower_check.isChecked(),
                "lower": self._lower_bound,
                "upper_clip": self.upper_check.isChecked(),
                "upper": self._upper_bound}

    @staticmethod
    def createinstance(params):
        lower, upper = None, None
        if params.get("lower_clip", StandardizeEditor.DEFAULT_LOWER_CLIP):
            lower = params.get("lower", StandardizeEditor.DEFAULT_LOWER_BOUND)
        if params.get("upper_clip", StandardizeEditor.DEFAULT_UPPER_CLIP):
            upper = params.get("upper", StandardizeEditor.DEFAULT_UPPER_BOUND)
        return Standardize(lower, upper)

    def __repr__(self):
        clips = []
        if self.lower_check.isChecked():
            clips.append("Lower Bound: {}".format(self.lower_spin.value()))
        if self.upper_check.isChecked():
            clips.append("Upper Bound: {}".format(self.upper_spin.value()))
        return ", ".join(clips) if clips else "No Clipping"


class SelectGenesEditor(ScBaseEditor):
    DEFAULT_N_GENS = 1000
    DEFAULT_METHOD = SelectMostVariableGenes.Dispersion
    DEFAULT_COMPUTE_STATS = True
    DEFAULT_N_GROUPS = 20

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.setLayout(QVBoxLayout())
        self._n_genes = self.DEFAULT_N_GENS
        self._n_groups = self.DEFAULT_N_GROUPS

        form = QFormLayout()
        self.n_genes_spin = QSpinBox(minimum=1, maximum=10 ** 6,
                                     value=self._n_genes)
        self.n_genes_spin.valueChanged[int].connect(self._set_n_genes)
        self.n_genes_spin.editingFinished.connect(self.edited)
        form.addRow("Number of genes:", self.n_genes_spin)
        self.layout().addLayout(form)

        disp_b = QRadioButton("Dispersion", checked=True)
        vari_b = QRadioButton("Variance")
        mean_b = QRadioButton("Mean")
        self.group = QButtonGroup()
        self.group.buttonClicked.connect(self._on_button_clicked)
        for i, button in enumerate([disp_b, vari_b, mean_b]):
            index = index_to_enum(SelectMostVariableGenes.Method, i).value
            self.group.addButton(button, index - 1)
            form.addRow(button)

        self.stats_check = QCheckBox("Compute statistics for",
                                     checked=self.DEFAULT_COMPUTE_STATS)
        self.stats_check.clicked.connect(self.edited)
        self.n_groups_spin = QSpinBox(minimum=1, value=self._n_groups)
        self.n_groups_spin.valueChanged[int].connect(self._set_n_groups)
        self.n_groups_spin.editingFinished.connect(self.edited)

        box = QHBoxLayout()
        box.addWidget(self.stats_check)
        box.addWidget(self.n_groups_spin)
        box.addWidget(QLabel("gene groups."))
        box.addStretch()
        self.layout().addLayout(box)

    def _set_n_genes(self, n):
        if self._n_genes != n:
            self._n_genes = n
            self.n_genes_spin.setValue(n)
            self.changed.emit()

    def _set_n_groups(self, n):
        if self._n_groups != n:
            self._n_groups = n
            self.n_groups_spin.setValue(n)
            self.changed.emit()

    def _on_button_clicked(self):
        self.changed.emit()
        self.edited.emit()

    def setParameters(self, params):
        self._set_n_genes(params.get("n_genes", self.DEFAULT_N_GENS))
        method = params.get("method", self.DEFAULT_METHOD)
        index = enum_to_index(SelectMostVariableGenes.Method, method)
        self.group.buttons()[index].setChecked(True)
        compute_stats = params.get("compute_stats", self.DEFAULT_COMPUTE_STATS)
        self.stats_check.setChecked(compute_stats)
        self._set_n_groups(params.get("n_groups", self.DEFAULT_N_GROUPS))

    def parameters(self):
        method = index_to_enum(SelectMostVariableGenes.Method,
                               self.group.checkedId())
        return {"n_genes": self._n_genes, "method": method,
                "compute_stats": self.stats_check.isChecked(),
                "n_groups": self._n_groups}

    @staticmethod
    def createinstance(params):
        method = params.get("method", SelectGenesEditor.DEFAULT_METHOD)
        n_genes = params.get("n_genes", SelectGenesEditor.DEFAULT_N_GENS)
        compute_stats = params.get(
            "compute_stats", SelectGenesEditor.DEFAULT_COMPUTE_STATS)
        n_groups = params.get("n_groups", SelectGenesEditor.DEFAULT_N_GROUPS) \
            if compute_stats else None
        return SelectMostVariableGenes(method, n_genes, n_groups)

    def __repr__(self):
        method = self.group.button(self.group.checkedId()).text()
        text = "Method: {}, Number of Genes: {}".format(method, self._n_genes)
        if self.stats_check.isChecked():
            text += ", Number of Groups: {}".format(self._n_groups)
        return text


class DropoutEditor(ScBaseEditor):
    DEFAULT_N_GENES = 1000

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.setLayout(QVBoxLayout())
        self._n_genes = self.DEFAULT_N_GENES

        form = QFormLayout()
        self.n_genes_spin = QSpinBox(minimum=1, maximum=10 ** 6,
                                     value=self._n_genes)
        self.n_genes_spin.valueChanged[int].connect(self._set_n_genes)
        self.n_genes_spin.editingFinished.connect(self.edited)
        form.addRow("Number of genes:", self.n_genes_spin)
        self.layout().addLayout(form)

    def _set_n_genes(self, n):
        if self._n_genes != n:
            self._n_genes = n
            self.n_genes_spin.setValue(n)
            self.changed.emit()

    def setParameters(self, params):
        self._set_n_genes(params.get("n_genes", self.DEFAULT_N_GENES))

    def parameters(self):
        return {"n_genes": self._n_genes}

    @staticmethod
    def createinstance(params):
        n_genes = params.get("n_genes", DropoutEditor.DEFAULT_N_GENES)
        return DropoutGeneSelection(n_genes)

    def __repr__(self):
        return "Number of Genes: {}".format(self._n_genes)


PREPROCESS_ACTIONS = [
    PreprocessAction(
        "Logarithmic Scale", "preprocess.log_scale", "Value-Based",
        Description("Logarithmic Scale",
                    icon_path("LogarithmicScale.svg")),
        LogarithmicScaleEditor
    ),
    PreprocessAction(
        "Binarize Expression", "preprocess.binarize", "Value-Based",
        Description("Binarize Expression",
                    icon_path("Binarize.svg")),
        BinarizeEditor
    ),
    PreprocessAction(
        "Normalize Samples", "preprocess.normalize", "Row-Based",
        Description("Normalize Samples",
                    icon_path("Normalize.svg")),
        NormalizeEditor
    ),
    PreprocessAction(
        "Standardize Genes", "preprocess.standardize", "Column-Based",
        Description("Standardize Genes",
                    icon_path("Standardize.svg")),
        StandardizeEditor
    ),
    PreprocessAction(
        "Select Most Variable Genes", "preprocess.select_genes",
        "Column-Based",
        Description("Select Most Variable Genes",
                    icon_path("SelectGenes.svg")),
        SelectGenesEditor
    ),
    PreprocessAction(
        "Dropout Gene Selection", "preprocess.dropout", "Column-Based",
        Description("Dropout Gene Selection", icon_path("Dropout.svg")),
        DropoutEditor
    )
]


class ScController(Controller):
    def __init__(self, view, model=None, parent=None):
        super().__init__(view, model, parent)
        self._master = parent

    def createWidgetFor(self, index):
        definition = index.data(DescriptionRole)
        widget = definition.viewclass(master=self._master)
        return widget


class OWscPreprocess(Orange.widgets.data.owpreprocess.OWPreprocess):
    name = "Single Cell Preprocess"
    description = "Preprocess Single Cell data set"
    icon = "icons/SingleCellPreprocess.svg"
    priority = 220

    class Inputs:
        data = Input("Data", Orange.data.Table)

    class Outputs:
        preprocessed_data = Output("Preprocessed Data", Orange.data.Table)

    class Error(Orange.widgets.data.owpreprocess.OWPreprocess.Error):
        unknown_error = Msg("{}")
        discrete_attributes = Msg("Data with discrete attributes "
                                  "can not be preprocessed.")

    class Warning(Orange.widgets.data.owpreprocess.OWPreprocess.Warning):
        missing_values = Msg("Missing values have been replaced with 0.")
        dropout_warning = Msg("{}")
        bad_pp_combination = Msg("Dropout should not be used after "
                                 "Normalization/Log-transformation.")

    PREPROCESSORS = PREPROCESS_ACTIONS
    DEFAULT_PP = {"preprocessors": [("preprocess.normalize", {}),
                                    ("preprocess.log_scale", {}),
                                    ("preprocess.select_genes", {}),
                                    ("preprocess.standardize", {})]}
    CONTROLLER = ScController
    storedsettings = Setting(DEFAULT_PP)
    group_var = ContextSetting(None)
    settingsHandler = DomainContextHandler()

    input_data_changed = Signal()

    @Inputs.data
    @check_sql_input
    def set_data(self, data=None):
        """Set the input dataset."""
        self.closeContext()
        self.data = data
        self.openContext(data)
        self.check_data()
        self.input_data_changed.emit()
        self.load_group_var()

    def check_data(self):
        self.Error.discrete_attributes.clear()
        if self.data and self.data.domain.has_discrete_attributes():
            self.data = None
            self.Error.discrete_attributes()

        self.Warning.missing_values.clear()
        if self.data and np.isnan(self.data.X).any():
            self.data.X = np.nan_to_num(self.data.X)
            self.Warning.missing_values()

    def load_group_var(self):
        for index in range(self.preprocessormodel.rowCount()):
            item = self.preprocessormodel.item(index)
            params = item.data(ParametersRole)
            if "group_var" in params:
                params["group_var"] = self.group_var
                item.setData(params, ParametersRole)

    def save(self, model):
        d = {"name": ""}
        preprocessors = []
        for i in range(model.rowCount()):
            item = model.item(i)
            pp_def = item.data(DescriptionRole)
            params = item.data(ParametersRole)
            group_var = params.get("group_var")
            if group_var is not None:
                self.group_var = group_var
                params = dict(params)
                params["group_var"] = None
            preprocessors.append((pp_def.qualname, params))

        d["preprocessors"] = preprocessors
        return d

    def apply(self):
        self.storeSpecificSettings()
        preprocessor = self.buildpreproc()
        data = None
        self.Warning.dropout_warning.clear()
        self.Warning.bad_pp_combination.clear()
        self.Error.unknown_error.clear()
        if self.data is not None:
            try:
                with warnings.catch_warnings(record=True) as w:
                    data = preprocessor(self.data)
                for warning in w:
                    if issubclass(warning.category, DropoutWarning):
                        self.Warning.dropout_warning(warning.message)
                    else:
                        warnings.warn(warning.message)
                if self.is_bad_combination(preprocessor):
                    self.Warning.bad_pp_combination()
            except (ValueError, ZeroDivisionError) as e:
                self.Error.unknown_error(str(e))
        self.Outputs.preprocessed_data.send(data)

    @staticmethod
    def is_bad_combination(pp) -> bool:
        if not isinstance(pp, PreprocessorList):
            return False

        index_norm, index_log, index_dropout = None, None, None
        for i in range(len(pp.preprocessors)):
            if isinstance(pp.preprocessors[i], NormalizeSamples):
                index_norm = i
            elif isinstance(pp.preprocessors[i], LogarithmicScale):
                index_log = i
            elif isinstance(pp.preprocessors[i], DropoutGeneSelection):
                index_dropout = i
        return index_dropout is not None and \
               (index_norm is not None and index_dropout > index_norm or
                index_log is not None and index_dropout > index_log)


def main(args=None):
    from Orange.data import Table

    app = QApplication(args or [])
    w = OWscPreprocess()
    w.set_data(Table("iris"))
    w.show()
    w.raise_()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()


if __name__ == "__main__":
    sys.exit(main())

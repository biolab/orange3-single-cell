import sys
import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import kneighbors_graph

from AnyQt.QtCore import Qt
from AnyQt.QtGui import QIntValidator
from AnyQt.QtWidgets import QApplication, QComboBox, QLineEdit

from Orange.data import Table, ContinuousVariable
from Orange.widgets import gui
from Orange.widgets.settings import Setting, DomainContextHandler
from Orange.widgets.widget import Input, Output, Msg, OWWidget
from Orange.preprocess import score
from Orange.widgets.utils.concurrent import ConcurrentWidgetMixin, TaskState


def morans_i(x: np.ndarray, adj: sp.spmatrix) -> np.ndarray:
    N = x.shape[0]
    W = adj.sum()

    x_centered = x - x.mean(axis=0)

    n = np.sum(x_centered * (adj.tocsr().dot(x_centered)), axis=0)
    d = np.sum(x_centered**2, axis=0)

    return N / W * n / (d + 1e-16)


def gearys_c(x: np.ndarray, adj: sp.spmatrix) -> np.ndarray:
    adj = adj.tocoo()

    N = x.shape[0]
    W = adj.sum()

    diff = (x[adj.row, :] - x[adj.col, :]) ** 2
    n = adj.data.dot(diff)
    d = np.sum((x - np.mean(x, axis=0)) ** 2, axis=0)

    return (N - 1) / (2 * W) * n / (d + 1e-16)


class SpatialScorer(score.Scorer):
    name = "Spatial Autocorrelation"

    def __init__(self, adjacency_matrix, method="Moran", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.adjacency_matrix = adjacency_matrix
        self.method = method

    def score_data(self, data: Table) -> np.ndarray:
        scores = []
        for i in range(data.X.shape[1]):
            if self.method == "Moran":
                score = morans_i(data.X[:, i].reshape(-1, 1), self.adjacency_matrix)
            else:
                score = gearys_c(data.X[:, i].reshape(-1, 1), self.adjacency_matrix)
            scores.append(
                [score[0]],
            )
        scores = np.array(scores).T

        return scores


class OWSpatialAutocorrelation(OWWidget, ConcurrentWidgetMixin):
    name = "Spatial Autocorrelation Scorer"
    description = "Calculate Moran's I or Geary's C spatial autocorrelation score."
    icon = "icons/SpatialAutocorrelation.svg"
    priority = 420

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        scorer = Output("Scorer", SpatialScorer)

    class Error(OWWidget.Error):
        general_error = Msg({})
        discrete_attributes = Msg(
            "Data with discrete attributes " "can not be processed."
        )

    class Warning(OWWidget.Warning):
        missing_values = Msg("Missing values have been replaced with 0.")
        negative_values = Msg(
            "Unable to use current settings due " "to negative values in data."
        )
        invalid_k = Msg("Invalid value for k. Using default value of 30.")

    resizing_enabled = False
    want_main_area = False

    feature_x = Setting(None)
    feature_y = Setting(None)
    method_name = Setting("Moran")  # Default method is Moran
    k_neighbors = Setting(30)
    auto_commit = Setting(True)

    settingsHandler = DomainContextHandler()

    def __init__(self, parent=None):
        super().__init__(parent)
        ConcurrentWidgetMixin.__init__(self)
        self.data = None
        self.adjacency_matrix = None

        # Info
        infobox = gui.widgetBox(self.controlArea, "Info")
        self.info_label = gui.widgetLabel(infobox, "No data on input.")

        # Variables box
        variables_box = gui.widgetBox(self.controlArea, "Variables")
        
        # Combo boxes for feature selection
        self.feature_x_combo = QComboBox(self)
        self.feature_y_combo = QComboBox(self)
        self.feature_x_combo.currentIndexChanged.connect(self._on_feature_x_changed)
        self.feature_y_combo.currentIndexChanged.connect(self._on_feature_y_changed)

        gui.widgetLabel(variables_box, "Select X:")
        variables_box.layout().addWidget(self.feature_x_combo)
        gui.widgetLabel(variables_box, "Select Y:")
        variables_box.layout().addWidget(self.feature_y_combo)

        # Method box
        method_box = gui.widgetBox(self.controlArea, "Method")
        method_box.layout().setAlignment(Qt.AlignTop)
        
        method_combo = gui.comboBox(
            method_box,
            self,
            "method_name",
            items=["Moran I", "Geary C"],
            callback=self._on_method_changed,
        )
        method_combo.setCurrentText(
            self.method_name
        )  # Set current selection to the saved method

        # Input for k-neighbors
        self.k_input = QLineEdit(self)
        self.k_input.setText(str(self.k_neighbors))
        self.k_input.setValidator(QIntValidator(1, 100, self))
        self.k_input.editingFinished.connect(self._on_k_changed)
        gui.widgetLabel(method_box, "Radius (k):")
        method_box.layout().addWidget(self.k_input)

        # Auto commit
        gui.auto_commit(
            self.controlArea, self, "auto_commit", "Apply", "Apply Automatically"
        )

    def _on_feature_x_changed(self):
        self.feature_x = self.feature_x_combo.currentText()
        self.commit()

    def _on_feature_y_changed(self):
        self.feature_y = self.feature_y_combo.currentText()
        self.commit()

    def _on_method_changed(self):
        self.commit()

    def _on_k_changed(self):
        try:
            self.k_neighbors = int(self.k_input.text())
        except ValueError:
            self.k_neighbors = 30  # Default value if conversion fails
            self.Warning.invalid_k()
        self.commit()

    @Inputs.data
    def set_data(self, data):
        self.closeContext()
        self.clear()
        self.data = data
        self._setup_info_label()
        self._check_data()
        self.openContext(data)
        if self.data is not None:
            self._populate_combos()
        self.commit()

    def clear(self):
        self.feature_x_combo.clear()
        self.feature_y_combo.clear()
        self.adjacency_matrix = None

    def _setup_info_label(self):
        text = "No data on input."
        if self.data is not None:
            domain, attrs = self.data.domain, self.data.domain.attributes
            text = "{} cells, {} genes\n".format(len(self.data), len(attrs))
            text += (
                "{} meta features".format(len(domain.metas))
                if len(domain.metas)
                else "(no meta features)"
            )
        self.info_label.setText(text)

    def _check_data(self):
        self.clear_messages()
        if self.data and self.data.domain.has_discrete_attributes():
            self.data = None
            self.Error.discrete_attributes()
        if self.data is not None and np.isnan(self.data.X).any():
            # copy data to not to edit input data
            self.data = self.data.copy()
            with self.data.unlocked_reference(self.data.X):
                self.data.X = np.nan_to_num(self.data.X)
            self.Warning.missing_values()

    def _populate_combos(self):
        domain = self.data.domain
        continuous_attrs = [
            attr
            for attr in domain.metas + domain.attributes
            if isinstance(attr, ContinuousVariable)
        ]
        self.feature_x_combo.addItems([attr.name for attr in continuous_attrs])
        self.feature_y_combo.addItems([attr.name for attr in continuous_attrs])
        self.feature_x_combo.setCurrentIndex(0)
        self.feature_y_combo.setCurrentIndex(1)

    def calculate(self, state: TaskState):
        if self.data is None:
            self.Error.general_error("No data to process.")
            return
        x_coords = np.nan_to_num(self.data.get_column(self.feature_x))
        y_coords = np.nan_to_num(self.data.get_column(self.feature_y))
        graph = np.column_stack((x_coords, y_coords))
        state.set_status("Calculating adjacency matrix...")
        self.adjacency_matrix = kneighbors_graph(
            graph, self.k_neighbors, mode="connectivity", include_self=False
        )

    def on_exception(self, ex):
        self.Error.general_error(str(ex))

    def on_done(self, result):
        scorer = SpatialScorer(self.adjacency_matrix, method=self.method_name)
        self.Outputs.scorer.send(scorer)

    def commit(self):
        self.Error.general_error.clear()
        self.Warning.negative_values.clear()
        self.Warning.missing_values.clear()
        self._calculate_task()

    def _calculate_task(self):
        self.start(self.calculate)


def main(args=None):
    app = QApplication(args or [])
    w = OWSpatialAutocorrelation()
    w.set_data(Table("iris"))
    w.show()
    w.raise_()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()


if __name__ == "__main__":
    sys.exit(main())

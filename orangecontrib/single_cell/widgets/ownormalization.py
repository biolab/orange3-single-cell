from AnyQt.QtCore import Qt, QTimer

from Orange.data import Table, DiscreteVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.preprocess.preprocess import Preprocess

from orangecontrib.single_cell.preprocess.scnormalize import SCNormalizer
from orangecontrib.single_cell.preprocess.scnorm import ScNormPreprocessor


class OWNormalization(widget.OWWidget):
    name = 'Single Cell Normalization'
    description = ('Normalization of single cell count data')
    icon = 'icons/Normalization.svg'
    priority = 110

    # Default grouping model
    GROUPING_CELL_GROUP = "(One group per cell)"

    # Normalization models
    MODEL_NONE = 0
    MODEL_CELL_MEDIAN = 1
    MODEL_SCNORM = 2

    inputs = [("Data", Table, 'set_data')]
    outputs = [("Data", Table), ("Preprocessor", Preprocess)]

    want_main_area = False
    resizing_enabled = False

    settingsHandler = settings.DomainContextHandler(metas_in_res=True)
    selected_attr = settings.Setting(GROUPING_CELL_GROUP)
    autocommit = settings.Setting(True)

    model_select = settings.Setting(MODEL_CELL_MEDIAN)
    log_check = settings.Setting(True)
    log_base = settings.Setting(2)

    # ScNorm-specific settings
    scnorm_num_groups = settings.Setting(10)
    scnorm_p_subgroup = settings.Setting(5)

    def __init__(self):
        self.data = None
        self.info = gui.label(self.controlArea, self,
                              "No data on input", box="Info")

        # Library / group variable
        box0 = gui.vBox(self.controlArea, "Normalization model")
        box0.layout().setSpacing(7)
        model_selection = gui.radioButtons(
            widget=box0,
            master=self,
            value="model_select",
            callback=self.on_changed,
            btnLabels=("None",
                       "Equalize group median expression",
                       "Gene group quantile regression (scNorm)"))
        model_selection.layout().setSpacing(7)

        # Cell grouping
        attrs_model = self.attrs_model = DomainModel(
            placeholder=self.GROUPING_CELL_GROUP,
            order=(DomainModel.CLASSES, DomainModel.METAS),
            valid_types=DiscreteVariable)
        combo_attrs = self.combo_attrs = gui.comboBox(
            box0, self, 'selected_attr',
            label="Group by:",
            callback=self.on_changed,
            sendSelectedValue=True)
        combo_attrs.setModel(attrs_model)

        # ScNorm specific
        # Steps and parameters
        self.scnorm_spin_n = gui.spin(box0, self, "scnorm_num_groups", 1, 20, step=1, label="Max. gene groups: ",
                                        alignment=Qt.AlignRight, callback=self.on_changed,
                                        checkCallback=self.on_changed, controlWidth=60)
        self.scnorm_spin_p = gui.spin(box0, self, "scnorm_p_subgroup", 0, 100, step=1, label="Prototype genes (%): ",
                                      tooltip="Fraction of prototype genes to fit expression in each group. " +
                                              "Smaller values imply faster execution, but lower accuracy",
                                      alignment=Qt.AlignRight, callback=self.on_changed,
                                      checkCallback=self.on_changed, controlWidth=60)

        # Steps and parameters
        box1 = gui.widgetBox(self.controlArea, 'Further steps and parameters')
        gui.spin(box1, self, "log_base", 2.0, 10.0, label="Log(1 + x) transform, base: ",
                 checked="log_check", alignment=Qt.AlignRight,
                 callback=self.on_changed,
                 checkCallback=self.on_changed, controlWidth=60)

        gui.auto_commit(self.controlArea, self, 'autocommit', '&Apply')
        QTimer.singleShot(0, self.commit)

    def set_data(self, data):
        self.closeContext()
        self.data = data

        if self.data is None:
            self.attrs_model.clear()
            self.commit()
            self.info.setText("No data on input")
            return

        self.info.setText("%d cells, %d features." %
                          (len(data), len(data.domain.attributes)))

        self.attrs_model.set_domain(data.domain)

        self.send("Data", None)
        self.openContext(self.data.domain)
        self.on_changed()

    def on_changed(self):
        self.scnorm_spin_n.setEnabled(self.model_select == self.MODEL_SCNORM)
        self.scnorm_spin_p.setEnabled(self.model_select == self.MODEL_SCNORM)
        self.commit()

    def commit(self):
        log_base = self.log_base if self.log_check else None
        library_var = None
        if self.data is not None and self.selected_attr in self.data.domain:
            library_var = self.data.domain[self.selected_attr]

        if self.model_select == self.MODEL_NONE:
            pp = SCNormalizer(equalize_var=None,
                              normalize_cells=False,
                              log_base=log_base)
        elif self.model_select == self.MODEL_CELL_MEDIAN:
            pp = SCNormalizer(equalize_var=library_var,
                              normalize_cells=True,
                              log_base=log_base)
        elif self.model_select == self.MODEL_SCNORM:
            pp = ScNormPreprocessor(equalize_var=library_var,
                                    p_subgroup=self.scnorm_p_subgroup / 100.0,
                                    K=self.scnorm_p_subgroup,
                                    log_base=log_base)

        data = None
        if self.data is not None:
            data = pp(self.data)

        self.send("Data", data)
        self.send("Preprocessor", pp)


if __name__ == "__main__":
    from sys import argv
    from AnyQt.QtWidgets import QApplication

    app = QApplication([])
    ow = OWNormalization()

    # Load test file from arguments
    test_file = argv[1] if len(argv) >= 2 else "matrix_counts_sample.tab"
    table = Table(test_file)[:100, :100]
    ow.set_data(table)

    ow.show()
    app.exec()
    ow.saveSettings()
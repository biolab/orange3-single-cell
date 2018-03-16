from AnyQt.QtCore import Qt, QTimer
from AnyQt.QtWidgets import QListView

from Orange.data import Table, DiscreteVariable, ContinuousVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel, PyListModel
from Orange.widgets.widget import Input, Output
from Orange.preprocess.preprocess import PreprocessorList, Preprocess

from orangecontrib.single_cell.preprocess.scnormalize import SCNormalizer
from orangecontrib.single_cell.preprocess.scbnorm import LINKS, LINK_IDENTITY, LINK_LOG, \
    SCBatchNormalizer


class OWNormalization(widget.OWWidget):
    name = 'Normalize'
    description = 'Basic normalization of single cell count data'
    icon = 'icons/Normalization.svg'
    priority = 160

    DEFAULT_CELL_NORM = "(One group per cell)"

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)
        preprocessor = Output("Preprocessor", Preprocess)

    want_main_area = False
    resizing_enabled = False

    settingsHandler = settings.DomainContextHandler()
    selected_attr = settings.Setting(DEFAULT_CELL_NORM)
    autocommit = settings.Setting(True)

    normalize_cells = settings.Setting(True)
    log_check = settings.Setting(True)
    log_base = 2

    batch_link = LINK_IDENTITY
    batch_vars = ()

    def __init__(self):
        self.data = None
        self.info = gui.label(self.controlArea, self,
                              "No data on input", box="Info")

        # Library / group variable
        box0 = gui.vBox(
            self.controlArea, "Data from multiple libraries")
        self.normalize_check = gui.checkBox(box0,
                                self, "normalize_cells",
                                "Normalize cell profiles on:",
                                callback=self.on_changed_normalize,
                                addSpace=True)
        attrs_model = self.attrs_model = DomainModel(
            placeholder=self.DEFAULT_CELL_NORM,
            order=(DomainModel.CLASSES, DomainModel.METAS),
            valid_types=DiscreteVariable)
        combo_attrs = self.combo_attrs = gui.comboBox(
            box0, self, 'selected_attr',
            callback=self.on_changed,
            sendSelectedValue=True)
        combo_attrs.setModel(attrs_model)

        # Steps and parameters
        box1 = gui.widgetBox(self.controlArea, 'Further steps and parameters')
        gui.spin(box1, self, "log_base", 2.0, 10.0, label="Log(1 + x) transform, base: ",
                 checked="log_check", alignment=Qt.AlignRight,
                 callback=self.on_changed,
                 checkCallback=self.on_changed, controlWidth=60)

        # Batch effects
        box2 = gui.vBox(self.controlArea, "Variables to regress out (batch effects)")
        self.batch_link_combo = gui.comboBox(
            box2, self, 'batch_link',
            callback=self._batch_changed,
            sendSelectedValue=True)
        self.batch_link_combo.setModel(PyListModel(LINKS.keys()))

        self.batch_attrs = DomainModel(order=(DomainModel.CLASSES, DomainModel.METAS),
                                       valid_types=(DiscreteVariable, ContinuousVariable))
        self.varview = view = QListView(selectionMode=QListView.MultiSelection)
        view.setModel(self.batch_attrs)
        view.selectionModel().selectionChanged.connect(self._batch_changed)
        box2.layout().addWidget(view)

        gui.auto_commit(self.controlArea, self, 'autocommit', '&Apply')
        QTimer.singleShot(0, self.commit)

    @Inputs.data
    def set_data(self, data):
        self.closeContext()
        self.data = data

        if self.data is None:
            self.attrs_model.set_domain(None)
            self.batch_attrs.set_domain(None)
            self.commit()
            self.info.setText("No data on input")
            return

        self.info.setText("%d cells, %d features." %
                          (len(data), len(data.domain.attributes)))

        self.attrs_model.set_domain(data.domain)
        self.batch_attrs.set_domain(data.domain)
        self.normalize_check.setEnabled(len(self.attrs_model) > 0)
        self.combo_attrs.setEnabled(self.normalize_cells)

        self.Outputs.data.send(None)
        self.openContext(self.data.domain)
        self.on_changed()

    def _batch_changed(self):
        self.batch_vars = tuple(self.batch_attrs[ind.row()].name
                                for ind in self.varview.selectionModel().selectedRows())

    def on_changed(self):
        self.commit()

    def on_changed_normalize(self):
        self.combo_attrs.setEnabled(self.normalize_cells)
        self.commit()

    def commit(self):
        log_base = self.log_base if self.log_check else None
        library_var = None
        if self.data is not None and \
                self.normalize_cells and \
                self.selected_attr in self.data.domain:
            library_var = self.data.domain[self.selected_attr]

        pp = SCNormalizer(equalize_var=library_var,
                          normalize_cells=self.normalize_cells,
                          log_base=log_base)

        pp_batch = SCBatchNormalizer(link=self.batch_link,
                                     nonzero_only=self.batch_link == LINK_LOG,
                                     batch_vars=self.batch_vars)

        data = None
        if self.data is not None:
            data = pp_batch(pp(self.data))

        self.Outputs.data.send(data)
        self.Outputs.preprocessor.send(PreprocessorList([pp, pp_batch]))


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication

    app = QApplication([])
    ow = OWNormalization()

    # Load test file from arguments
    table = Table("iris")
    ow.set_data(table)

    ow.show()
    app.exec()
    ow.saveSettings()

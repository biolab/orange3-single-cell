import numpy as np
import logging

from AnyQt.QtCore import Qt, QItemSelection, QItemSelectionRange, QItemSelectionModel
from Orange.data import Table, DiscreteVariable, ContinuousVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.widget import Input, Output
from Orange.preprocess.preprocess import PreprocessorList, Preprocess


from orangecontrib.single_cell.preprocess.scnormalize import SCNormalizer
from orangecontrib.single_cell.widgets.owscoregenes import TableModel, TableView

from orangecontrib.single_cell.preprocess.scbnorm import LINKS, LINK_LOG, SCBatchNormalizer, ScBatchScorer

log = logging.getLogger(__name__)


class OWNormalization(widget.OWWidget):
    name = 'Normalize'
    description = 'Normalization of single cell count data'
    icon = 'icons/Normalization.svg'
    priority = 160

    DEFAULT_CELL_NORM = "(One group per cell)"
    SCORERS = (ScBatchScorer, )
    LINK_FUNCTIONS = sorted(LINKS.keys())

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)
        preprocessor = Output("Preprocessor", Preprocess)

    # Widget settings
    want_main_area = False
    resizing_enabled = False
    settingsHandler = settings.DomainContextHandler()
    autocommit = settings.Setting(True, schema_only=True)

    # Settings (basic preprocessor)
    normalize_cells = settings.Setting(True, schema_only=True)
    selected_attr_index = settings.Setting(0, schema_only=True)
    log_check = settings.Setting(True, schema_only=True)
    log_base = settings.Setting(2, schema_only=True)

    # Settings (batch preprocessor)
    batch_link_index = settings.Setting(0, schema_only=True)
    batch_vars_selected = settings.Setting([], schema_only=True)
    batch_vars_all = []
    batch_vars_names = []

    # Preprocessors
    pp = None
    pp_batch = None

    # ranksView global settings
    sorting = settings.Setting((0, Qt.DescendingOrder))

    def __init__(self):
        self.data = None
        self.info = gui.label(self.controlArea, self,
                              "No data on input", box="Info")

        # Library / group variable
        box0 = gui.vBox(
            self.controlArea, "Data from multiple libraries")
        self.normalize_check = gui.checkBox(box0,
                                self, "normalize_cells",
                                "Normalize median expressions on cell groups:",
                                callback=self.on_changed,
                                addSpace=True)

        self.attrs_model = DomainModel(
            placeholder=self.DEFAULT_CELL_NORM,
            order=(DomainModel.CLASSES, DomainModel.METAS),
            valid_types=DiscreteVariable)

        self.combo_attrs = gui.comboBox(
            box0, self, 'selected_attr_index',
            callback=self.on_changed)

        # Steps and parameters
        box1 = gui.widgetBox(self.controlArea, 'Further steps and parameters')
        gui.spin(box1, self, "log_base", 2.0, 10.0, label="Log(1 + x) transform, base: ",
                 checked="log_check", alignment=Qt.AlignRight,
                 callback=self.on_changed,
                 checkCallback=self.on_changed, controlWidth=60)

        # Batch effects - link function
        box2 = gui.vBox(self.controlArea, "Variables to regress out (batch effects)")
        self.batch_link_combo = gui.comboBox(
            box2, self, 'batch_link_index',
            callback=self.on_changed,
            items=self.LINK_FUNCTIONS)

        # Batch effects - variables
        self.ranksModel = model = TableModel(parent=self)  # type: TableModel
        self.ranksView = view = TableView(self)            # type: TableView
        box2.layout().addWidget(view)
        view.setModel(model)
        view.setColumnWidth(0, 30)
        view.setSelectionMode(TableView.MultiSelection)
        view.selectionModel().selectionChanged.connect(self.on_select)
        view.horizontalHeader().sectionClicked.connect(self.on_header_click)

        # Autocommit
        gui.auto_commit(self.controlArea, self, 'autocommit', '&Apply')

    ### Called once at the arrival of new data ###

    @Inputs.data
    def set_data(self, data):
        self.data = data

        if self.data is None:
            self.selected_attr_index = 0
            self.attrs_model.set_domain(None)

            self.batch_vars_selected.clear()
            self.ranksModel.clear()
            self.ranksModel.resetSorting(True)
            self.info.setText("No data on input")
            self.commit()
            return

        self.info.setText("%d cells, %d features." %
                          (len(data), len(data.domain.attributes)))


        self.attrs_model.set_domain(data.domain)
        self.normalize_check.setEnabled(len(self.attrs_model) > 0)
        self.combo_attrs.setEnabled(self.normalize_cells)
        self.combo_attrs.setModel(self.attrs_model)
        self.set_batch_variables()

        # Implicit commit
        self.update_state()

    def set_batch_variables(self):
        """ Search for meta variables and classes in new data. """
        self.batch_vars_all = [var
                               for var in self.data.domain.metas + self.data.domain.class_vars
                               if isinstance(var, ContinuousVariable) or isinstance(var, DiscreteVariable)]
        self.batch_vars_names = tuple(a.name for a in self.batch_vars_all)
        if self.data is not None and len(self.batch_vars_all) == 0:
            return

        self.ranksModel.setVerticalHeaderLabels(self.batch_vars_all)
        self.ranksView.setVHeaderFixedWidthFromLabel(max(self.batch_vars_names, key=len))

    ### Event handlers ###

    def on_changed(self):
        """ Update graphics, model parameters and commit. """
        self.combo_attrs.setEnabled(self.normalize_cells)
        self.update_state()
        self.commit()

    def on_select(self):
        """ Save indices of attributes in the original, unsorted domain.
            Warning: this method must not call update_scores; """
        selected_rows = self.ranksModel.mapToSourceRows([
            i.row() for i in self.ranksView.selectionModel().selectedRows(0)])
        self.batch_vars_selected.clear()
        self.batch_vars_selected.extend([self.batch_vars_all[i].name for i in selected_rows])
        self.update_preprocessors()
        self.commit()

    def on_header_click(self):
        """ Store the header states. """
        sort_order = self.ranksModel.sortOrder()
        sort_column = self.ranksModel.sortColumn() - 1  # -1 for '#' (discrete count) column
        self.sorting = (sort_column, sort_order)

    ### Updates to model parameters / scores ###

    def update_state(self):
        """ Updates preprocessors and scores in the correct order. """
        self.update_preprocessors()
        self.update_scores()

    def update_selection(self):
        """ Update selected rows (potentially from loaded scheme) at once."""
        sel_model = self.ranksView.selectionModel()
        ncol = self.ranksModel.columnCount()
        model = self.ranksModel
        selection = QItemSelection()
        selected_rows = [self.batch_vars_names.index(b)
                         for b in self.batch_vars_selected.copy()
                         if b in self.batch_vars_names]
        if len(selected_rows):
            for row in model.mapFromSourceRows(selected_rows):
                selection.append(QItemSelectionRange(
                    model.index(row, 0), model.index(row, ncol - 1)))
            sel_model.select(selection, QItemSelectionModel.ClearAndSelect)
        else:
            self.commit()

    def update_preprocessors(self):
        """ Update parameters of processors. """
        log_base = self.log_base if self.log_check else None
        library_var = None
        selected_attr = self.attrs_model[self.selected_attr_index]
        batch_link = self.LINK_FUNCTIONS[self.batch_link_index]

        if self.data is not None and \
                self.normalize_cells and \
                selected_attr in self.data.domain:
            library_var = self.data.domain[selected_attr]

        self.pp = SCNormalizer(equalize_var=library_var,
                               normalize_cells=self.normalize_cells,
                               log_base=log_base)

        self.pp_batch = SCBatchNormalizer(link=batch_link,
                                          nonzero_only=batch_link == LINK_LOG,
                                          batch_vars=self.batch_vars_selected)

    def update_scores(self):
        """ Update scores for current data and preprocessors. """
        if self.data is None:
            self.ranksModel.clear()
            return

        method_scores = tuple(self.get_method_scores(method)
                              for method in self.SCORERS)

        labels = tuple(method.friendly_name for method in self.SCORERS)

        model_array = np.column_stack(
            ([len(a.values) if a.is_discrete else np.nan
              for a in self.batch_vars_all],) +
            (method_scores if method_scores else ())
        )

        # Set fixed extreme values
        self.ranksModel.setExtremesFrom(column=1, values=[0, 1])

        # Update, but retain previous selection
        self.ranksModel.wrap(model_array.tolist())
        self.ranksModel.setHorizontalHeaderLabels(('#',) + labels)
        self.ranksView.setColumnWidth(0, 40)

        # Rows must be reselected again as ranksModel.wrap resets the selection
        self.update_selection()

        # Re-apply sort
        try:
            sort_column, sort_order = self.sorting
            if sort_column < len(labels):
                self.ranksModel.sort(sort_column + 1, sort_order)  # +1 for '#' (discrete count) column
                self.ranksView.horizontalHeader().setSortIndicator(sort_column + 1, sort_order)
        except ValueError:
            pass

    def get_method_scores(self, method):
        """ Compute scores for all batch variables.
            Scores must be computed after applying the first pre-processor. """
        assert self.pp is not None
        estimator = method()
        data = self.pp(self.data)
        try:
            scores = np.array([estimator.score_data(data=data, feature=attr)
                               for attr in self.batch_vars_all])
        except ValueError:
            log.error(
                "Scorer %s wasn't able to compute scores at all",
                estimator.name)
            scores = np.full(len(self.batch_vars_all), np.nan)
        return scores

    ### Set output ###

    def commit(self):
        """ Update parameters to preprocessors and set output signals. """
        data = None
        if self.data is not None:
            data = self.pp_batch(self.pp(self.data))

        self.Outputs.data.send(data)
        self.Outputs.preprocessor.send(PreprocessorList([self.pp, self.pp_batch]))


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

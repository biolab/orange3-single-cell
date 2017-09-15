import numpy as np
import scipy.sparse as sp

from AnyQt.QtCore import Qt
from Orange.data import Table, DiscreteVariable
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import DomainModel


def normalize(X, y=None, equalize=True, normalize_cells=True, log_base=2):
    """
    A simple ad-hoc normalization to provide basic facilities.
    :param X: (Sparse) matrix with expression values as counts.
            Columns are genes and rows are cells.
    :param y: Iterable with unique values for each library (default: None).
        If None, data is assumed to be from a single library.
    :param equalize: Equalize libraries (default: True)
    :param normalize_cells: Normalize w.r.t cells
    :param log_base: log(1 + x) transform.
    :return:
    """
    # Result in expected number of reads
    Xeq = X.copy()

    # Equalize based on read depth per library / match mean read count per cell
    if equalize and y is not None:
        lib_sizes = dict()
        libraries = dict([(lib, np.where(y == lib)[0]) for lib in set(y)])
        for lib, inxs in sorted(libraries.items()):
            lib_sizes[lib] = X[inxs, :].sum(axis=1).mean()
        target = min(lib_sizes.values())
        size_factors = dict([(lib, target/float(size)) for lib, size in lib_sizes.items()])
        for lib, inxs in sorted(libraries.items()):
            Xeq[inxs, :] *= size_factors[lib]

    # Normalize by cells, sweep columns by means / median
    if normalize_cells:
        rs  = np.array(Xeq.sum(axis=1).reshape((Xeq.shape[0], 1)))
        rsm = np.median(rs)
        Xd = sp.dia_matrix(((rsm / rs).ravel(), 0), shape=(len(rs), len(rs)))
        Xeq = Xd.dot(Xeq)

    # Log transform log(1 + x)
    if log_base is not None:
        if sp.isspmatrix(Xeq):
            Xeq = Xeq.log1p() / np.log(log_base)
        else:
            Xeq = np.log(1 + Xeq) / np.log(log_base)

    # Preserve sparsity
    return Xeq.tocsr() if sp.isspmatrix(Xeq) else Xeq


class OWNormalization(widget.OWWidget):
    name = 'Single Cell Normalization'
    description = ('Basic normalization of single cell count data')
    icon = 'icons/Normalization.svg'
    priority = 110

    inputs = [("Data", Table, 'set_data')]
    outputs = [("Data", Table)]

    want_main_area = False
    resizing_enabled = False

    settingsHandler = settings.DomainContextHandler(metas_in_res=True)
    selected_attr = settings.ContextSetting("")
    autocommit = settings.Setting(True)

    equalize_lib = settings.Setting(True)
    normalize_cells = settings.Setting(True)
    log_check = settings.Setting(True)
    log_base = 2

    def __init__(self):
        self.data = None
        self.info = gui.label(self.controlArea, self,
                              "No data on input", box="Info")

        # Library / group variable
        box0 = gui.vBox(
            self.controlArea, "Data from multiple libraries")
        self.equalize_check = gui.checkBox(box0,
                     self, "equalize_lib",
                     "Equalize library read depth on:",
                     callback=self.on_changed_equalize,
                     addSpace=True)
        attrs_model = self.attrs_model = DomainModel(
            order=(DomainModel.CLASSES, DomainModel.METAS),
            valid_types=DiscreteVariable)
        combo_attrs = self.combo_attrs = gui.comboBox(
            box0, self, 'selected_attr',
            callback=self.on_changed,
            sendSelectedValue=True)
        combo_attrs.setModel(attrs_model)

        # Steps and parameters
        box1 = gui.widgetBox(self.controlArea, 'Steps and parameters')
        gui.checkBox(gui.vBox(box1), self, "normalize_cells", "Normalize cell profiles",
                     callback=self.on_changed)
        gui.spin(box1, self, "log_base", 2.0, 10.0, label="Log(1 + x) transform, base: ",
                 checked="log_check", alignment=Qt.AlignRight,
                 callback=self.on_changed,
                 checkCallback=self.on_changed, controlWidth=60)

        gui.auto_commit(self.controlArea, self, 'autocommit', '&Apply')

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
        self.equalize_check.setEnabled(len(self.attrs_model) > 0)
        self.combo_attrs.setEnabled(self.equalize_lib)
        if len(self.attrs_model) != 0:
            self.selected_attr = str(self.attrs_model[0])

        self.send("Data", None)
        self.openContext(self.data.domain)
        self.on_changed()

    def on_changed(self):
        self.commit()

    def on_changed_equalize(self):
        self.combo_attrs.setEnabled(self.equalize_lib)
        self.commit()

    def commit(self):
        if self.data is None:
            self.send("Data", None)
            return

        if self.equalize_lib and self.selected_attr in self.data.domain:
            library_var = self.data.domain[self.selected_attr]
            library_values, _ = self.data.get_column_view(library_var)
        else:
            library_values = None
        log_base = self.log_base if self.log_check else None
        X_normed = normalize(self.data.X,
                             y=library_values,
                             equalize=self.equalize_lib,
                             normalize_cells=self.normalize_cells,
                             log_base=log_base)

        new_data = Table.from_table(source=self.data, domain=self.data.domain)
        new_data.X = X_normed
        self.send("Data", new_data)


if __name__ == "__main__":
    from sys import argv
    from AnyQt.QtWidgets import QApplication

    app = QApplication([])
    ow = OWNormalization()

    # Load test file from arguments
    test_file = argv[1] if len(argv) >= 2 else "matrix_counts_sample.tab"
    table = Table(test_file)
    ow.set_data(table)

    ow.show()
    app.exec()
    ow.saveSettings()
from math import isnan, isinf

import unicodedata
from AnyQt.QtCore import Qt, QItemSelection, QItemSelectionModel
from AnyQt.QtGui import QStandardItem, QColor, QFont, QBrush
from AnyQt.QtWidgets import QTableView, QSizePolicy, QHeaderView, QStyledItemDelegate
from Orange.widgets import gui


BorderRole = next(gui.OrangeUserRole)
BorderColorRole = next(gui.OrangeUserRole)


class BorderedItemDelegate(QStyledItemDelegate):
    """Item delegate that paints border at the specified sides

    Data for `BorderRole` is a string containing letters t, r, b and/or l,
    which defines the sides at which the border is drawn.

    Role `BorderColorRole` sets the color for the cell. If not color is given,
    `self.color` is used as default.

    Args:
        color (QColor): default color (default default is black)
    """
    def __init__(self, color=Qt.black):
        super().__init__()
        self.color = color

    def paint(self, painter, option, index):
        """Overloads `paint` to draw borders"""
        QStyledItemDelegate.paint(self, painter, option, index)
        borders = index.data(BorderRole)
        if borders:
            color = index.data(BorderColorRole) or self.color
            painter.save()
            painter.setPen(color)
            rect = option.rect
            for side, p1, p2 in (("t", rect.topLeft(), rect.topRight()),
                                 ("r", rect.topRight(), rect.bottomRight()),
                                 ("b", rect.bottomLeft(), rect.bottomRight()),
                                 ("l", rect.topLeft(), rect.bottomLeft())):
                if side in borders:
                    painter.drawLine(p1, p2)
            painter.restore()


class ContingencyTable(QTableView):
    """
    A contingency table widget which can be used wherever ``QTableView`` could be used.

    Parameters
    ----------
    parent : Orange.widgets.widget.OWWidget as a
        The containing widget to which the table is connected.
    tablemodel : AnyQt.QtGui.QtStandardItemModel
        The model in which table data will be stored.

    Attributes
    ----------
    classesv : :obj:`list` of :obj:`str`
        Vertical class headers.
    classesh : :obj:`list` of :obj:`str`
        Horizontal class headers.
    headerv : :obj:`str`, optional
        Vertical top header.
    headerh : :obj:`str`, optional
        Horizontal top header.
    corner_string : str
        String that is top right and bottom left corner of the table.
        Default is ``unicodedata.lookup("N-ARY SUMMATION")``.
    """

    def __init__(self, parent, tablemodel):
        super().__init__(editTriggers=QTableView.NoEditTriggers)

        self.classesv = None
        self.classesh = None
        self.headerv = None
        self.headerh = None
        self.parent = parent
        self.tablemodel = tablemodel

        self.corner_string = unicodedata.lookup("N-ARY SUMMATION")

        self.setModel(self.tablemodel)
        self.horizontalHeader().hide()
        self.verticalHeader().hide()
        self.horizontalHeader().setMinimumSectionSize(60)
        self.selectionModel().selectionChanged.connect(self.parent._invalidate)
        self.setShowGrid(False)
        self.setItemDelegate(BorderedItemDelegate(Qt.white))
        self.setSizePolicy(QSizePolicy.MinimumExpanding,
                           QSizePolicy.MinimumExpanding)
        self.clicked.connect(self._cell_clicked)

    def _cell_clicked(self, model_index):
        """Handle cell click event"""
        i, j = model_index.row(), model_index.column()
        if not i or not j:
            return
        n = self.tablemodel.rowCount()
        m = self.tablemodel.columnCount()
        index = self.tablemodel.index
        selection = None
        if i == j == 1 or i == n - 1 and j == m - 1:
            selection = QItemSelection(index(2, 2), index(n - 1, m - 1))
        elif i in (1, n - 1):
            selection = QItemSelection(index(2, j), index(n - 1, j))
        elif j in (1, m - 1):
            selection = QItemSelection(index(i, 2), index(i, m - 1))

        if selection is not None:
            self.selectionModel().select(
                selection, QItemSelectionModel.ClearAndSelect)

    def _item(self, i, j):
        return self.tablemodel.item(i, j) or QStandardItem()

    def _set_item(self, i, j, item):
        self.tablemodel.setItem(i, j, item)

    def set_variables(self, variablev, variableh):
        """
        Sets class headers and top headers and initializes table structure.

        Parameters
        ----------
        variablev : Orange.data.variable.DiscreteVariable
            Class headers are set to ``variablev.values``, top header is set to ``variablev.name``.
        variableh : Orange.data.variable.DiscreteVariable
            Class headers are set to ``variableh.values``, top header is set to ``variableh.name``.
        """
        self.classesv = variablev.values
        self.classesh = variableh.values
        self.headerv = variablev.name
        self.headerh = variableh.name
        self.initialize()

    def set_headers(self, classesv, classesh, headerv=None, headerh=None):
        """
        Sets class headers and top headers and initializes table structure.

        Parameters
        ----------
        classesv : :obj:`list` of :obj:`str`
            Vertical class headers.
        classesh : :obj:`list` of :obj:`str`
            Horizontal class headers.
        headerv : :obj:`str`, optional
            Vertical top header.
        headerh : :obj:`str`, optional
            Horizontal top header.
        """
        self.classesv = classesv
        self.classesh = classesh
        self.headerv = headerv
        self.headerh = headerh
        self.initialize()

    def initialize(self):
        """
        Initializes table structure. Class headers must be set beforehand.
        """
        assert self.classesv is not None and self.classesh is not None

        item = self._item(0, 2)
        item.setData(self.headerh, Qt.DisplayRole)
        item.setTextAlignment(Qt.AlignCenter)
        item.setFlags(Qt.NoItemFlags)

        self._set_item(0, 2, item)
        item = self._item(2, 0)
        item.setData(self.headerv, Qt.DisplayRole)
        item.setTextAlignment(Qt.AlignHCenter | Qt.AlignBottom)
        item.setFlags(Qt.NoItemFlags)
        self.setItemDelegateForColumn(0, gui.VerticalItemDelegate())
        self._set_item(2, 0, item)
        self.setSpan(0, 2, 1, len(self.classesh)+1)
        self.setSpan(2, 0, len(self.classesv)+1, 1)

        font = self.tablemodel.invisibleRootItem().font()
        bold_font = QFont(font)
        bold_font.setBold(True)

        for i in (0, 1):
            for j in (0, 1):
                item = self._item(i, j)
                item.setFlags(Qt.NoItemFlags)
                self._set_item(i, j, item)

        for headers, ix in ((self.classesv + [self.corner_string], lambda p: (p + 2, 1)),
                            (self.classesh + [self.corner_string], lambda p: (1, p + 2))):
            for p, label in enumerate(headers):
                i, j = ix(p)
                item = self._item(i, j)
                item.setData(label, Qt.DisplayRole)
                item.setFont(bold_font)
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                item.setFlags(Qt.ItemIsEnabled)
                if p < len(headers) - 1:
                    item.setData("br"[j == 1], BorderRole)
                    item.setData(QColor(192, 192, 192), BorderColorRole)
                else:
                    item.setData("", BorderRole)
                self._set_item(i, j, item)

        hor_header = self.horizontalHeader()
        if len(' '.join(self.classesh + [self.corner_string])) < 120:
            hor_header.setSectionResizeMode(QHeaderView.ResizeToContents)
        else:
            hor_header.setDefaultSectionSize(60)
        self.tablemodel.setRowCount(len(self.classesv) + 3)
        self.tablemodel.setColumnCount(len(self.classesh) + 3)

    def get_selection(self):
        """
        Get indexes of selected cells.

        Returns
        -------
        :obj:`set` of :obj:`tuple` of :obj:`int`
            Set of pairs of indexes.
        """
        return {(ind.row() - 2, ind.column() - 2) for ind in self.selectedIndexes()}

    def set_selection(self, indexes):
        """
        Set indexes of selected cells.

        Parameters
        ----------
        indexes : :obj:`set` of :obj:`tuple` of :obj:`int`
            Set of pairs of indexes.
        """
        selection = QItemSelection()
        index = self.model().index
        for row, col in indexes:
            sel = index(row + 2, col + 2)
            selection.select(sel, sel)
        self.selectionModel().select(
            selection, QItemSelectionModel.ClearAndSelect)

    def update_table(self, matrix, colsum=None, rowsum=None, colors=None, formatstr="{}", tooltip=None):
        """
        Sets ``matrix`` as data of the table.

        Parameters
        ----------
        matrix : numpy.array
            2D array to be set as data.
        colsum : :obj:`numpy.array`, optional
            1D optional array with aggregate values of columns, defaults to sum.
        rowsum : :obj:`numpy.array`, optional
            1D optional array with aggregate values of rows, defaults to sum.
        colors : :obj:`numpy.array`, optional
            2D array with color values, defaults to no color.
        formatstr : :obj:`str`, optional
            Format string for cell data, defaults to ``"{}"``.
        tooltip : :obj:`(str, str) -> str`, optional
            Function which takes vertical class and horizontal class (strings) as arguments and returns
            desired tooltip as a string. Defaults to no tooltips.
        """
        def _isinvalid(x):
            return isnan(x) or isinf(x)

        if colsum is None:
            colsum = matrix.sum(axis=0)
        if rowsum is None:
            rowsum = matrix.sum(axis=1)

        selected_indexes = {(ind.row() - 2, ind.column() - 2) for ind in self.selectedIndexes()}

        for i in range(len(self.classesv)):
            for j in range(len(self.classesh)):
                val = matrix[i, j]
                col_val = float('nan') if colors is None else colors[i, j]
                item = QStandardItem()
                item.setData(
                    "NA" if _isinvalid(val) else formatstr.format(val),
                    Qt.DisplayRole)
                bkcolor = QColor.fromHsl(
                    [0, 240][i == j], 160,
                    255 if _isinvalid(col_val) else int(255 - 30 * col_val))
                item.setData(QBrush(bkcolor), Qt.BackgroundRole)
                item.setData("trbl", BorderRole)
                if tooltip is not None:
                    item.setToolTip(tooltip(self.classesv[i], self.classesh[j]))
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self._set_item(i + 2, j + 2, item)

        bold_font = self.tablemodel.invisibleRootItem().font()
        bold_font.setBold(True)

        def _sum_item(value, border=""):
            item = QStandardItem()
            item.setData(value, Qt.DisplayRole)
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            item.setFlags(Qt.ItemIsEnabled)
            item.setFont(bold_font)
            item.setData(border, BorderRole)
            item.setData(QColor(192, 192, 192), BorderColorRole)
            return item

        for i in range(len(self.classesh)):
            self._set_item(len(self.classesv) + 2, i + 2, _sum_item(int(colsum[i]), "t"))
        for i in range(len(self.classesv)):
            self._set_item(i + 2, len(self.classesh) + 2, _sum_item(int(rowsum[i]), "l"))
        self._set_item(len(self.classesv) + 2, len(self.classesh) + 2, _sum_item(int(rowsum.sum())))

        self.set_selection(selected_indexes)

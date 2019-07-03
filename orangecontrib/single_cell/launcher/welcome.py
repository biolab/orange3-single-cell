"""
Welcome Screen Dialog

Copied from https://github.com/biolab/orange3/pull/2755
"""
from types import SimpleNamespace
from typing import List  # pylint: disable=unused-import

import pkg_resources
from AnyQt.QtWidgets import (
    QDialog, QWidget, QToolButton, QCheckBox, QAction,
    QHBoxLayout, QVBoxLayout, QSizePolicy, QLabel,
    QListView, QDialogButtonBox, QStackedWidget,
    QStyle, QStyledItemDelegate, QStyleOption, QStyleOptionViewItem,
    QFrame,
    QPushButton, QShortcut)
from AnyQt.QtGui import (
    QFont, QIcon, QPixmap, QImage, QPainter, QColor, QBrush,
    QStandardItemModel, QStandardItem, QDesktopServices,
    QKeySequence)

from AnyQt.QtCore import (
    Qt, QEvent, QRect, QSize, QPoint, QModelIndex, QItemSelectionModel, QUrl,
    QSettings)
from AnyQt.QtCore import pyqtSignal as Signal, pyqtProperty as Property
from orangecanvas.application import examples
from orangecanvas.application.canvasmain import canvas_icons, CanvasMainWindow

from orangecanvas.gui.dropshadow import DropShadowFrame
from orangecanvas.preview import previewbrowser, previewmodel
from orangecanvas import config

from orangecontrib.single_cell.launcher.config import welcome_screen_specs
from orangecontrib.single_cell.launcher.iconview import LinearIconView


def decorate_welcome_icon(icon, background_color):
    """Return a `QIcon` with a circle shaped background.
    """
    welcome_icon = QIcon()
    sizes = [32, 48, 64, 80, 128, 256]
    background_color = QColor(background_color)
    for size in sizes:
        icon_size = QSize(5 * size / 8, 5 * size / 8)
        icon_rect = QRect(QPoint(0, 0), icon_size)
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.transparent)
        p = QPainter(pixmap)
        p.setRenderHint(QPainter.Antialiasing, True)
        p.setBrush(QBrush(background_color))
        p.setPen(Qt.NoPen)
        ellipse_rect = QRect(0, 0, size, size)
        p.drawEllipse(ellipse_rect)
        icon_rect.moveCenter(ellipse_rect.center())
        icon.paint(p, icon_rect, Qt.AlignCenter, )
        p.end()

        welcome_icon.addPixmap(pixmap)

    return welcome_icon


class PagedWidget(QFrame):
    class Page(SimpleNamespace):
        icon = ...  # type: QIcon
        text = ...  # type: str
        toolTip = ...  # type: str
        widget = ...   # type: QWidget

    class TabView(LinearIconView):
        def __init__(self, *args, focusPolicy=Qt.TabFocus, **kwargs):
            super().__init__(*args, focusPolicy=focusPolicy, **kwargs)

        def viewOptions(self):
            # type: () -> QStyleOptionViewItem
            option = super().viewOptions()
            # by default items in views are active only if the view is focused
            if self.isActiveWindow():
                option.state |= QStyle.State_Active
            return option

        def selectionCommand(self, index, event=None):
            # type: (QModelIndex, QEvent) -> QItemSelectionModel.SelectionFlags
            command = super().selectionCommand(index, event)
            if not index.isValid():
                # Prevent deselection on click/drag in an empty view part
                return QItemSelectionModel.NoUpdate
            else:
                # Prevent deselect on click + ctrl modifier
                return command & ~QItemSelectionModel.Deselect

    class TabViewDelegate(QStyledItemDelegate):
        def sizeHint(self, option, index):
            # type: (QStyleOptionViewItem, QModelIndex) -> QSize
            sh = super().sizeHint(option, index)
            widget = option.widget
            if isinstance(widget, PagedWidget.TabView):
                if widget.flow() == QListView.TopToBottom:
                    return sh.expandedTo(QSize(82, 100))
                else:
                    return sh.expandedTo(QSize(100, 82))
            else:
                return sh

        def initStyleOption(self, option, index):
            # type: (QStyleOptionViewItem, QModelIndex) -> None
            super().initStyleOption(option, index)
            widget = option.widget
            if isinstance(widget, PagedWidget.TabView):
                # extend the item rect to cover the whole viewport
                # (probably not a good idea).
                if widget.flow() == QListView.TopToBottom:
                    option.rect.setLeft(0)
                    option.rect.setRight(widget.viewport().width())
                else:
                    option.rect.setTop(0)
                    option.rect.setBottom(widget.viewport().height())

            if option.state & QStyle.State_Selected:
                # make sure the selection highlights cover the whole area
                option.showDecorationSelected = True

    #: Signal emitted when the current displayed widget changes
    currentIndexChanged = Signal(int)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__pages = []  # type: List[PagedWidget.Page]
        self.__currentIndex = -1

        self.setContentsMargins(0, 0, 0, 0)
        self.setLayout(QHBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        self.layout().setSpacing(0)

        self.__tabview = PagedWidget.TabView(
            viewMode=QListView.IconMode,
            flow=QListView.TopToBottom,
            editTriggers=QListView.NoEditTriggers,
            uniformItemSizes=True,
            horizontalScrollBarPolicy=Qt.ScrollBarAlwaysOff
        )
        self.__tabview.setAttribute(Qt.WA_LayoutUsesWidgetRect)
        self.__tabview.setContentsMargins(0, 0, 0, 0)
        self.__tabview.setSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Expanding)

        self.__tabview.setItemDelegate(PagedWidget.TabViewDelegate())
        self.__tabview.setModel(QStandardItemModel(self))
        self.__tabview.selectionModel().selectionChanged.connect(
            self.__on_activated, Qt.UniqueConnection
        )
        iconsize = self.style().pixelMetric(QStyle.PM_LargeIconSize) * 3 // 2
        self.__tabview.setIconSize(QSize(iconsize, iconsize))
        self.__tabview.setAttribute(Qt.WA_MacShowFocusRect, False)

        self.__stack = QStackedWidget(objectName="contents")
        self.__stack.setContentsMargins(0, 0, 0, 0)
        self.__stack.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.layout().addWidget(self.__tabview)
        self.layout().addWidget(self.__stack)

    def currentIndex(self):
        # type: () -> int
        return self.__currentIndex

    def setCurrentIndex(self, index):
        # type: (int) -> None
        assert index < self.count()
        if self.__currentIndex != index:
            self.__currentIndex = index
            if index < 0:
                self.__tabview.selectionModel().clearSelection()
            else:
                self.__tabview.selectionModel().select(
                    self.__tabview.model().index(index, 0),
                    QItemSelectionModel.ClearAndSelect
                )
            self.__stack.setCurrentIndex(index)
            self.currentIndexChanged.emit(index)

    def count(self):
        # type: () -> int
        return len(self.__pages)

    def addPage(self, icon, text, widget):
        # type: (QIcon, str, QWidget) -> int
        return self.insertPage(len(self.__pages), icon, text, widget)

    def insertPage(self, index, icon, text, widget):
        # type: (int, QIcon, str, QWidget) -> int
        if not 0 <= index < self.count():
            index = self.count()

        page = PagedWidget.Page(
            icon=QIcon(icon), text=text, toolTip="", widget=widget
        )
        item = QStandardItem()
        item.setIcon(icon)
        item.setText(text)

        self.__pages.insert(index, page)
        self.__tabview.model().insertRow(index, item)
        self.__stack.insertWidget(index, page.widget)

        if len(self.__pages) == 1:
            self.setCurrentIndex(0)
        elif index <= self.__currentIndex:
            self.__currentIndex += 1
        return index

    def removePage(self, index):
        # type: (int) -> None
        if 0 <= index < len(self.__pages):
            page = self.__pages[index]
            model = self.__tabview.model()  # type: QStandardItemModel
            currentIndex = self.__currentIndex
            if index < currentIndex:
                newCurrent = currentIndex - 1
            else:
                newCurrent = currentIndex
            selmodel = self.__tabview.selectionModel()
            selmodel.selectionChanged.disconnect(self.__on_activated)
            model.removeRow(index)
            del self.__pages[index]
            self.__stack.removeWidget(page.widget)
            selmodel.selectionChanged.connect(
                self.__on_activated, Qt.UniqueConnection)
            self.setCurrentIndex(newCurrent)

    def widget(self, index):
        # type: (int) -> QWidget
        return self.__pages[index].widget

    def setPageEnabled(self, index, enabled):
        # type: (int, bool) -> None
        item = self.__tabview.model().item(index)  # type: QStandardItem
        if item is not None:
            flags = item.flags()
            if enabled:
                flags = flags | Qt.ItemIsEnabled | Qt.ItemIsSelectable
            else:
                flags = flags & ~(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            item.setFlags(flags)

    def isPageEnabled(self, index):
        # type: (int) -> bool
        item = self.__tabview.model().item(index)
        return bool(item.flags() & Qt.ItemIsEnabled)

    def setPageToolTip(self, index, toolTip):
        # type: (int, str) -> None
        if 0 <= index < self.count():
            model = self.__tabview.model()  # type: QStandardItemModel
            item = model.item(index, 0)
            item.setToolTip(toolTip)

    def pageToolTip(self, index):
        model = self.__tabview.model()  # type: QStandardItemModel
        return model.item(index, 0).toolTip()

    def __on_activated(self, selected, deselected):
        indexes = selected.indexes()
        if len(indexes) == 1:
            self.setCurrentIndex(indexes[0].row())
        elif len(indexes) == 0:
            self.setCurrentIndex(-1)
        else:
            assert False, "Invalid selection mode"


class PagedDialog(QDialog):
    """
    A paged dialog widget.
    A paged widget dialog displays a tabbed paged interface
    """
    currentIndexChanged = Signal(int)

    class BottomBar(QWidget):
        def paintEvent(self, event):
            style = self.style()  # type: QStyle
            option = QStyleOption()
            option.initFrom(self)
            p = QPainter(self)
            style.drawPrimitive(QStyle.PE_PanelStatusBar, option, p, self)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setContentsMargins(0, 0, 0, 0)
        self.setLayout(QVBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        self.layout().setSpacing(0)
        self.__pageview = PagedWidget()
        self.__pageview.currentIndexChanged.connect(self.currentIndexChanged)

        self.__bottom = PagedDialog.BottomBar(objectName="bottom-area")
        self.__bottom.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.__bottom.setLayout(QHBoxLayout())

        self.__buttons = QDialogButtonBox(objectName="dialog-buttons")
        self.__buttons.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        self.__buttons.setVisible(False)
        self.__buttons.rejected.connect(self.reject)

        self.__bottom.layout().addWidget(self.__buttons)

        self.layout().addWidget(self.__pageview)
        self.layout().addWidget(self.__bottom)

    def currentIndex(self):
        # type: () -> int
        return self.__pageview.currentIndex()

    def setCurrentIndex(self, index):
        # type: (int) -> None
        self.__pageview.setCurrentIndex(index)

    def count(self):
        # type: () -> int
        return self.__pageview.count()

    def addPage(self, icon, text, widget):
        # type: (QIcon, str, QWidget) -> int
        return self.__pageview.addPage(icon, text, widget)

    def insertPage(self, index, icon, text, widget):
        # type: (int, QIcon, str, QWidget) -> int
        return self.__pageview.insertPage(index, icon, text, widget)

    def removePage(self, index):
        # type: (int) -> None
        return self.__pageview.removePage(index)

    def widget(self, index):
        # type: (int) -> QWidget
        return self.__pageview.widget(index)

    def setPageEnabled(self, index, enabled):
        # type: (int, bool) -> None
        self.__pageview.setPageEnabled(index, enabled)

    def isPageEnabled(self, index):
        # type: (int) -> bool
        return self.__pageview.isPageEnabled(index)

    def buttonBox(self):
        # type: () -> QDialogButtonBox
        """
        Return a QDialogButtonBox instance.
        """
        return self.__buttons


def pixmap_from_image(image):
    # type: (QImage) -> QPixmap
    pixmap = QPixmap.fromImage(image)  # type: QPixmap
    if hasattr(pixmap, "setDevicePixelRatio"):
        pixmap.setDevicePixelRatio(image.logicalDpiX() / 72)
    else:
        pixmap = pixmap.scaled(
            (image.size() * 72) / image.logicalDpiX(),
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )
    return pixmap


class FancyWelcomeScreen(QWidget):
    """
    Fancy welcome screen.
    +-----------+
    |  Welcome  |
    +-----------+
    | A | B | C |
    +---+---+---+
    The upper part consist of static image while the lower items select some
    prespecified action.
    """
    class StartItem(QWidget):
        """
        An active item in the bottom row of the welcome screen.
        """
        def __init__(self, *args, text="", icon=QIcon(), iconSize=QSize(),
                     iconActive=QIcon(), **kwargs):
            self.__iconSize = QSize()
            self.__icon = QIcon()
            self.__icon_active = QIcon()
            self.__text = ""
            self.__active = False
            super().__init__(*args, **kwargs)
            self.setAutoFillBackground(True)
            font = self.font()
            font.setPointSize(18)
            self.setFont(font)
            self.setAttribute(Qt.WA_SetFont, False)
            self.setText(text)
            self.setIcon(icon)
            self.setIconSize(iconSize)
            self.setIconActive(iconActive)
            self.installEventFilter(self)

        def iconSize(self):
            if not self.__iconSize.isValid():
                size = self.style().pixelMetric(
                    QStyle.PM_LargeIconSize, None, self) * 2
                return QSize(size, size)
            else:
                return QSize(self.__iconSize)

        def setIconSize(self, size):
            if size != self.__iconSize:
                self.__iconSize = QSize(size)
                self.updateGeometry()

        iconSize_ = Property(QSize, iconSize, setIconSize, designable=True)

        def icon(self):
            if self.__active:
                return QIcon(self.__icon_active)
            else:
                return QIcon(self.__icon)

        def setIcon(self, icon):
            self.__icon = QIcon(icon)
            self.update()

        icon_ = Property(QIcon, icon, setIcon, designable=True)

        def iconActive(self):
            return QIcon(self.__icon_active)

        def setIconActive(self, icon):
            self.__icon_active = QIcon(icon)
            self.update()

        icon_active_ = Property(QIcon, iconActive, setIconActive, designable=True)

        def sizeHint(self):
            return QSize(200, 150)

        def setText(self, text):
            if self.__text != text:
                self.__text = text
                self.updateGeometry()
                self.update()

        def text(self):
            return self.__text

        text_ = Property(str, text, setText, designable=True)

        def initStyleOption(self, option):
            # type: (QStyleOptionViewItem) -> None
            option.initFrom(self)
            option.backgroundBrush = option.palette.brush(self.backgroundRole())
            option.font = self.font()
            option.text = self.text()
            option.icon = self.icon()

            option.decorationPosition = QStyleOptionViewItem.Top
            option.decorationAlignment = Qt.AlignCenter
            option.decorationSize = self.iconSize()
            option.displayAlignment = Qt.AlignCenter
            option.features = (
                QStyleOptionViewItem.WrapText |
                QStyleOptionViewItem.HasDecoration |
                QStyleOptionViewItem.HasDisplay
            )
            option.showDecorationSelected = True
            option.widget = self

        def paintEvent(self, event):
            style = self.style()  # type: QStyle
            painter = QPainter(self)
            option = QStyleOption()
            option.initFrom(self)
            style.drawPrimitive(QStyle.PE_Widget, option, painter, self)

            option = QStyleOptionViewItem()
            self.initStyleOption(option)
            style.drawControl(QStyle.CE_ItemViewItem, option, painter, self)

        def eventFilter(self, obj, event):
            try:
                if event.type() == QEvent.Enter:
                    self.__active = True
                    self.setCursor(Qt.PointingHandCursor)
                    self.update()
                    return True
                elif event.type() == QEvent.Leave:
                    self.__active = False
                    self.unsetCursor()
                    self.update()
                    return True
            except Exception as ex:
                pass
            return False

    #: Signal emitted when the current selected item in changes.
    currentChanged = Signal(int)

    #: Signal emitted when the item is double clicked.
    activated = Signal(int)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFocusPolicy(Qt.TabFocus)
        vlayout = QVBoxLayout(spacing=0)
        vlayout.setContentsMargins(0, 0, 0, 0)
        self.__currentIndex = -1

        self.__contents = QLabel(
            sizePolicy=QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        )
        vlayout.addWidget(self.__contents)

        hlayout = QHBoxLayout(spacing=0)
        hlayout.setContentsMargins(0, 0, 0, 0)

        self.__items = items = [
            FancyWelcomeScreen.StartItem(objectName="item"),
            FancyWelcomeScreen.StartItem(objectName="item"),
            FancyWelcomeScreen.StartItem(objectName="item"),
        ]
        items[0].setProperty("-position", QStyleOptionViewItem.Beginning)
        items[1].setProperty("-position", QStyleOptionViewItem.Middle)
        items[-1].setProperty("-position", QStyleOptionViewItem.End)

        for item in items:
            hlayout.addWidget(item)

        vlayout.addLayout(hlayout)
        self.setLayout(vlayout)
        self.setCurrentIndex(0)

    def setImage(self, image):
        # type: (QImage) -> None
        """
        Set the welcome image.
        Parameters
        ----------
        image : QImage
        """
        pixmap = pixmap_from_image(image)
        self.__contents.setPixmap(pixmap)

    def setItem(self, index, image):
        item = self.layout().itemAt(1).layout().itemAt(index)
        widget = item.widget()  # type: FancyWelcomeScreen.StartItem
        widget.setIcon(image)

    def item(self, index):
        item = self.layout().itemAt(1).layout().itemAt(index)
        return item.widget()

    def setItemText(self, index, text):
        item = self.layout().itemAt(1).layout().itemAt(index)
        widget = item.widget()  # type: FancyWelcomeScreen.StartItem
        widget.setText(text)

    def setItemIcon(self, index, icon):
        item = self.item(index)
        item.setIcon(icon)

    def setItemActiveIcon(self, index, icon):
        item = self.item(index)  # type: FancyWelcomeScreen.StartItem
        item.setIconActive(icon)

    def setItemToolTip(self, index, tip):
        item = self.item(index)
        item.setToolTip(tip)

    def setCurrentIndex(self, index):
        self.__currentIndex = index

    def currentIndex(self):
        return self.__currentIndex

    def __indexAtPos(self, pos):
        # type: (QPoint) -> int
        for i, item in enumerate(self.__items):
            if item.geometry().contains(pos):
                return i
        return -1

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            index = self.__indexAtPos(event.pos())
            if index != -1:
                self.setCurrentIndex(index)
            event.accept()
        else:
            event.ignore()

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.LeftButton:
            index = self.__indexAtPos(event.pos())
            if index != -1:
                self.setCurrentIndex(index)
            event.accept()
        else:
            event.ignore()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            index = self.__indexAtPos(event.pos())
            if index != -1:
                self.activated.emit(index)
            event.accept()
        else:
            event.ignore()

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.LeftButton:
            index = self.__indexAtPos(event.pos())
            if index != -1:
                self.activated.emit(index)
            event.accept()
        else:
            event.ignore()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Right:
            direction = 1
        elif event.key() == Qt.Key_Left:
            direction = -1
        else:
            super().keyPressEvent(event)
            return

        event.accept()
        if len(self.__items):
            index = self.__currentIndex + direction
            self.setCurrentIndex(max(0, min(index, len(self.__items) - 1)))


class SingleLinkPage(QFrame):
    """
    An simple (overly) large image with a external link
    """
    def __init__(self, *args, image=QImage(), heading="", link=QUrl(), **kwargs):
        super().__init__(*args, **kwargs)
        self.__link = QUrl()

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        self.setLayout(layout)
        self.__heading = QLabel()
        self.__content = QLabel()

        self.layout().addWidget(self.__heading)
        self.layout().addWidget(self.__content, 10, Qt.AlignCenter)

        self.__shadow = DropShadowFrame()
        self.__shadow.setWidget(self.__content)

        self.setImage(image)
        self.setHeading(heading)
        self.setLink(link)

    def setHeading(self, heading):
        self.__heading.setText("<h2>{}</h2>".format(heading))

    def setImage(self, image):
        pm = pixmap_from_image(image)
        self.__content.setPixmap(pm)

    def setLink(self, url):
        self.__link = QUrl(url)
        if not self.__link.isEmpty():
            self.__content.setCursor(Qt.PointingHandCursor)
        else:
            self.__content.unsetCursor()
        self.__content.setToolTip(self.__link.toString())

    def mousePressEvent(self, event):
        if self.__content.geometry().contains(event.pos()) and \
                not self.__link.isEmpty():
            QDesktopServices.openUrl(self.__link)
            event.accept()
        else:
            super().mousePressEvent(event)


def resource_path(name):
    return pkg_resources.resource_filename(
        "orangecontrib.single_cell.launcher", name)


def sc_icon(filename):
    return QIcon(pkg_resources.resource_filename(
        "orangecontrib.single_cell.launcher", "icons/"+filename))


def welcome_dialog_paged(self):
    # type: (CanvasMainWindow) -> None
    """
    Show a modal multipaged welcome screen.
    """
    dlg = PagedDialog(
        self, windowTitle=self.tr("Orange Data Mining"),
    )
    dlg.setWindowModality(Qt.ApplicationModal)
    dlg.setAttribute(Qt.WA_DeleteOnClose)
    dlg.layout().setSizeConstraint(QVBoxLayout.SetFixedSize)
    dlg.setStyleSheet("""
        TabView::item:selected {
               background: rgb(243, 171, 86);
        }
    """)
    main = FancyWelcomeScreen()
    spec = welcome_screen_specs()
    if spec.image:
        background = QImage(spec.image)
    else:
        background = QImage("canvas_icons:orange-start-background.png")
    main.setImage(background)

    if spec.css:
        main.setStyleSheet(spec.css)
    else:
        main.setStyleSheet(
            "StartItem { background-color: rgb(123, 164, 214) }"
        )

    def decorate_icon(icon):
        return decorate_welcome_icon(icon, "#6dacb2")

    for i, item in zip(range(3), spec.items):
        main.setItemText(i, item.text)
        main.setItemToolTip(i, item.tooltip)
        main.setItemIcon(i, decorate_icon(QIcon(item.icon)))
        main.setItemActiveIcon(i, decorate_icon(QIcon(item.active_icon)))
        main.item(i).setProperty("path", item.path)

    main.setCurrentIndex(0)
    main.activated.connect(lambda: openselectedbutton.click())

    PageWelcome = dlg.addPage(
        sc_icon("Welcome.svg"), "Welcome", main
    )
    examples_ = examples.workflows(config.default)
    items = [previewmodel.PreviewItem(path=t.abspath()) for t in examples_]
    model = previewmodel.PreviewModel(items=items)
    model.delayedScanUpdate()
    browser = previewbrowser.PreviewBrowser()
    browser.setModel(model)

    PageTemplates = dlg.addPage(
        sc_icon("Templates.svg"), "Templates", browser
    )
    browser.activated.connect(lambda: openselectedbutton.click())

    recent = [previewmodel.PreviewItem(name=item.title, path=item.path)
              for item in self.recent_schemes]
    model = previewmodel.PreviewModel(items=recent)
    browser = previewbrowser.PreviewBrowser()
    browser.setModel(model)
    model.delayedScanUpdate()

    PageRecent = dlg.addPage(
        self.recent_action.icon(), "Recent", browser
    )
    browser.activated.connect(lambda: openselectedbutton.click())
    dlg.setPageEnabled(PageRecent, model.rowCount() > 0)

    page = SingleLinkPage(
        image=QImage(resource_path("icons/getting-started-video-tutorials.png")),
        heading="Getting Started",
        link=QUrl("https://www.youtube.com/watch?v=3nMcI4Hxm7c"),
    )
    page.setContentsMargins(25, 25, 25, 25)
    PageGetStarted = dlg.addPage(
        canvas_icons("YouTube.svg"), "Get Started", page,
    )
    buttons = dlg.buttonBox()
    buttons.setVisible(True)
    buttons.setStandardButtons(QDialogButtonBox.Open |
                               QDialogButtonBox.Cancel)
    # choose the selected workflow button
    openselectedbutton = buttons.button(QDialogButtonBox.Open)
    openselectedbutton.setText(self.tr("Open"))
    openselectedbutton.setToolTip("Open the selected workflow")
    openselectedbutton.setDefault(True)

    newbutton = QPushButton(
        "New", toolTip="Create a new workflow")
    s = QShortcut(QKeySequence.New, newbutton)
    s.activated.connect(newbutton.click)
    buttons.addButton(newbutton, QDialogButtonBox.AcceptRole)

    openexisting = QPushButton(
        "Open Existing\N{HORIZONTAL ELLIPSIS}",
        toolTip="Open an existing workflow file"
    )
    s = QShortcut(QKeySequence.Open, dlg)
    s.activated.connect(openexisting.click)
    buttons.addButton(openexisting, QDialogButtonBox.AcceptRole)

    settings = QSettings()

    show_start_key = "startup/show-welcome-screen"
    show_start = QCheckBox(
        "Show at startup",
        checked=settings.value(show_start_key, True, type=bool)
    )
    # Abusing ResetRole to push the check box to the left in all button
    # layouts.
    buttons.addButton(show_start, QDialogButtonBox.ResetRole)

    def update_show_at_startup(value):
        settings.setValue(show_start_key, value)

    show_start.toggled.connect(update_show_at_startup)

    def on_page_changed(index):
        if index == PageWelcome:
            openselectedbutton.setEnabled(True)
        elif index == PageTemplates:
            openselectedbutton.setEnabled(bool(examples))
        elif index == PageRecent:
            openselectedbutton.setEnabled(bool(recent))
        elif index == PageGetStarted:
            openselectedbutton.setEnabled(False)
        else:
            openselectedbutton.setEnabled(False)

    dlg.currentIndexChanged.connect(on_page_changed)

    def open_example_workflow(path):
        # open a workflow without filename/directory tracking.
        wf = self.new_scheme_from(path)
        if self.is_transient():
            window = self
        else:
            window = self.create_new_window()
            window.show()
            window.raise_()
            window.activateWindow()
        window.set_new_scheme(wf)

    def open_url(url):
        return QDesktopServices.openUrl(QUrl(url))

    def on_clicked(button):
        current = dlg.currentIndex()
        path = None
        open_workflow_file = None

        if current == PageWelcome:
            open_workflow_file = open_url
        elif current == PageTemplates:
            open_workflow_file = open_example_workflow
        elif current == PageRecent:
            open_workflow_file = self.open_scheme_file

        if button is openselectedbutton and \
                        current in {PageTemplates, PageRecent}:
            w = dlg.widget(current)
            assert isinstance(w, previewbrowser.PreviewBrowser)
            assert w.currentIndex() != -1
            model = w.model()
            item = model.item(w.currentIndex())
            path = item.path()
        elif button is openselectedbutton and current == PageWelcome:
            w = dlg.widget(current)
            assert isinstance(w, FancyWelcomeScreen)
            assert w.currentIndex() != -1
            path = w.item(w.currentIndex()).property("path")

        if path is not None:
            open_workflow_file(path)
            dlg.accept()

    buttons.clicked.connect(on_clicked)

    def on_open_existing():
        filedlg = self._open_workflow_dialog()
        filedlg.fileSelected.connect(self.open_scheme_file)
        filedlg.accepted.connect(dlg.accept)
        filedlg.exec()

    openexisting.clicked.connect(on_open_existing)

    def new_window():
        if not self.is_transient():
            self.new_workflow_window()
        dlg.accept()

    newbutton.clicked.connect(new_window)

    dlg.show()

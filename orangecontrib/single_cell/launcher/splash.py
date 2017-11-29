import pkg_resources
from PyQt5.QtCore import QCoreApplication, QPoint, Qt, QRect
from PyQt5.QtGui import QPixmap, QFont, QFontMetrics, QPainter, QColor


def splash_screen():
    """
    """
    pm = QPixmap(
        pkg_resources.resource_filename(
            __name__, "icons/splash.png")
    )

    version = QCoreApplication.applicationVersion()
    size = 21 if len(version) < 5 else 16
    font = QFont("Helvetica")
    font.setPixelSize(size)
    font.setLetterSpacing(QFont.AbsoluteSpacing, 2)
    metrics = QFontMetrics(font)
    br = metrics.boundingRect(version).adjusted(-5, 0, 5, 0)
    br.moveCenter(QPoint(388, 215))

    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setRenderHint(QPainter.TextAntialiasing)
    p.setFont(font)
    p.setPen(QColor("#ffffff"))
    p.drawText(br, Qt.AlignCenter, version)
    p.end()
    return pm, QRect(36, 327, 250, 348)
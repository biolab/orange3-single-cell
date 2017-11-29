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
    size = 21 if len(version) < 6 else 16
    font = QFont("Arial")
    font.setPixelSize(size)
    metrics = QFontMetrics(font)
    br = metrics.boundingRect(version).adjusted(-5, 0, 5, 0)
    br.moveCenter(QPoint(367, 228))

    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setRenderHint(QPainter.TextAntialiasing)
    p.setFont(font)
    p.setPen(QColor("#4c85c5"))
    p.drawText(br, Qt.AlignCenter, version)
    p.end()
    return pm, QRect(33, 314, 250, 328)

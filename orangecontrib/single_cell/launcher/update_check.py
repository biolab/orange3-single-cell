import logging
import time
import sys
import pkg_resources

from AnyQt.QtCore import QSettings, QThread, pyqtSignal, Qt, QUrl
from AnyQt.QtGui import QDesktopServices
from AnyQt.QtWidgets import QMessageBox
from urllib.request import urlopen, Request
from Orange.version import version as current


log = logging.getLogger(__name__)
VERSION_URL = 'https://singlecell.biolab.si/version/'
DOWNLOAD_URL = 'https://singlecell.biolab.si/download/'
USER_AGENT = 'scOrange{orange_version}:Python{py_version}:{platform}:{conda}'
UPDATE_AVAILABLE_TITLE = 'scOrange Update Available'
UPDATE_AVAILABLE_MESSAGE = (
    'A newer version of scOrange is available.<br><br>'
    '<b>Current version:</b> {current_version}<br>'
    '<b>Latest version:</b> {latest_version}'
)


def check_for_updates():
    settings = QSettings()
    check_updates = settings.value('startup/check-updates', True, type=bool)
    last_check_time = settings.value('startup/last-update-check-time', 0, type=int)
    ONE_DAY = 86400

    if check_updates and time.time() - last_check_time > ONE_DAY:
        settings.setValue('startup/last-update-check-time', int(time.time()))

        thread = GetLatestVersion()
        thread.resultReady.connect(compare_versions)
        thread.start()
        return thread


class GetLatestVersion(QThread):
    resultReady = pyqtSignal(str)

    def run(self):
        try:
            request = Request(VERSION_URL,
                              headers={
                                  'Accept': 'text/plain',
                                  'Accept-Encoding': 'gzip, deflate',
                                  'Connection': 'close',
                                  'User-Agent': self.ua_string()})
            contents = urlopen(request, timeout=10).read().decode()
        # Nothing that this fails with should make Orange crash
        except Exception:  # pylint: disable=broad-except
            log.exception('Failed to check for updates')
        else:
            self.resultReady.emit(contents)

    @staticmethod
    def ua_string():
        is_anaconda = 'Continuum' in sys.version or 'conda' in sys.version
        return USER_AGENT.format(
            orange_version=current,
            py_version='.'.join(sys.version[:3]),
            platform=sys.platform,
            conda='Anaconda' if is_anaconda else '',
        )


def current_version():
    return pkg_resources.get_distribution("Orange3-SingleCell").version


def compare_versions(latest):
    current = current_version()
    version = pkg_resources.parse_version
    if version(latest) <= version(current):
        return
    question = QMessageBox(
        QMessageBox.Information,
        UPDATE_AVAILABLE_TITLE,
        UPDATE_AVAILABLE_MESSAGE.format(
            current_version=current,
            latest_version=latest),
        textFormat=Qt.RichText)
    ok = question.addButton('Download Now', question.AcceptRole)
    question.setDefaultButton(ok)
    question.addButton('Remind Later', question.RejectRole)
    question.finished.connect(
        lambda:
        question.clickedButton() == ok and
        QDesktopServices.openUrl(QUrl(DOWNLOAD_URL)))
    question.show()

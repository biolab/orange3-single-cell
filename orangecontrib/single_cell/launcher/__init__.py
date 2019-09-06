from typing import Iterable, List

import pkg_resources

from AnyQt.QtCore import Qt
from AnyQt.QtGui import QColor

from orangecanvas.gui.splashscreen import SplashScreen

from Orange.canvas import config as oconfig, __main__ as main, mainwindow

from orangecontrib.single_cell.launcher.splash import splash_screen
from orangecontrib.single_cell.launcher.update_check import check_for_updates
from orangecontrib.single_cell.launcher.welcome import welcome_dialog_paged


class SCConfig(oconfig.Config):
    ApplicationName = "scOrange"
    dist = pkg_resources.get_distribution("Orange3-SingleCell")
    ApplicationVersion = dist.version
    del dist

    @staticmethod
    def examples_entry_points():
        # type: () -> Iterable[pkg_resources.EntryPoint]
        return filter(
            lambda ep: ep.dist.project_name.lower() == "orange3-singlecell",
            super(SCConfig, SCConfig).examples_entry_points()
        )

    @staticmethod
    def splash_screen():
        return splash_screen()

    @staticmethod
    def core_packages():
        # type: () -> List[str]
        return super(SCConfig, SCConfig).core_packages() + [
            "orange3-singlecell",
            "orange3-bioinformatics",
        ]


class SCMainWindow(mainwindow.MainWindow):
    def welcome_dialog(self):
        return welcome_dialog_paged(self, )


class SCSplashScreen(SplashScreen):
    def showMessage(self, message, alignment=Qt.AlignLeft, color=Qt.black):
        super().showMessage(message, alignment, color=QColor("#4c85c5"))


class SCOrangeLauncher:
    def launch(self):
        oconfig.Config = SCConfig
        main.SplashScreen = SCSplashScreen
        main.MainWindow = SCMainWindow

        self.fix_application_dirs()
        self.replace_update_check()

        self.main()

    def fix_application_dirs(self):
        import os, sys
        from Orange.misc import environ

        def data_dir(versioned=True):
            """
            Return the platform dependent Orange data directory.

            This is ``data_dir_base()``/scOrange/ directory.
            """
            base = environ.data_dir_base()
            return os.path.join(base, "scOrange")
        environ.data_dir = data_dir

        def cache_dir(*args):
            """
            Return the platform dependent Orange cache directory.
            """
            if sys.platform == "darwin":
                base = os.path.expanduser("~/Library/Caches")
            elif sys.platform == "win32":
                base = os.getenv("APPDATA", os.path.expanduser("~/AppData/Local"))
            elif os.name == "posix":
                base = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
            else:
                base = os.path.expanduser("~/.cache")

            base = os.path.join(base, "scOrange")
            if sys.platform == "win32":
                # On Windows cache and data dir are the same.
                # Microsoft suggest using a Cache subdirectory
                return os.path.join(base, "Cache")
            else:
                return base
        environ.cache_dir = cache_dir

    def replace_update_check(self):
        main.check_for_updates = check_for_updates

    def replace_splash_screen(self):
        main.SplashScreen = SCSplashScreen

    def main(self):
        main.main()

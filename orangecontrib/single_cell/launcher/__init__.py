from Orange.canvas.application.workflows import list_schemes, ExampleWorkflow

from orangecontrib.single_cell.launcher.splash import splash_screen
from orangecontrib.single_cell.launcher.update_check import check_for_updates


class SCOrangeLauncher:
    def launch(self):
        self.fix_application_name()
        self.fix_application_dirs()
        self.replace_splash_screen()
        self.replace_welcome_screen()
        self.replace_example_workflows()
        self.replace_update_check()

        self.fix_widget_categories()

        self.main()

    def fix_application_name(self):
        from PyQt5.QtCore import QCoreApplication, QSettings
        from Orange.canvas import config

        def init():
            """
            Initialize the QCoreApplication.organizationDomain, applicationName,
            applicationVersion and the default settings format. Will only run once.

            .. note:: This should not be run before QApplication has been initialized.
                      Otherwise it can break Qt's plugin search paths.

            """
            import pkg_resources

            dist = pkg_resources.get_distribution("Orange3-SingleCell")
            version = dist.version
            # Use only major.minor
            #version = ".".join(version.split(".", 2)[:2])

            QCoreApplication.setOrganizationDomain("biolab.si")
            QCoreApplication.setApplicationName("scOrange")
            QCoreApplication.setApplicationVersion(version)
            QSettings.setDefaultFormat(QSettings.IniFormat)

            # Make it a null op.
            config.init = lambda: None
        config.init = init

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

    def replace_welcome_screen(self):
        from Orange.canvas.application.canvasmain import CanvasMainWindow
        from orangecontrib.single_cell.launcher.welcome import welcome_dialog_paged

        CanvasMainWindow.welcome_dialog = welcome_dialog_paged

    def replace_example_workflows(self):
        from Orange.canvas.application import workflows

        def example_workflows():
            from orangecontrib.single_cell import tutorials
            workflows = list_schemes(tutorials)
            workflows = [ExampleWorkflow(wf, tutorials, "scOrange")
                         for wf in workflows]
            return workflows
        workflows.example_workflows = example_workflows

    def replace_update_check(self):
        from Orange.canvas import __main__

        __main__.check_for_updates = check_for_updates

    def replace_splash_screen(self):
        from AnyQt.QtCore import Qt
        from AnyQt.QtGui import QColor
        from Orange.canvas import config
        from Orange.canvas.gui.splashscreen import SplashScreen

        config.splash_screen = splash_screen

        sm = SplashScreen.showMessage

        def showMessage(self, message, alignment=Qt.AlignLeft, color=Qt.black):
            sm(self, message, alignment=alignment, color=QColor("#4c85c5"))
        SplashScreen.showMessage = showMessage

    def fix_widget_categories(self):
        from Orange.canvas import registry

        wd = registry.WidgetDiscovery.widget_description
        def widget_description(self, module, widget_name=None,
                               category_name=None, distribution=None):
            """Move scOrange widgets to more appropriate methods

            t-SNE -> Unsupervised
            Louvain Clustering -> Unsupervised

            A better way to do this would be to change the behaviour of the
            widget_description method to *not* overwrite the category when
            manually specified on widget.
            """

            desc = wd(
                self, module, widget_name, category_name, distribution)

            if desc.qualified_name == 'orangecontrib.single_cell.widgets.owtsne.OWtSNE':
                desc.category = 'Unsupervised'
            if desc.qualified_name == 'orangecontrib.single_cell.widgets.owlouvainclustering.OWLouvainClustering':
                desc.category = 'Unsupervised'
            return desc

        registry.WidgetDiscovery.widget_description = widget_description

    def main(self):
        from Orange.canvas.__main__ import main

        main()

#!/usr/bin/env python3

from os import path, walk

import sys
from setuptools import setup, find_packages, Command
import subprocess

NAME = "Orange3-SingleCell"

VERSION = "0.8.2"

DESCRIPTION = "Add-on for bioinformatics analysis of single cell data"
LONG_DESCRIPTION = open(path.join(path.dirname(__file__), 'README.md')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'info@biolab.si'
URL = "https://github.com/biolab/orange3-single-cell"
DOWNLOAD_URL = "https://github.com/biolab/orange3-single-cell/tarball/{}".format(VERSION)
LICENSE = 'GPLv3+'

KEYWORDS = (
    # [PyPi](https://pypi.python.org) packages with keyword "orange3 add-on"
    # can be installed using the Orange Add-on Manager
    'orange3 add-on',
)

PACKAGES = find_packages()

PACKAGE_DATA = {
    'orangecontrib.single_cell.tutorials': ['*.ows', '*.tsv', '*.mtx', '*.tab.gz', '*.tab'],
    'orangecontrib.single_cell.tests': ['*'],
    'orangecontrib.single_cell.launcher.icons': ['*.ows', '*.png'],
    'orangecontrib.single_cell.widgets.icons': ['*.ows']
}

DATA_FILES = [
    # Data files that will be installed outside site-packages folder
]

INSTALL_REQUIRES = sorted(set(
    line.partition('#')[0].strip()
    for line in open(path.join(path.dirname(__file__), 'requirements.txt'))
) - {''})

ENTRY_POINTS = {
    # Entry points that marks this package as an orange add-on. If set, addon will
    # be shown in the add-ons manager even if not published on PyPi.
    'orange3.addon': (
        'single_cell = orangecontrib.single_cell',
    ),
    # Entry point used to specify packages containing tutorials accessible
    # from welcome screen. Tutorials are saved Orange Workflows (.ows files).
    'orange.widgets.tutorials': (
        # Syntax: any_text = path.to.package.containing.tutorials
        'vcftutorials = orangecontrib.single_cell.tutorials',
    ),

    # Entry point used to specify packages containing widgets.
    'orange.widgets': (
        # Syntax: category name = path.to.package.containing.widgets
        # Widget category specification can be seen in
        #    orangecontrib/example/widgets/__init__.py
        'Single Cell = orangecontrib.single_cell.widgets',
    ),

    # Register widget help
    "orange.canvas.help": (
        'html-index = orangecontrib.single_cell.widgets:WIDGET_HELP_PATH',)
}

NAMESPACE_PACKAGES = ["orangecontrib"]

TEST_SUITE = "orangecontrib.single_cell.tests.suite"


def include_documentation(local_dir, install_dir):
    global DATA_FILES
    if (('bdist_wheel' in sys.argv) or ('install' in sys.argv))\
            and not path.exists(local_dir):
        print("Directory '{}' does not exist. "
              "Please build documentation before running bdist_wheel."
              .format(path.abspath(local_dir)))
        sys.exit(0)

    doc_files = []
    for dirpath, dirs, files in walk(local_dir):
        doc_files.append((dirpath.replace(local_dir, install_dir),
                          [path.join(dirpath, f) for f in files]))
    DATA_FILES.extend(doc_files)


class CoverageCommand(Command):
    """A setup.py coverage subcommand developers can run locally."""
    description = "run code coverage"
    user_options = []
    initialize_options = finalize_options = lambda self: None

    def check_requirements(self):

        try:
            import coverage
        except ImportError as e:
            raise e

    def run(self):
        """Check coverage on current workdir"""
        self.check_requirements()

        sys.exit(subprocess.call(r'''
        coverage run setup.py test
        echo; echo
        coverage report
        coverage html &&
            { echo; echo "See also: file://$(pwd)/coverage_html_report/index.html"; echo; }
        ''', shell=True, cwd=path.dirname(path.abspath(__file__))))


if __name__ == '__main__':
    include_documentation('doc/build/htmlhelp', 'help/orange3-single_cell')

    cmdclass = {
        'coverage': CoverageCommand,
    }

    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        packages=PACKAGES,
        package_data=PACKAGE_DATA,
        data_files=DATA_FILES,
        install_requires=INSTALL_REQUIRES,
        entry_points=ENTRY_POINTS,
        keywords=KEYWORDS,
        namespace_packages=NAMESPACE_PACKAGES,
        test_suite=TEST_SUITE,
        include_package_data=True,
        zip_safe=False,
        cmdclass=cmdclass
    )

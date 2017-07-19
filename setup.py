#!/usr/bin/env python3

from os import path, walk

import sys
from setuptools import setup, find_packages

NAME = "Orange3-Variants"

VERSION = "0.0.2"

DESCRIPTION = "Add-on for SNV calling from NGS data"
LONG_DESCRIPTION = open(path.join(path.dirname(__file__), 'README.md')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'info@biolab.si'
URL = "https://github.com/biolab/orange3-variants"
DOWNLOAD_URL = "https://github.com/biolab/orange3-variants/tarball/{}".format(VERSION)
LICENSE = 'GPLv3+'

KEYWORDS = (
    # [PyPi](https://pypi.python.org) packages with keyword "orange3 add-on"
    # can be installed using the Orange Add-on Manager
    'orange3 add-on',
)

PACKAGES = find_packages()

PACKAGE_DATA = {
    'orangecontrib.variants': ['tutorials/*.ows'],
    'orangecontrib.variants.widgets': ['icons/*'],
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
        'variants = orangecontrib.variants',
    ),
    # Entry point used to specify packages containing tutorials accessible
    # from welcome screen. Tutorials are saved Orange Workflows (.ows files).
    'orange.widgets.tutorials': (
        # Syntax: any_text = path.to.package.containing.tutorials
        'vcftutorials = orangecontrib.variants.tutorials',
    ),

    # Entry point used to specify packages containing widgets.
    'orange.widgets': (
        # Syntax: category name = path.to.package.containing.widgets
        # Widget category specification can be seen in
        #    orangecontrib/example/widgets/__init__.py
        'Variants = orangecontrib.variants.widgets',
    ),

    # Register widget help
    "orange.canvas.help": (
        'html-index = orangecontrib.variants.widgets:WIDGET_HELP_PATH',)
}

NAMESPACE_PACKAGES = ["orangecontrib"]

TEST_SUITE = "orangecontrib.variants.tests.suite"


def include_documentation(local_dir, install_dir):
    global DATA_FILES
    if 'bdist_wheel' in sys.argv and not path.exists(local_dir):
        print("Directory '{}' does not exist. "
              "Please build documentation before running bdist_wheel."
              .format(path.abspath(local_dir)))
        sys.exit(0)

    doc_files = []
    for dirpath, dirs, files in walk(local_dir):
        doc_files.append((dirpath.replace(local_dir, install_dir),
                          [path.join(dirpath, f) for f in files]))
    DATA_FILES.extend(doc_files)

if __name__ == '__main__':
    include_documentation('doc/build/html', 'help/orange3-variants')
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
    )

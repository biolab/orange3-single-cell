Orange3 Single Cell
===================

[![Discord Chat](https://img.shields.io/discord/633376992607076354)](https://discord.gg/FWrfeXV)
![main](https://github.com/biolab/orange3-single-cell/workflows/main/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/biolab/orange3-single-cell/branch/master/graph/badge.svg)](https://codecov.io/gh/biolab/orange3-single-cell)
[![PyPI](https://img.shields.io/pypi/v/orange3-singlecell.svg)](https://pypi.org/project/Orange3-SingleCell/)
[![Conda](https://img.shields.io/conda/v/conda-forge/orange3-singlecell.svg)](https://anaconda.org/conda-forge/orange3-singlecell)
[![PyPI - License](https://img.shields.io/pypi/l/orange3-singlecell.svg)](https://pypi.org/project/Orange3-SingleCell/)

The Single Cell add-on for [Orange3](http://orange.biolab.si) adds functionality for analysis of single cell data.

Installation
------------

To install the add-on, run

    pip install .

To register this add-on with Orange, but keep the code in the development directory (do not copy it to 
Python's site-packages directory), run

    pip install -e .

Documentation / widget help can be built by running

    make html htmlhelp

from the doc directory.

Usage
-----

After the installation, the widget from this add-on is registered with Orange. To run Orange from the terminal,
use

    python -m Orange.canvas

The new widget appears in the toolbox bar under the section Single Cell.

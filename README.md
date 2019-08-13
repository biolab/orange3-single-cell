Orange3 Single Cell
======================

<a href="https://dev.azure.com/orange-biolab/orange3-singlecell/">
  <img src="https://dev.azure.com/orange-biolab/orange3-singlecell/_apis/build/status/CI%20Pipeline?branchName=master" />
</a>

<a href="https://codecov.io/gh/biolab/orange3-single-cell">
  <img src="https://codecov.io/gh/biolab/orange3-single-cell/branch/master/graph/badge.svg" />
</a>

<a href="https://pypi.org/project/Orange3-SingleCell/">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/orange3-singlecell.svg" />
</a>

<a href="https://anaconda.org/conda-forge/orange3-singlecell">
  <img alt="Conda" src="https://img.shields.io/conda/v/conda-forge/orange3-singlecell.svg" />
</a>

<a href="https://pypi.org/project/Orange3-SingleCell/">
  <img alt="PyPI - License" src="https://img.shields.io/pypi/l/orange3-singlecell.svg" />
</a>

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

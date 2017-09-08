Orange3 Single Cell
======================

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

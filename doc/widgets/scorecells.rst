Score Cells
===========

.. figure:: icons/ScoreCells.png

Score cells in the single cell expression data according to a set of biomarker genes.

Signals
-------

**Inputs**:

-  **Data**

   A single cell expression data containing cells in rows and genes in columns. Column headers (attribute names) include gene names.

-  **Genes**

   Data table where one of the columns includes names of genes that match those from expression data.

**Outputs**:

-  **Data**

   Single cell gene expression data from the input that has been augmented with a meta column reporting on cell scores.

Description
-----------

**Score Cells** scores the cells (rows) in the input data based on expression of input marker genes. The score is computed independently for each cell and is equal to the maximum expression of the marker genes. If expressions for all the marker genes are missing, the cell obtains a score of zero.

.. figure:: images/ScoreCells-stamped.png

1. A column name from table on genes that provides names of genes that are included in the input gene expression table.
2. Tick to automatically process input data and send the result of scoring to the output. If left unchecked, processing must be triggered manually.

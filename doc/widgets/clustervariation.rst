Cluster Variation
=================

.. figure:: icons/ClusterVariation.png

Computes cluster-specific differential expression for all the genes in the input gene expression data set. The input data has to contain a meta feature (column) with labels of clusters.

Signals
-------

**Inputs**:

-  **Data**

   Single cell gene expression data with cells in rows and genes in columns. The data has to contain a meta feature (column) with a cluster identifier, where each cell has to a member of exactly one cluster.

**Outputs**:

-  **Data**

   A data matrix with genes in rows and clusters in columns. The entry of this matrix is a differential expression estimated as a logarithm of mean expression of a gene within the cluster of cells divided by mean expression of the same gene in cells outside the cluster.

Description
-----------

Input to a **Cluster Variation** widget is a data where cells have been clustered or provided some domain-specific label. The widget computes the differential expression of all the genes in the data set, where expression within the cluster of cells is compared to the expression outside the cluster.

.. figure:: images/ClusterVariation-stamped.png

1. Selection of a meta feature (column) in the input data that reports on cluster labels.
2. Tick to automatically process input data and send the result of scoring to the output. If left unchecked, processing must be triggered manually.

Change Log
==========

Version 0.8 (2018-08-07)
------------------------
##### New Features and Widgets
* Single Cell Preprocess: takes raw input data and does normalization, gene expression scaling, standardization, or binarization.
* Batch Effect Removal: uses model-based removal of batch effects, where batches are identified through possible combinations of continuous or discrete variables. Note that the widget transforms the input data.
* Align Datasets: uses canonical correlation analysis to embed cells within space possibly absent of batch effects. A faithful reimplementation of original batch effect removal as implemented in Seurat. Note that the widget performs embedding, not data transformation.
* Cluster Analysis: this widget is from Bioinformatics add-on. Provided the data clusters, it lists differentially expressed genes and associated pathways.
* Dot Matrix: plots subset of genes and their cluster-based expressions.

##### Notes
This is a major update of Single Cell Orange. Widgets can now mimic the pipeline from Seurat while retaining unique features such as combination with machine learning, gene selection, visualization, and differential expression analysis. Future releases will improve on documentation and example workflows.

Version 0.7 (2018-07-11)
------------------------
##### New Features and Widgets
* Load Data: a major upgrade of a Load Data widget from 0.6, now supports simultaneous loading of multiple data sets that are merged in a single data table. Tab or comma-separated and 10X Genomics Data files are supported.
* Gene Set Enrichment: now supports custom lists and classifications of gene markers. Use File to read the data with gene name and group marker, standardize gene names with Gene Name Matcher, and then feed the custom gene list into Gene Set Enrichment.
* Score Cells: includes new scoring methods that result in the clearer identification of cell types in embedding maps.
* Other improvements stem from improvements in Orange suite. Workflows can now be edited in separate windows. Saving of the data has been speed-up.

#### Performance Improvements
* t-SNE: improved checking of the input data to avoid unnecessary recomputation.

Version 0.6 (2018-05-31)
------------------------
##### New Features and Widgets
* Single Cell Datasets: now serves the data that includes NCBI gene IDs.
* Gene Name Matcher: a new widget that outputs the data table and associates genes with their NCBI IDs, either by adding a new meta column or providing meta information to column headers. The widget provides information on the success of matching the provided gene labels (names) with NCBI database.
* Simplified gene name matching across the platform: matching uses NCBI gene IDs. Used by widgets like Score Cells and Cluster Analysis, and all the widgets from Bioinformatics add-on including GO Browser and KEGG Pathways. Matching assumes that the data includes NCBI IDs. 
* Simplified interface in widgets like GO Browser and KEGG Pathways: datasets now carry information on organism and location of gene IDs, so widgets no longer need the part of the dialog to provide this type of information. We have removed this part of the dialog and simplified the interface of several widgets.
* Cluster Analysis: can now accept a list of marker genes to be displayed in the analysis.
* Gene Markers: all markers now include NCBI IDs, we have also extended a list of markers.

#### Bug Fixes
* Cluster Analysis: we fixed a bug in computation of over-expression scores. The widget now outputs the same scores as Differential Expression widget when using a hypergeometric test for scoring.

Version 0.5 (2018-05-08)
------------------------
##### New Features and Widgets
* Cluster Analysis: a prototype widget with a plot of most enriched genes per cluster. Bi-clustering of the gene-cluster association chart.
* Marker Genes: improved rendering of the marker lists, links to the data source, and synching of the markers with google docs spreadsheet to simplify modifications of the list.
* Score Cells: improved gene name matching in the widget, matching is no longer performed by name but through a library of gene name synonyms.
* A new example workflow for batch normalization is now included in scOrange, find it on the splash screen when starting scOrange or choose Help -> Welcome.

Version 0.4 (2018-04-16)
------------------------
##### New Features and Widgets
 * Normalize: batch normalization using a linear model.

Version 0.3 (2018-02-19)
------------------------
##### New Features and Widgets
 * Data Loader: loads gene expression data from tab-delimed files and separate gene and cell annotation files
 * Marker Genes: shows a list of marker genes for selected organism.


Version 0.2 (2018-01-25)
------------------------
##### New Features and Widgets
 * Louvain Clustering: cluster cells using a Louven clustering on cell similarity network.
 * Cell/Gene: filters cells or genes based on the frequency of non-zero entries in the data.


Version 0.1 (2017-12-01)
------------------------
##### New Features and Widgets
 * Single Cell Datasets: provides access to precompiled repository of the single cell data.

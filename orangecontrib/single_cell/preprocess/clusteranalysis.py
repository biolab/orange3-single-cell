import Orange
import numpy as np
import time
from functools import lru_cache

from orangecontrib.bioinformatics.widgets.utils.data import GENE_ID_ATTRIBUTE
from scipy.stats import hypergeom
from sklearn.cluster.bicluster import SpectralBiclustering
from Orange.data import Domain, DiscreteVariable, ContinuousVariable, Table


class ClusterAnalysis:
    """
    Analysis of single cell clusters based on enriched genes.

    Parameters
    ----------
    data : Orange table of cells with genes as attributes and each row has a defined cluster
    cluster_var : string of preferred row class_var


    Attributes
    ----------
    data : Orange table of cells with genes as attributes and each cell belongs to a cluster
    clusters_rows : list of indexes in which cluster belongs each cell
    clusters_ind : list of unique indexes that appear in data
    clusters_names : ordered list of all cluster names that appear in data
    genes : list of genes for which we want the output
    model : result table of percentage expressing
    o_model : result Orange.data.Table of percentage expressing
    row_order_ : order of clusters after spectral biclustering
    column_order_ : order of genes after spectral biclustering
    enriched_matrix_low : percentages for under-expression of genes
    enriched_matrix_high : percentages for over-expression of genes
    enriched_matrix : percentages for gene expression according to given 'enrichment' parameter
    enrichment : set of available enrichment parameters
    """

    def __init__(self, data, cluster_var='Cluster', callback=None):
        self.data = data
        self.gene_id_attribute = self.data.attributes.get(GENE_ID_ATTRIBUTE, None)

        self.X = self.data.X > 0

        if cluster_var is not None:
            self.class_var = self.data.domain[cluster_var]
        else:
            self.class_var = self.data.domain.class_var

        self.clusters_rows = self.data.get_column_view(self.class_var)[0]
        self.clusters_names = self.class_var.values
        self.clusters_ind = np.unique(self.clusters_rows)

        # take only those that appear in data
        self.clusters_names = [self.clusters_names[int(i)] for i in self.clusters_ind]
        self.columns = self.data.domain.attributes

        self.enriched_matrix = None
        self.enrichment = {'high', 'low', 'either'}
        self.enriched_matrix_low = None
        self.enriched_matrix_high = None
        self.model = None
        self.o_model = None
        self.row_order_ = None
        self.column_order_ = None
        self.genes = None

        self._create_enriched_matrix(callback=callback)

    def _create_enriched_matrix(self, callback=None):
        """
        Create matrix of how much each gene enriches each cluster.
        shape: (clusters_names, genes)
        """
        # get a list of indexes of all genes (NOTE: can cause bugs)
        genes = range(len(self.columns))

        Z = self.X

        # n, M for hypergeometric distribution
        count_all_pos = np.sum(Z, axis=0)  # n - all positive (number of success states in the population)
        count_all = Z.shape[0]  # M - all (population size)

        low = np.empty(shape=(len(self.clusters_names), len(genes)), dtype='float64')
        high = np.empty(shape=(len(self.clusters_names), len(genes)), dtype='float64')

        # for each cluster we calculate the n_enriched most enriched
        for c in range(len(self.clusters_names)):
            if callback is not None:
                callback(c / len(self.clusters_names))
            # get all cells that belong to cluster c
            cells = Z[self.clusters_rows == self.clusters_ind[c]]

            # k, N for hypergeometric distribution
            count_cluster_pos = np.sum(cells, axis=0)  # k - positive in cluster c ) number of observed successes)
            count_cluster = cells.shape[0]  # N - all in cluster c (number of draws)

            # calculate cdf for every gene
            low[c] = np.array(
                [hypergeom.cdf(count_cluster_pos[i], count_all, count_all_pos[i], count_cluster) for i in
                 genes])
            high[c] = hypergeom.sf(count_cluster_pos, count_all, count_all_pos, count_cluster) + \
                      hypergeom.pmf(count_cluster_pos, count_all, count_all_pos, count_cluster)

        self.enriched_matrix_low = low
        self.enriched_matrix_high = high

    def intersection(self, gene_list):
        """
        Get intersection of genes in the dataset and in the iterable.

        Parameters
        ----------
        gene_list : iterable
            Genes to find in the dataset.

        Returns
        -------
        list
            List of genes in the intersection.
        """
        gene_set = set()

        if gene_list:
            gene_set = set([str(gene) for gene in gene_list])

        return [gene.attributes[self.gene_id_attribute] for gene in self.columns
                if self.gene_id_attribute in gene.attributes
                and str(gene.attributes[self.gene_id_attribute]) in gene_set]

    @lru_cache(maxsize=3)
    def enriched_genes(self, gene_list, enrichment=None, biclustering=True, callback=None):
        """
        Cluster-based enrichment scores for a tuple of genes

        Parameters
        ----------
        gene_list : tuple of ints, gene ids
        enrichment : ignored, used for consistency in function call
        biclustering : boolean, return biclustering model
        """
        # fix enrichment at 'either', to calculate the shortest tail p-value
        enrichment = 'either'
        self.enriched_matrix = np.hstack((self.enriched_matrix_low, self.enriched_matrix_high))
        enumerated_attributes = enumerate(gene.attributes[self.gene_id_attribute] for gene in self.columns
                                          if self.gene_id_attribute in gene.attributes)
        self.genes = [attribute[0] for attribute in enumerated_attributes if attribute[1] in gene_list]

        return self._create_model(enrichment, biclustering, callback=callback)

    @lru_cache(maxsize=3)
    def enriched_genes_per_cluster(self, n=3, enrichment='high', biclustering=True, callback=None):
        """
        n genes that are most enriched for each cluster.

        Parameters
        ----------
        n : int, number of enriched genes per cluster
        enrichment : string, type of enrichment (high, low, either)
        biclustering : boolean, return biclustering model
        """
        res = list()

        # get a list of indexes of all genes
        genes = list(range(len(self.columns)))
        genes_enum = genes

        # handle enrichment matrix according to enrichment parameter
        if enrichment not in self.enrichment:
            raise ValueError("enrichment should be either 'high', 'low' or 'either'"
                             ", but a value %r was passed" %
                             enrichment)
        elif enrichment == 'high':
            self.enriched_matrix = self.enriched_matrix_high
        elif enrichment == 'low':
            self.enriched_matrix = self.enriched_matrix_low
        elif enrichment == 'either':
            self.enriched_matrix = np.hstack((self.enriched_matrix_low, self.enriched_matrix_high))
            genes.extend(genes)
            genes_enum = list(range(len(genes)))

        # for each cluster we calculate the n most enriched genes
        for c, enriched_percent in enumerate(self.enriched_matrix):
            # search amongst the genes that haven't bee selected yet
            enriched_percent = enriched_percent[genes_enum]

            # sorting
            zipped = list(zip(*sorted(zip(enriched_percent, genes), key=lambda x: x[0])))
            enriched_genes = list(zipped[1][:n])

            # remove found enriched genes from list of all genes (ensure every cluster has n unique enriched genes)
            to_remove_enum = enriched_genes + [g + len(self.columns) for g in enriched_genes]
            genes = [g for g in genes if g not in enriched_genes]
            genes_enum = [g for g in genes_enum if g not in to_remove_enum]

            res.extend(enriched_genes)

        self.genes = res
        return self._create_model(enrichment, biclustering, callback=callback)

    @lru_cache(maxsize=3)
    def enriched_genes_data(self, n=20, enrichment='high', biclustering=True, callback=None):
        """
        n top enriched genes, where "top" means for any cluster

        Parameters
        ----------
        n : int, number of enriched genes overall
        enrichment : string, type of enrichment (high, low, either)
        biclustering : boolean, return biclustering model
        """
        res = list()

        # get a list of indexes of all genes
        genes = list(range(len(self.columns)))

        # handle enrichment matrix according to enrichment parameter
        if enrichment not in self.enrichment:
            raise ValueError("enrichment should be either 'high', 'low' or 'either'"
                             ", but a value %r was passed" %
                             enrichment)
        elif enrichment == 'high':
            self.enriched_matrix = self.enriched_matrix_high
        elif enrichment == 'low':
            self.enriched_matrix = self.enriched_matrix_low
        elif enrichment == 'either':
            self.enriched_matrix = np.hstack((self.enriched_matrix_low, self.enriched_matrix_high))
            genes.extend(genes)

        # select n enriched genes for each cluster
        for enriched_percent in self.enriched_matrix:
            zipped = sorted(zip(enriched_percent, genes), key=lambda x: x[0])
            enriched_genes_tuple = zipped[:n]

            res.extend(enriched_genes_tuple)

        # sort whole list
        res = sorted(res, key=lambda x: x[0])

        # remove duplicate genes that are less enriched
        seen = set()
        res = [item for item in res if item[1] not in seen and not seen.add(item)]
        res = list(zip(*res))

        # take only n genes
        self.genes = res[1][:n]

        return self._create_model(enrichment, biclustering, callback=callback)

    def _fraction_expressing(self, enrichment, callback=None):
        """
        Expression percentage and p-values for each enriched gene for each cluster.

        Params
        ----------
        enrichment : string, used internally for p-value
        """
        Z = self.X
        res = list()
        pvalues = list()

        for c in range(len(self.clusters_names)):
            if callback is not None:
                callback(c/len(self.clusters_names) * 0.2)
            cells = Z[self.clusters_rows == self.clusters_ind[c]]
            # calculate fraction expressing
            res.append([sum(cells[:, gene] > 0) / cells.shape[0] if cells.shape[0] > 0 else 0 for gene in self.genes])

            # if 'either' - we choose the lower p-value one (higher wouldn't have been chosen beforehand)
            if enrichment == 'either':
                pvalues.append([min(self.enriched_matrix[c, gene],
                                    self.enriched_matrix[c, gene + len(self.columns)]) for gene in
                                self.genes])
            # else we take p-value normally
            else:
                pvalues.append([self.enriched_matrix[c, gene] for gene in self.genes])

        self.model = np.array(res)
        self.pvalues = np.array(pvalues)

    @staticmethod
    def biclustering(matrix, distance, callback=None):
        if min(matrix.shape) <= 2:
            return np.arange(matrix.shape[0]), np.arange(matrix.shape[1])

        best_score = np.iinfo(np.dtype('uint16')).max
        best_model = None

        # find the best biclusters (needs revision)
        limit = int(min(matrix.shape) / 2) - 1
        limit = 3 if limit < 3 else limit
        for i in range(2, limit):
            if callback is not None:
                callback(0.2 + (i - 2) / (limit - 2) * 0.8)
            # perform biclustering
            model = SpectralBiclustering(
                n_clusters=i, method='log', random_state=0)
            model.fit(matrix)
            fit_data = matrix[np.argsort(model.row_labels_)]
            fit_data = fit_data[:, np.argsort(model.column_labels_)]

            # calculate score and save the lowest one
            score = distance(fit_data)
            if score < best_score:
                best_score = score
                best_model = model

        return np.argsort(best_model.row_labels_), np.argsort(best_model.column_labels_)

    @staticmethod
    def reorder(matrix, row_order, column_order):
        return

    def _sort_fraction_expressing(self, callback=None):
        """
        Sort rows and columns based on expressed percentages.
        """

        self.row_order_, self.column_order_ = self.biclustering(self.model,
                                                                self.neighbor_distance,
                                                                callback)

        self.model = self.model[self.row_order_]
        self.model = self.model[:, self.column_order_]

        # order p-values matrix too
        self.pvalues = self.pvalues[self.row_order_]
        self.pvalues = self.pvalues[:, self.column_order_]

    def _create_model(self, enrichment, biclustering, callback):
        """
        Create cluster analysis model.

        Params
        ----------
        enrichment : string, used internally

        Return
        ----------
        res_rows : list of strings, clusters_names
        res_genes : list of strings, genes
        fraction_expressing : numpy matrix, fraction of expressing for enriched genes per cluster
        p-value : numpy matrix, p-values for enriched genes for each cluster
        """

        self._fraction_expressing(enrichment)

        if not biclustering:
            self.row_order_ = list(range(len(self.clusters_names)))
            self.column_order_ = list(range(len(self.genes)))
        else:
            self._sort_fraction_expressing(callback=callback)

        res_genes = list(np.ravel([self.columns[self.genes[i]].name for i in self.column_order_]))
        res_rows = [self.clusters_names[i] for i in self.row_order_]
        return res_rows, res_genes, self.model, self.pvalues

    @staticmethod
    def neighbor_distance(X):
        """
        Calculate euclicean distance between neighbors in rows and columns

        Params
        ----------
        X : (n_clusters, n_genes) : biclustered matrix of expressing percentage

        Return
        ----------
        res : sum of distances
        """
        # rows
        res = np.sum([np.linalg.norm(X[i - 1] - X[i]) for i in range(1, len(X))])
        # columns
        res += np.sum([np.linalg.norm(X[:, i - 1] - X[:, i]) for i in range(1, len(X[0]))])
        return res

    def draw_pyplot(self):
        """
        Draw with pyplot clusters_names on y and enriched genes on x, with size of circle representing
        percentage of expressed.
        """
        import matplotlib.pyplot as plt

        circle_mult = .4
        x_len = range(len(self.column_order_))
        y_len = range(len(self.row_order_))

        ## Plot
        fig, ax = plt.subplots()

        for i in range(len(self.model)):
            for j in range(len(self.model[0])):
                ax.add_artist(plt.Circle((j, i), self.model[i][j] * circle_mult))

        plt.xlim(-1, len(self.column_order_))
        plt.ylim(-1, len(self.row_order_))

        ## Y axis - clusters_names
        clusters_list = self.clusters_names
        clusters = [clusters_list[i] for i in self.row_order_]
        plt.yticks(y_len, clusters, rotation='horizontal')
        plt.margins(0.2)

        ## X axis - genes
        genes_list = self.genes
        genes = np.ravel(list([self.columns[genes_list[i]].name for i in self.column_order_]))

        plt.xticks(x_len, genes, rotation='vertical')
        plt.margins(0.2)

        plt.subplots_adjust(bottom=0.2)
        # plt.savefig('fig1.png')
        plt.show()
        return plt

    @staticmethod
    def contingency_table(data: np.array, rows: DiscreteVariable, columns: list, metas) -> Table:
        return Table.from_numpy(Domain([ContinuousVariable.make(column) for column in columns], metas=[rows]),
                                data, metas=metas)

    def create_contingency_table(self):
        """
        Create Orange.table from results

        Return
        --------
        o_model : Orange.Table
        """
        # create Orange.Table for calculated model
        dmn = [self.columns[self.genes[i]].name for i in self.column_order_]
        mts = DiscreteVariable.make(self.class_var.name, values=self.clusters_names)
        self.o_model = self.contingency_table(self.model, mts, dmn, np.array([self.row_order_]).T)
        return self.o_model


if __name__ == '__main__':
    # Example usages

    # data = Orange.data.Table('data/bone_marrow_louvain.pickle')
    data = Orange.data.Table('data/bone_marrow_louvain.pickle')

    start = time.time()
    print("START")
    CA = ClusterAnalysis(data, cluster_var='Cluster')
    print(time.time() - start)
    res = CA.enriched_genes(gene_list=('HBG1', 'S100A9', 'HBG2'))
    print(time.time() - start)
    CA.draw_pyplot()

    res = CA.enriched_genes_per_cluster(n=2, enrichment='either', biclustering=False)
    print(time.time() - start)
    CA.draw_pyplot()

    res = CA.enriched_genes_data(enrichment='low')
    print(time.time() - start)
    CA.draw_pyplot()

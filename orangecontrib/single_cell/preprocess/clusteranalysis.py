import Orange
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import hypergeom
from Orange.data import Domain, DiscreteVariable, ContinuousVariable, Table
from sklearn.cluster.bicluster import SpectralBiclustering


class ClusterAnalysis:
    """
    Analysis of single cell clusters based on enriched genes.

    Parameters
    ----------
    data : Orange table of cells with genes as attributes and each row has a defined cluster
    clustering_var : string of preferred row class_var


    Attributes
    ----------
    data : Orange table of cells with genes as attributes and each cell belongs to a cluster
    clusters_all : list of; in which cluster belongs each cell
    clusters : ordered list of all clusters
    genes : list of genes for which we want the output
    model : result table of percentage expressing
    o_model : result Orange.data.Table of percentage expressing
    row_order_ : order of clusters after spectral biclustering
    column_order_ : order of genes after spectral biclustering
    """

    def __init__(self, data, cluster_var='Cluster'):
        self.data = data

        self.data.X = self.data.X > 0

        if cluster_var is not None:
            self.class_var = self.data.domain[cluster_var]
        else:
            self.class_var = self.data.domain.class_var

        self.clusters_all = self.data.get_column_view(self.class_var)[0]
        self.clusters = self.class_var.values

        self.enriched_matrix = None
        self.model = None
        self.o_model = None
        self.row_order_ = None
        self.column_order_ = None
        self.genes = None

        self._create_enriched_matrix()

    def _create_enriched_matrix(self):
        """
        Create matrix of how much each gene enriches each cluster.
        shape: (clusters, genes)
        """
        # get a list of indexes of all genes
        genes = range(len(self.data.domain.attributes))

        Z = self.data.X

        # n, M for hypergeometric
        count_all_pos = np.sum(Z, axis=0)  # n - all positive
        count_all = Z.shape[0]  # M - all

        res = np.empty(shape=(len(self.clusters), len(genes)), dtype='float32')

        # for each cluster we calculate the n_enriched most enriched
        for i,c in enumerate(range(len(self.clusters))):
            # get all cells that belong to cluster c
            cells = Z[self.clusters_all == c]

            # k, N for hypergeometric
            count_cluster_pos = np.sum(cells, axis=0)  # k - positive in cluster c
            count_cluster = cells.shape[0]  # N - all in cluster c

            # calculate cdf for every gene
            res[i] = [hypergeom.cdf(count_cluster_pos[i], count_all, count_all_pos[i], count_cluster) for i in
                                genes]

        self.enriched_matrix = res

    def enriched_genes(self, gene_list):
        """
        Cluster-based enrichment scores for a list of genes

        Parameters
        ----------
        gene_list : list of strings, names of genes
        """

        enumereted_attributes = enumerate([gene.name for gene in self.data.domain.attributes])
        res = [attribute[0] for attribute in enumereted_attributes if attribute[1] in gene_list]

        """
        res = list()
        for gene in gene_list:
            for index, att_gene in enumerate(self.data.domain.attributes):
                if gene == att_gene.name:
                    res.append(index)
        """

        self.genes = res
        self.percentage_expressing()
        self.sort_percentage_expressing()


    def enriched_genes_per_cluster(self, n=3):
        """
        n genes that are most enriched for each cluster.

        Parameters
        ----------
        m : int, number of enriched genes per cluster
        """
        res = list()

        # get a list of indexes of all genes
        genes = range(len(self.data.domain.attributes))

        # for each cluster we calculate the n most enriched genes
        for enriched_percent in self.enriched_matrix:
            # search amongst the genes that haven't been selected yet
            enriched_percent = enriched_percent[genes]

            # slower sorting
            zipped = list(zip(*sorted(zip(enriched_percent, genes), key=lambda x: x[0])))
            enriched_genes = zipped[1][:n]

            # remove found enriched genes from list of all genes (ensure every cluster has 2 unique enriched genes)
            genes = list(set(genes) - set(enriched_genes))

            res.extend(enriched_genes)

        self.genes = res
        self.percentage_expressing()
        self.sort_percentage_expressing()

    def enriched_genes_data(self, n=20):
        """
        n top enriched genes, where "top" means for any cluster

        Parameters
        ----------
        m : int, number of enriched genes overall
        """
        res = list()

        # get a list of indexes of all genes
        genes = range(len(self.data.domain.attributes))

        # for each cluster we calculate the n most enriched genes
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

        self.genes = res[1][:n]
        self.percentage_expressing()
        self.sort_percentage_expressing()


    def percentage_expressing(self):
        """
        Expression percentage for each enriched gene for each cluster.
        """
        Z = self.data.X
        res = list()

        for c in range(len(self.clusters)):
            cells = Z[self.clusters_all == c]
            res.append([sum(cells[:, gene] > 0) / cells.shape[0] for gene in self.genes])
        self.model = np.array(res)

    def sort_percentage_expressing(self):
        """
        Sort clusters based on expressed percentages.

        Returns
        --------
        o_model : Orange table of biclustered expressed percentages.
        """
        if (len(self.genes) < 2):
            self.column_order_ = range(len(self.genes))
            self.row_order_ = range(len(self.clusters))
            return

        best_score = np.iinfo(np.dtype('uint16')).max
        best_fit = None
        best_model = None

        # find the best biclusters (needs revision)
        limit = int(min(len(self.genes), len(self.clusters)) / 2) - 1
        limit = 3 if limit < 3 else limit
        for i in range(2, limit):
            # perform biclustering
            model = SpectralBiclustering(
                n_clusters=i, method='log', random_state=0)
            model.fit(self.model)
            fit_data = self.model[np.argsort(model.row_labels_)]
            fit_data = fit_data[:, np.argsort(model.column_labels_)]

            # calculate score and safe the lowest one
            score = self._neighbor_distance(fit_data)
            if score < best_score:
                print(i, score)
                best_score = score
                best_fit = fit_data
                best_model = model

        self.model = best_fit
        self.row_order_ = np.argsort(best_model.row_labels_)
        self.column_order_ = np.argsort(best_model.column_labels_)

        # create orange table (robust)
        dmn = np.ravel(list([self.data.domain.attributes[self.genes[i]].name for i in self.column_order_]))
        mts = DiscreteVariable.make(self.class_var.name, values=self.clusters)
        self.o_model = Table.from_numpy(Domain([ContinuousVariable.make(column) for column in dmn], metas=[mts]),
                                        self.model, metas=np.array([self.row_order_]).T)
        return self.o_model

    def _neighbor_distance(self, X):
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
        Draw with pyplot clusters on y and enriched genes on x, with size of circle representing
        percentage of expressed.
        """
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

        ## Y axis - clusters
        clusters_list = self.clusters
        clusters = [clusters_list[i] for i in self.row_order_]
        plt.yticks(y_len, clusters, rotation='horizontal')
        plt.margins(0.2)

        ## X axis - genes
        genes_list = self.genes
        genes = np.ravel(list([self.data.domain.attributes[genes_list[i]].name for i in self.column_order_]))

        plt.xticks(x_len, genes, rotation='vertical')
        plt.margins(0.2)

        plt.subplots_adjust(bottom=0.2)
        #plt.savefig('fig1.png')
        plt.show()


if __name__ == '__main__':
    # Example usages
    data = Orange.data.Table('data/bone_marrow_louvain.pickle')

    start = time.time()
    print("START")
    CA = ClusterAnalysis(data, cluster_var='Cluster')
    #CA.enriched_genes(['HBG1', 'S100A9', 'HBG2'])
    #CA.enriched_genes_per_cluster(n=3)
    CA.enriched_genes_data()
    CA.draw_pyplot()

    end = time.time()
    print(end - start)

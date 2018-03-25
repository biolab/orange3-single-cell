import Orange
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.stats import hypergeom
from fastcluster import linkage
from scipy.spatial.distance import pdist, squareform


class SortDistanceMatrix:
    def seriation(Z, N, cur_index):
        '''
            Source: https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
            input:
                - Z is a hierarchical tree (dendrogram)
                - N is the number of points given to the clustering process
                - cur_index is the position in the tree for the recursive traversal
            output:
                - order implied by the hierarchical tree Z

            seriation computes the order implied by a hierarchical tree (dendrogram)
        '''
        if cur_index < N:
            return [cur_index]
        else:
            left = int(Z[cur_index - N, 0])
            right = int(Z[cur_index - N, 1])
            return (SortDistanceMatrix.seriation(Z, N, left) + SortDistanceMatrix.seriation(Z, N, right))

    def compute_serial_matrix(dist_mat, method="ward"):
        '''
            input:
                - dist_mat is a distance matrix
                - method = ["ward","single","average","complete"]
            output:
                - seriated_dist is the input dist_mat,
                  but with re-ordered rows and columns
                  according to the seriation, i.e. the
                  order implied by the hierarchical tree
                - res_order is the order implied by
                  the hierarhical tree
                - res_linkage is the hierarhical tree (dendrogram)

            compute_serial_matrix transforms a distance matrix into
            a sorted distance matrix according to the order implied
            by the hierarchical tree (dendrogram)
        '''
        N = len(dist_mat)
        flat_dist_mat = squareform(dist_mat)
        res_linkage = linkage(flat_dist_mat, method=method, preserve_input=True)
        res_order = SortDistanceMatrix.seriation(res_linkage, N, N + N - 2)
        seriated_dist = np.zeros((N, N))
        a, b = np.triu_indices(N, k=1)
        seriated_dist[a, b] = dist_mat[[res_order[i] for i in a], [res_order[j] for j in b]]
        seriated_dist[b, a] = seriated_dist[a, b]

        return seriated_dist, res_order, res_linkage


class ClusterAnalysis:
    """
    Analysis of single cell clusters based on enriched genes.

    Attributes
    ----------
    data : Orange table of already clustered cells with genes as attributes

    n_enriched : number of enriched genes per cluster
    """

    def __init__(self, data, n_enriched=1):
        self.data = data
        self.clusters = np.array([d["Cluster"].value for d in self.data])  # get a cluster for every row
        self.n_enriched = n_enriched  # number of most enriched genes per cluster
        self.enriched_genes = None
        self.enriched_genes = None
        self.model = None

    def gene_selection(self):

        res = OrderedDict()  # most enriched genes

        genes = np.array([d.name for d in self.data.domain.attributes])

        Z = self.data.X > 0  # Is gene expressed - discretize data
        count_all_pos = np.sum(Z, axis=0)  # n - all positive
        count_all = Z.shape[0]  # M - all

        for c in sorted(set(self.clusters)):
            cells = Z[self.clusters == c]  # get all cells that belong to cluster c

            count_cluster_pos = np.sum(cells, axis=0)  # k - positive in cluster
            count_cluster = cells.shape[0]  # N - all in cluster

            enriched_percent = [hypergeom.cdf(count_cluster_pos[i], count_all, count_all_pos[i], count_cluster) for i in
                                range(len(genes))]
            zipped = list(zip(*sorted(zip(enriched_percent, range(len(enriched_percent))), key=lambda x: x[0])))

            enriched_genes = zipped[1][:self.n_enriched]
            res[c] = enriched_genes

        self.enriched_genes = res

    def percentage_expressing(self):
        """
        Percentage of expressed for each enriched gene for each cluster.
        """
        Z = data.X > 0
        res = list()

        for c, _ in self.enriched_genes.items():
            # print(c)
            cells = Z[self.clusters == c]
            res.append(np.ravel([list(sum(cells[:, gene]) / cells.shape[0] for gene in genes) for _, genes in
                                 self.enriched_genes.items()]))
        self.model = np.array(res)

    def sort_percentage_expressing(self):
        """
        Sort clusters based on expressed percentages.
        """
        _, permutation, _ = SortDistanceMatrix.compute_serial_matrix(squareform(pdist(CA.model)))
        self.model = self.model[np.argsort(permutation)]

        temp = OrderedDict()
        enriched_genes = list(self.enriched_genes.items())
        for i in permutation:
            key, value = enriched_genes[i]
            temp[key] = value

        self.enriched_genes = temp

    def draw_pyplot(self):
        """
        Draw with pyplot clusters on y and enriched genes on x, with size of circle representing
        percentage of expressed.
        """
        circle_mult = .4
        x_len = range(len(self.model[0]))
        y_len = range(len(self.enriched_genes))

        genes = np.ravel([list(self.data.domain.attributes[gene].name for gene in genes) for _, genes in
                          self.enriched_genes.items()])

        ## Plot
        fig, ax = plt.subplots()

        for i in range(len(self.model)):
            for j in range(len(self.model[0])):
                ax.add_artist(plt.Circle((j, i), self.model[i][j] * circle_mult))

        plt.xlim(-1, len(self.model[0]))
        plt.ylim(-1, len(self.enriched_genes))

        ## Y axis
        plt.yticks(y_len, list(self.enriched_genes.keys()), rotation='horizontal')
        plt.margins(0.2)

        ## X axis
        plt.xticks(x_len, genes, rotation='vertical')
        plt.margins(0.2)

        plt.subplots_adjust(bottom=0.2)
        plt.savefig('fig1.png')
        plt.show()


if __name__ == '__main__':
    # np.set_printoptions(linewidth=160)
    # Example usage
    data = Orange.data.Table('a.pickle')
    CA = ClusterAnalysis(data, n_enriched=2)
    CA.gene_selection()
    CA.percentage_expressing()
    CA.sort_percentage_expressing()

    CA.draw_pyplot()
    print('DONE')

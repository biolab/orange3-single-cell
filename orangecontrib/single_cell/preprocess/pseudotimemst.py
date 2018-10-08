from math import ceil
import random

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr

import Orange
from Orange.distance import Euclidean


def add_tuple(t1, t2):
    return tuple(a+b for a,b in zip(t1,t2))

def farthest(x, tab_centroids, adj, p=-1, d=0):
    f = [farthest(z, tab_centroids, adj, d+dist(tab_centroids[x], tab_centroids[z])) for z in adj[x] if z!=p]
    if f: return max(f)
    else: return (d,x)


def order(x, adjp, pseudotime, tab_centroids, proj, proj_sign, adj, p=-1, d=0):
    for i in adjp[x]:
        if pseudotime[i] is not None: continue
        dxi = dist(tab_centroids[x], proj[i])
        if p==-1 and proj_sign[i]==-1:
            pseudotime[i] = d - dxi
        else:
            pseudotime[i] = d + dxi
    for z in adj[x]:
        if z==p: continue
        order(z, x, d+dist(tab_centroids[x], tab_centroids[z]))
    return pseudotime

def dist(coord1, coord2):
    return sum((x1-x2)**2 for x1, x2 in zip(coord1, coord2))**0.5


def mst_prim(data, dist=Euclidean()):
    n = len(data)
    dist = dist.fit(data)
    distance_to_tree = [float("inf") for i in range(n)]
    neighbour_in_tree = [None for i in range(n)]
    added = [False for i in range(n)]
    edges = []
    for r in range(n):  # repeat n times (add a new node each time)
        # find the node closest to the tree
        a = -1
        for i in range(n):
            if not added[i] and (a == -1 or distance_to_tree[i] < distance_to_tree[a]):
                a = i
        # add node a to the tree
        if neighbour_in_tree[a] is not None:
            edges.append((a, neighbour_in_tree[a]))
        added[a] = True
        # update distances of nodes to tree
        distance_a = dist(data[a], data)[0]  # distance function
        for i in range(n):
            if not added[i] and distance_a[i] < distance_to_tree[i]:
                distance_to_tree[i] = distance_a[i]
                neighbour_in_tree[i] = a

    return edges


class Pseudotimemst:

    def __init__(self, n_components=2, n_clusters=5,
                 projection=Orange.projection.PCA, clustering=Orange.clustering.KMeans): # can add projection, clustering alg
        self.ncomp = n_components
        self.data = None
        self.transformed_data = None
        self.nclust = n_clusters
        self.projection = projection
        self.clustering = clustering
        self.clustering_model = None
        self.pseudotime = None
        self.proj = None

    def fit_transform(self, data):
        self.data = data
        # transform
        projection_obj = self.projection(n_components=self.ncomp)
        self.transformed_data = projection_obj(self.data)(self.data)

        # cluster
        if self.nclust is None:
            self.nclust = ceil(0.05*len(self.transformed_data))
        clustering_model = self.clustering(n_clusters=self.nclust)(self.transformed_data)
        self.clustering_model = clustering_model

        self.clusters = [[] for _ in range(clustering_model.k)]
        for i, c in enumerate(clustering_model.labels_):
            self.clusters[c].append(i)

        # minimum spanning tree
        edges = mst_prim(Orange.data.Table.from_numpy(None, clustering_model.centroids))
        adj = [[] for _ in range(clustering_model.k)]
        for a,b in edges:
            adj[a].append(b)
            adj[b].append(a)
        mst_coordinates = [[clustering_model.centroids[a], clustering_model.centroids[b]] for a,b in edges]
        self.mst_coordinates = mst_coordinates

        # project
        proj = np.zeros((len(self.data), self.ncomp))
        proj_sign = [1 for _ in range(len(self.data))]
        adjp = [[] for _ in range(clustering_model.k)]

        tab_centroids = Orange.data.Table.from_numpy(None, clustering_model.centroids)
        for c in range(clustering_model.k):  # for every cluster
            for i in self.clusters[c]:  # for each instance in cluster c
                # find nearest adjacent cluster c2 (second nearest cluster)
                c2 = min(adj[c], key=lambda c2: dist(self.transformed_data[i].x, clustering_model.centroids[c2]))
                # project instance i onto edge (c, c2)
                vi = [self.transformed_data[i][d] - tab_centroids[c][d] for d in range(self.ncomp)]
                v = [tab_centroids[c2][d] - tab_centroids[c][d] for d in range(self.ncomp)]
                a = np.dot(vi, v) / np.dot(v, v)
                proj_sign[i] = 1
                if a < 0:
                    if len(adj[c]) > 1:
                        a = 0
                    else:
                        proj_sign[i] = -1
                if a > 1: a = 1
                proj[i] = [tab_centroids[c][d] + a * v[d] for d in range(self.ncomp)]
                adjp[c].append(i)
                adjp[c2].append(i)
        self.proj = proj

        def farthest(x, p=-1, d=0):
            f = [farthest(z, x, d + dist(tab_centroids[x], tab_centroids[z])) for z in adj[x] if
                 z != p]
            if f:
                return max(f)
            else:
                return (d, x)

        def order(x, p=-1, d=0):
            for i in adjp[x]:
                if pseudotime[i] is not None: continue
                dxi = dist(tab_centroids[x], proj[i])
                if p == -1 and proj_sign[i] == -1:
                    pseudotime[i] = d - dxi
                else:
                    pseudotime[i] = d + dxi
            for z in adj[x]:
                if z == p: continue
                order(z, x, d + dist(tab_centroids[x], tab_centroids[z]))

        start = farthest(0)[1]
        pseudotime = [None for i in range(len(self.transformed_data))]
        order(start)
        ptm, ptM = min(pseudotime), max(pseudotime)
        pseudotime = np.array([(pt - ptm) / (ptM - ptm) for pt in pseudotime])
        self.pseudotime = pseudotime

    def visualize(self):
        import matplotlib.pyplot as plt
        import matplotlib.colors as clr

        ax = plt.axes()
        ax.set_aspect('equal')
        ax.set_adjustable('box')
        for a, b in self.mst_coordinates:
            ax.plot([a[0], b[0]], [a[1], b[1]], c="black")

        x = self.clustering_model.centroids[:, 0]
        y = self.clustering_model.centroids[:, 1]
        ax.scatter(x, y, s=10, c="black", marker='x')

        x = np.array([row[0] for row in self.transformed_data])
        y = np.array([row[1] for row in self.transformed_data])
        l = np.array([str(row.get_class()) for row in self.transformed_data])

        colors = list(clr.CSS4_COLORS)
        random.seed(1234)
        random.shuffle(colors)

        labels = sorted(set(l))
        for i, label in enumerate(labels):
            mask = (l == label)
            ax.scatter(x[mask], y[mask], s=5, color=colors[i], label=label)

        for i in range(len(self.transformed_data)):
            ax.annotate("%.2f" % self.pseudotime[i], (x[i], y[i]), fontsize=7)

        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", markerscale=3)
        ax.autoscale()

        # projection lines
        for i in range(len(self.transformed_data)):
            ax.plot([x[i], self.proj[i][0]], [y[i], self.proj[i][1]], linewidth=1, color='lightgray')

        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    pmst = Pseudotimemst()

    tab = Orange.data.Table("brown-selected")

    pmst.fit_transform(tab)
    pmst.visualize()

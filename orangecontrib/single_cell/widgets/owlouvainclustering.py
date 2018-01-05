from typing import Optional

import networkx as nx
import numpy as np
from AnyQt.QtCore import Qt, QThread

from Orange.data import Table, DiscreteVariable
from Orange.distance import Euclidean
from Orange.misc import DistMatrix
from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting
from Orange.widgets.utils.annotated_data import get_next_name, add_columns
from Orange.widgets.utils.concurrent import FutureWatcher, ThreadExecutor
from Orange.widgets.utils.signals import Input, Output
from orangecontrib.single_cell.widgets.louvain import best_partition


def jaccard(x, y):
    # type: (set, set) -> float
    """Compute the Jaccard similarity between two sets."""
    return len(x & y) / len(x | y)


def table_to_graph(data, k_neighbours, distances, progress_callback=None):
    """Convert tabular data to a graph using a nearest neighbours approach with
    the Jaccard similarity as the edge weights.

    Parameters
    ----------
    data : Table
    k_neighbours : int
    distances : DistMatrix
    progress_callback : Callable[[float], None]

    Returns
    -------
    nx.Graph

    """
    # We do k + 1 because each point is closest to itself, which is not useful
    nearest_neighbours = np.argsort(distances, axis=1)[:, 1:k_neighbours + 1]
    # Convert to list of sets so jaccard can be computed efficiently
    nearest_neighbours = list(map(set, nearest_neighbours))
    num_nodes = len(nearest_neighbours)

    progress = 0

    # Create an empty graph and add all the data ids as nodes for easy mapping
    graph = nx.Graph()
    graph.add_nodes_from(range(len(data)))

    for idx, node in enumerate(graph.nodes):
        # Make sure to update the progress only when there is a change
        if int(100 * (idx + 1) / num_nodes) > progress:
            progress += 1
            progress_callback(progress)

        for neighbour in nearest_neighbours[node]:
            # graph.add_edge(node, neighbour)
            graph.add_edge(node, neighbour, weight=jaccard(
                nearest_neighbours[node], nearest_neighbours[neighbour]))

    return graph


class OWLouvainClustering(widget.OWWidget):
    name = 'Louvain Clustering'

    class Inputs:
        data = Input('Data', Table, default=True)

    class Outputs:
        annotated_data = Output('Annotated Data', Table, default=True)

    k_neighbours = Setting(30)
    resolution = Setting(1.)
    auto_apply = Setting(True)

    def __init__(self):
        super().__init__()

        self.data = None  # type: Optional[Table]
        self.distances = None  # type: Optional[DistMatrix]
        self.graph = None  # type: Optional[nx.Graph]
        self.partition = None  # type: Optional[np.array]

        self.__executor = ThreadExecutor(parent=self)

        graph_box = gui.vBox(self.controlArea, 'Graph parameters')
        self.k_neighbours_spin = gui.spin(
            graph_box, self, 'k_neighbours', minv=1, maxv=200,
            label='k neighbours', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_graph)
        self.cls_epsilon_spin = gui.spin(
            graph_box, self, 'resolution', 0, 5., 1e-2, spinType=float,
            label='Resolution', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_partition)
        self.compute_btn = gui.button(
            graph_box, self, 'Run clustering', callback=self.compute_partition)

    def _compute_distances(self):
        def _distances():
            if self.distances is None:
                self.setBlocking(True)
                self.setStatusMessage('Computing distances...')
                self.distances = Euclidean(self.data)

        return self.__executor.submit(_distances)

    def _build_graph(self):
        def _process():
            self.setBlocking(True)
            self.setStatusMessage('Building graph...')
            return table_to_graph(self.data, k_neighbours=self.k_neighbours,
                                  distances=self.distances)

        return self.__executor.submit(_process)

    def _compute_partition(self):
        def _process():
            self.setBlocking(True)
            self.setStatusMessage('Detecting communities...')

            partition = best_partition(self.graph, resolution=self.resolution)
            return np.fromiter(list(zip(*sorted(partition.items())))[1], dtype=int)

        return self.__executor.submit(_process)

    def _invalidate_distances(self):
        self.distances = None
        self._invalidate_graph()

    def _invalidate_graph(self):
        self.graph = None
        self._invalidate_partition()

    def _invalidate_partition(self):
        self.partition = None

    def compute_partition(self):
        if not self.data:
            return

        def _graph_node_progress(percentage_done):
            self.progressBarSet(percentage_done / 2)

        def _process():
            self.setBlocking(True)

            # Only init the progress bar if there is anything to do
            if self.graph is None or self.partition is None:
                self.progressBarInit()

            if self.graph is None:
                self.setStatusMessage('Building graph...')
                self.graph = table_to_graph(
                    self.data, k_neighbours=self.k_neighbours,
                    distances=self.distances,
                    progress_callback=_graph_node_progress)

            if self.partition is None:
                self.setStatusMessage('Detecting communities...')
                partition = best_partition(self.graph, resolution=self.resolution)
                self.partition = np.fromiter(list(
                    zip(*sorted(partition.items())))[1], dtype=int)

            self.progressBarFinished()

        def _done(f):
            assert QThread.currentThread() is self.thread()
            assert f.done()

            self.setStatusMessage('')
            self.setBlocking(False)

            domain = self.data.domain
            cluster_var = DiscreteVariable(
                get_next_name(domain, 'Cluster'),
                values=['C%d' % i for i in np.unique(self.partition)]
            )

            new_domain = add_columns(domain, metas=[cluster_var])
            new_table = self.data.transform(new_domain)
            new_table.get_column_view(cluster_var)[0][:] = self.partition
            self.Outputs.annotated_data.send(new_table)

        future = self.__executor.submit(_process)
        watcher = FutureWatcher(future, parent=self)
        watcher.done.connect(_done)

    @Inputs.data
    def set_data(self, data):
        self._invalidate_distances()

        self.data = data

        def _done(f):
            assert QThread.currentThread() is self.thread()
            assert f.done()

            self.setStatusMessage('')
            self.setBlocking(False)

        future = self._compute_distances()
        watcher = FutureWatcher(future, parent=self)
        watcher.done.connect(_done)


if __name__ == '__main__':
    from AnyQt.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    ow = OWLouvainClustering()
    ow.resetSettings()

    ow.set_data(Table(sys.argv[1] if len(sys.argv) > 1 else 'iris'))
    ow.show()
    app.exec_()

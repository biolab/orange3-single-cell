from typing import Optional

import networkx as nx
import numpy as np
from AnyQt.QtCore import Qt, QThread, pyqtSignal as Signal
from AnyQt.QtWidgets import QSlider, QPushButton

from Orange.data import Table, DiscreteVariable
from Orange.distance import Euclidean
from Orange.misc import DistMatrix
from Orange.projection import PCA
from Orange.widgets import widget, gui
from Orange.widgets.gui import SpinBoxWFocusOut
from Orange.widgets.settings import Setting
from Orange.widgets.utils.annotated_data import get_next_name, add_columns, \
    ANNOTATED_DATA_SIGNAL_NAME
from Orange.widgets.utils.concurrent import FutureWatcher, ThreadExecutor
from Orange.widgets.utils.signals import Input, Output
from orangecontrib.single_cell.widgets.louvain import best_partition

try:
    from orangecontrib.network.network import Graph
except:
    Graph = None


_MAX_PCA_COMPONENTS = 50


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
            if progress_callback:
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
        annotated_data = Output(ANNOTATED_DATA_SIGNAL_NAME, Table, default=True)
        graph = Output('Network', Graph)

    pca_components = Setting(25)
    k_neighbours = Setting(30)
    resolution = Setting(1.)
    auto_apply = Setting(True)

    pca_projection_complete = Signal()
    distances_complete = Signal()
    build_graph_complete = Signal()
    find_partition_complete = Signal()

    def __init__(self):
        super().__init__()

        self.data = None  # type: Optional[Table]
        self.distances = None  # type: Optional[DistMatrix]
        self.graph = None  # type: Optional[nx.Graph]
        self.partition = None  # type: Optional[np.array]

        self.__executor = ThreadExecutor(parent=self)

        pca_box = gui.vBox(self.controlArea, 'PCA Preprocessing')
        self.pca_components_slider = gui.hSlider(
            pca_box, self, 'pca_components', label='Components: ', minValue=2,
            maxValue=_MAX_PCA_COMPONENTS,
            callback=self._invalidate_pca_projection,
        )  # type: QSlider

        graph_box = gui.vBox(self.controlArea, 'Graph parameters')
        self.k_neighbours_spin = gui.spin(
            graph_box, self, 'k_neighbours', minv=1, maxv=200,
            label='k neighbours', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_graph,
        )  # type: SpinBoxWFocusOut
        self.cls_epsilon_spin = gui.spin(
            graph_box, self, 'resolution', 0, 5., 1e-2, spinType=float,
            label='Resolution', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_partition,
        )  # type: SpinBoxWFocusOut
        self.compute_btn = gui.button(
            graph_box, self, 'Run clustering', callback=self.compute_partition,
        )  # type: QPushButton

        # Connect the pipeline together
        self.pca_projection_complete.connect(self._compute_distances)
        self.distances_complete.connect(self._compute_graph)
        self.build_graph_complete.connect(self._compute_partition)
        self.find_partition_complete.connect(self._send_data)
        self.find_partition_complete.connect(self._processing_complete)

    def _compute_pca_projection(self):
        def _process():
            if self.pca_projection is None:
                self.setStatusMessage('Computing PCA...')
                self.setBlocking(True)

                pca = PCA(n_components=self.pca_components, random_state=0)
                model = pca(self.data)
                self.pca_projection = model(self.data)

            self.pca_projection_complete.emit()

        self.__executor.submit(_process)

    def _compute_distances(self):
        def _process():
            if self.distances is None:
                self.setStatusMessage('Computing distances...')
                self.setBlocking(True)

                self.distances = Euclidean(self.data)

            self.distances_complete.emit()

        self.__executor.submit(_process)

    def _compute_graph(self):
        def _process():
            if self.graph is None:
                self.setStatusMessage('Building graph...')
                self.setBlocking(True)

                self.graph = table_to_graph(
                    self.data, k_neighbours=self.k_neighbours,
                    distances=self.distances)

            self.build_graph_complete.emit()

        self.__executor.submit(_process)

    def _compute_partition(self):
        def _process():
            if self.partition is None:
                self.setStatusMessage('Detecting communities...')
                self.setBlocking(True)

                partition = best_partition(self.graph, resolution=self.resolution)
                self.partition = np.fromiter(list(
                    zip(*sorted(partition.items())))[1], dtype=int)

            self.find_partition_complete.emit()

        self.__executor.submit(_process)

    def _processing_complete(self):
        self.setStatusMessage('')
        self.setBlocking(False)

    def _send_data(self):
        domain = self.data.domain
        cluster_var = DiscreteVariable(
            get_next_name(domain, 'Cluster'),
            values=['C%d' % i for i in np.unique(self.partition)]
        )

        new_domain = add_columns(domain, metas=[cluster_var])
        new_table = self.data.transform(new_domain)
        new_table.get_column_view(cluster_var)[0][:] = self.partition
        self.Outputs.annotated_data.send(new_table)

        if Graph is not None:
            graph = Graph(self.graph)
            graph.set_items(new_table)
            self.Outputs.graph.send(graph)

    def _invalidate_pca_projection(self):
        self.pca_projection = None
        self._invalidate_distances()

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

        self._compute_pca_projection()

    @Inputs.data
    def set_data(self, data):
        self._invalidate_pca_projection()

        self.data = data

        # Can't have more PCA components than the number of attributes
        n_attrs = len(data.domain.attributes)
        self.pca_components_slider.setMaximum(min(_MAX_PCA_COMPONENTS, n_attrs))


if __name__ == '__main__':
    from AnyQt.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    ow = OWLouvainClustering()
    ow.resetSettings()

    ow.set_data(Table(sys.argv[1] if len(sys.argv) > 1 else 'iris'))
    ow.show()
    app.exec_()

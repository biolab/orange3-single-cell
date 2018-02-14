from collections import deque
from concurrent.futures import Future
from enum import Enum
from types import SimpleNamespace as namespace
from typing import Optional

import networkx as nx
import numpy as np
from AnyQt.QtCore import Qt, pyqtSignal as Signal, QObject
from AnyQt.QtWidgets import QSlider, QPushButton, QCheckBox, QWidget
from sklearn.neighbors import NearestNeighbors

from Orange.data import Table, DiscreteVariable
from Orange.projection import PCA
from Orange.widgets import widget, gui
from Orange.widgets.settings import DomainContextHandler, ContextSetting, \
    Setting
from Orange.widgets.utils.annotated_data import get_next_name, add_columns, \
    ANNOTATED_DATA_SIGNAL_NAME
from Orange.widgets.utils.concurrent import ThreadExecutor
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.widget import Msg
from orangecontrib.single_cell.widgets.louvain import best_partition

try:
    from orangecontrib.network.network import Graph
except:
    Graph = None


_MAX_PCA_COMPONENTS = 50
_DEFAULT_PCA_COMPONENTS = 25
_MAX_K_NEIGBOURS = 200
_DEFAULT_K_NEIGHBOURS = 30


METRICS = [('Euclidean', 'l2'), ('Manhattan', 'l1')]


def jaccard(x, y):
    # type: (set, set) -> float
    """Compute the Jaccard similarity between two sets."""
    return len(x & y) / len(x | y)


def table_to_graph(data, k_neighbours, metric, progress_callback=None):
    """Convert tabular data to a graph using a nearest neighbours approach with
    the Jaccard similarity as the edge weights.

    Parameters
    ----------
    data : Table
    k_neighbours : int
    metric : str
        A distance metric supported by sklearn.
    progress_callback : Callable[[float], None]

    Returns
    -------
    nx.Graph

    """
    # We do k + 1 because each point is closest to itself, which is not useful
    knn = NearestNeighbors(n_neighbors=k_neighbours, metric=metric).fit(data.X)
    nearest_neighbours = knn.kneighbors(data.X, return_distance=False)
    # Convert to list of sets so jaccard can be computed efficiently
    nearest_neighbours = list(map(set, nearest_neighbours))
    num_nodes = len(nearest_neighbours)

    # Create an empty graph and add all the data ids as nodes for easy mapping
    graph = nx.Graph()
    graph.add_nodes_from(range(len(data)))

    for idx, node in enumerate(graph.nodes):
        if progress_callback:
            progress_callback(idx / num_nodes)

        for neighbour in nearest_neighbours[node]:
            graph.add_edge(node, neighbour, weight=jaccard(
                nearest_neighbours[node], nearest_neighbours[neighbour]))

    return graph


class TaskQueue(QObject):
    """Not really a task queue `per-se`. Running start will run the tasks in
    the current list and cannot handle adding other tasks while running."""
    on_exception = Signal(Exception)
    on_complete = Signal()
    on_progress = Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.__tasks = deque()
        self.__progress = 0

    def push(self, task):
        self.__tasks.append(task)

    def __set_progress(self, progress):
        # Only emit progress signal when the progress has changed sufficiently
        if int(progress * 100) > int(self.__progress * 100):
            self.on_progress.emit(progress)
        self.__progress = progress

    def start(self):
        num_tasks = len(self.__tasks)

        for idx, task_spec in enumerate(self.__tasks):

            def __task_progress(percentage):
                current_progress = idx / num_tasks
                # How much progress can each task contribute to the total
                # work to be done
                task_percentage = len(self.__tasks) ** -1
                # Convert the progress done by the task into the total
                # progress to the task
                relative_progress = task_percentage * percentage
                self.__set_progress(current_progress + relative_progress)

            try:
                if getattr(task_spec, 'progress_callback', False):
                    task_spec.task(progress_callback=__task_progress)
                else:
                    task_spec.task()
                self.__set_progress((idx + 1) / num_tasks)

            except Exception as e:
                self.on_exception.emit(e)
                break

        self.on_complete.emit()


class OWLouvainClustering(widget.OWWidget):
    name = 'Louvain Clustering'
    icon = 'icons/LouvainClustering.svg'
    priority = 2150

    want_main_area = False

    settingsHandler = DomainContextHandler()

    class Inputs:
        data = Input('Data', Table, default=True)

    if Graph is not None:
        class Outputs:
            annotated_data = Output(ANNOTATED_DATA_SIGNAL_NAME, Table, default=True)
            graph = Output('Network', Graph)
    else:
        class Outputs:
            annotated_data = Output(ANNOTATED_DATA_SIGNAL_NAME, Table, default=True)

    apply_pca = ContextSetting(True)
    pca_components = ContextSetting(_DEFAULT_PCA_COMPONENTS)
    metric_idx = ContextSetting(0)
    k_neighbours = ContextSetting(_DEFAULT_K_NEIGHBOURS)
    resolution = ContextSetting(1.)
    auto_commit = Setting(True)

    class Error(widget.OWWidget.Error):
        data_has_nans = Msg(
            'Data has missing values. Please impute the missing values before '
            'continuing\n'
        )
        general_error = Msg("Error occured during clustering\n{}")

    class State(Enum):
        Pending, Running = range(2)

    def __init__(self):
        super().__init__()

        self.data = None  # type: Optional[Table]
        self.graph = None  # type: Optional[nx.Graph]
        self.partition = None  # type: Optional[np.array]

        self.__executor = ThreadExecutor(parent=self)
        self.__future = None  # type: Optional[Future]
        self.__state = self.State.Pending

        pca_box = gui.vBox(self.controlArea, 'PCA Preprocessing')
        self.apply_pca_cbx = gui.checkBox(
            pca_box, self, 'apply_pca', label='Apply PCA preprocessing',
            callback=self._invalidate_graph,
        )  # type: QCheckBox
        self.pca_components_slider = gui.hSlider(
            pca_box, self, 'pca_components', label='Components: ', minValue=2,
            maxValue=_MAX_PCA_COMPONENTS,
            callback=self._invalidate_pca_projection,
        )  # type: QSlider

        graph_box = gui.vBox(self.controlArea, 'Graph parameters')
        self.metric_combo = gui.comboBox(
            graph_box, self, 'metric_idx', label='Distance metric',
            items=[m[0] for m in METRICS], callback=self._invalidate_graph,
            orientation=Qt.Horizontal,
        )  # type: gui.OrangeComboBox
        self.k_neighbours_spin = gui.spin(
            graph_box, self, 'k_neighbours', minv=1, maxv=_MAX_K_NEIGBOURS,
            label='k neighbours', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_graph,
        )  # type: gui.SpinBoxWFocusOut
        self.cls_epsilon_spin = gui.spin(
            graph_box, self, 'resolution', 0, 5., 1e-2, spinType=float,
            label='Resolution', controlWidth=80, alignment=Qt.AlignRight,
            callback=self._invalidate_partition,
        )  # type: gui.SpinBoxWFocusOut

        self.apply_button = gui.auto_commit(
            self.controlArea, self, 'auto_commit', 'Apply', box=False,
            commit=lambda: self.commit(force=True), callback=self.commit,
        )  # type: QWidget

    def _compute_pca_projection(self):
        if self.pca_projection is None and self.apply_pca:
            self.setStatusMessage('Computing PCA...')

            pca = PCA(n_components=self.pca_components, random_state=0)
            model = pca(self.data)
            self.pca_projection = model(self.data)

    def _compute_graph(self, progress_callback=None):
        if self.graph is None:
            self.setStatusMessage('Building graph...')

            data = self.pca_projection if self.apply_pca else self.data

            self.graph = table_to_graph(
                data, k_neighbours=self.k_neighbours,
                metric=METRICS[self.metric_idx][1],
                progress_callback=progress_callback,
            )

    def _compute_partition(self):
        if self.partition is None:
            self.setStatusMessage('Detecting communities...')
            self.setBlocking(True)

            partition = best_partition(self.graph, resolution=self.resolution)
            self.partition = np.fromiter(list(zip(*sorted(partition.items())))[1], dtype=int)

    def _processing_complete(self):
        self.setStatusMessage('')
        self.setBlocking(False)
        self.progressBarFinished()

    def _handle_exceptions(self, ex):
        self.Error.general_error(str(ex))

    def cancel(self):
        """Cancel any running jobs."""
        if self.__state == self.State.Running:
            assert self.__future is not None
            self.__future.cancel()
            self.__future = None

        self.__state = self.State.Pending

    def commit(self, force=False):
        self.Error.clear()
        # Kill any running jobs
        self.cancel()
        assert self.__state == self.State.Pending

        if self.data is None:
            return

        # We commit if auto_commit is on or when we force commit
        if not self.auto_commit and not force:
            return

        if np.any(np.isnan(self.data.X)):
            self.Error.data_has_nans()
            return

        # Prepare the tasks to run
        queue = TaskQueue(parent=self)

        if self.pca_projection is None and self.apply_pca:
            queue.push(namespace(task=self._compute_pca_projection))

        if self.graph is None:
            queue.push(namespace(task=self._compute_graph, progress_callback=True))

        if self.partition is None:
            queue.push(namespace(task=self._compute_partition))

        # Prepare callbacks
        queue.on_progress.connect(lambda val: self.progressBarSet(100 * val))
        queue.on_complete.connect(self._processing_complete)
        queue.on_complete.connect(self._send_data)
        queue.on_exception.connect(self._handle_exceptions)

        # Run the task queue
        self.progressBarInit()
        self.setBlocking(True)
        self.__future = self.__executor.submit(queue.start)
        self.__state = self.State.Running

    def _send_data(self):
        domain = self.data.domain
        # Compute the frequency of each cluster index
        counts = np.bincount(self.partition)
        indices = np.argsort(counts)[::-1]
        index_map = {n: o for n, o in zip(indices, range(len(indices)))}
        new_partition = list(map(index_map.get, self.partition))

        cluster_var = DiscreteVariable(
            get_next_name(domain, 'Cluster'),
            values=['C%d' % (i + 1) for i, _ in enumerate(np.unique(new_partition))]
        )

        new_domain = add_columns(domain, metas=[cluster_var])
        new_table = self.data.transform(new_domain)
        new_table.get_column_view(cluster_var)[0][:] = new_partition
        self.Outputs.annotated_data.send(new_table)

        if Graph is not None:
            graph = Graph(self.graph)
            graph.set_items(new_table)
            self.Outputs.graph.send(graph)

    def _invalidate_pca_projection(self):
        self.pca_projection = None
        self._invalidate_graph()

    def _invalidate_graph(self):
        self.graph = None
        self._invalidate_partition()

    def _invalidate_partition(self):
        self.partition = None
        self.commit()

    @Inputs.data
    def set_data(self, data):
        self.closeContext()
        self.Error.clear()
        self.data = data
        self._invalidate_pca_projection()
        self.Outputs.annotated_data.send(None)
        self.Outputs.graph.send(None)
        self.openContext(self.data)

        if self.data is None:
            return

        # Can't have more PCA components than the number of attributes
        n_attrs = len(data.domain.attributes)
        self.pca_components_slider.setMaximum(min(_MAX_PCA_COMPONENTS, n_attrs))
        self.pca_components_slider.setValue(min(_DEFAULT_PCA_COMPONENTS, n_attrs))
        # Can't have more k neighbours than there are data points
        self.k_neighbours_spin.setMaximum(min(_MAX_K_NEIGBOURS, len(data) - 1))
        self.k_neighbours_spin.setValue(min(_DEFAULT_K_NEIGHBOURS, len(data) - 1))

        self.commit()

    def onDeleteWidget(self):
        self.cancel()
        super().onDeleteWidget()


if __name__ == '__main__':
    from AnyQt.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    ow = OWLouvainClustering()
    ow.resetSettings()

    ow.set_data(Table(sys.argv[1] if len(sys.argv) > 1 else 'iris'))
    ow.show()
    app.exec_()

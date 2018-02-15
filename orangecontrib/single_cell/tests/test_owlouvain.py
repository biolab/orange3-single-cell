import numpy as np

from Orange.data import Table, Domain, ContinuousVariable
from Orange.widgets.tests.base import WidgetTest
from orangecontrib.single_cell.widgets.owlouvainclustering import \
    OWLouvainClustering


# Deterministic tests
np.random.seed(42)


class TestOWLouvain(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(
            OWLouvainClustering, stored_settings={'auto_commit': False}
        )

    def tearDown(self):
        self.widget.onDeleteWidget()
        super().tearDown()

    def test_clusters_ordered_by_size(self):
        """Cluster names should be sorted based on the number of instances."""
        x1 = np.array([[0, 0]] * 20)
        x2 = np.array([[1, 0]] * 15)
        x3 = np.array([[0, 1]] * 10)
        x4 = np.array([[1, 1]] * 5)
        data = np.vstack((x1, x2, x3, x4))
        # Remove any order depencence in data, not that this should affect it
        np.random.shuffle(data)

        table = Table.from_numpy(domain=Domain.from_numpy(X=data), X=data)
        self.send_signal(self.widget.Inputs.data, table)
        self.widget.k_neighbours = 4
        self.widget.commit(force=True)
        output = self.get_output(self.widget.Outputs.annotated_data, wait=1000)

        clustering = output.get_column_view('Cluster')[0].astype(int)
        counts = np.bincount(clustering)
        np.testing.assert_equal(counts, sorted(counts, reverse=True))

    def test_missing_values_with_no_pca_preprocessing(self):
        data = np.ones((5, 5))
        data[range(5), range(5)] = np.nan
        np.random.shuffle(data)

        table = Table.from_numpy(domain=Domain.from_numpy(X=data), X=data)
        self.send_signal(self.widget.Inputs.data, table)
        self.widget.apply_pca = False
        self.widget.commit(force=True)

        self.assertTrue(self.widget.Error.data_has_nans.is_shown())

    def test_empty_dataset(self):
        # Prepare a table with 5 rows with only meta attributes
        meta = np.array([0] * 5)
        meta_var = ContinuousVariable(name='meta_var')
        table = Table.from_domain(domain=Domain([], metas=[meta_var]), n_rows=5)
        table.get_column_view(meta_var)[0][:] = meta

        self.send_signal(self.widget.Inputs.data, table)
        self.widget.commit(force=True)
        self.show()

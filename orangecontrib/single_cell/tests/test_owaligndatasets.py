import unittest

from sklearn import datasets
from Orange.data import Table, Domain, ContinuousVariable, DiscreteVariable
import numpy as np

from Orange.data import Table, Domain
from Orange.widgets.data.owtranspose import OWTranspose
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.widgets.owaligndatasets import OWAlignDatasets


class TestOWMDA(WidgetTest):

    def generate_dataset(self, n=3, cols=10, samples=1000):
        cluster_std = 2.5
        f = lambda x, i: .5 * x * i + 2 * i

        X, Y = datasets.make_blobs(n_samples=samples, n_features=cols, centers=1,
                                   cluster_std=cluster_std, random_state=10)
        for i in range(samples):
            y = i % n
            Y[i] = y
            if y > 0:
                X[i, cols // 2:] = f(X[i, cols // 2:], i)
        dom = Domain([ContinuousVariable.make(str(i)) for i in range(cols)],
                     DiscreteVariable.make('class', values=[str(i) for i in range(n)]))
        return Table.from_numpy(dom, X, Y[np.newaxis].T)

    def setUp(self):
        self.widget = self.create_widget(OWAlignDatasets)

    def test_for_nans(self):
        """ Check if output transformed data contains any NaNs. """
        data = self.generate_dataset()
        self.send_signal(self.widget.Inputs.data, data)

        gc = self.get_output(self.widget.Outputs.genes_components)
        td = self.get_output(self.widget.Outputs.transformed_data)

        self.assertIsNotNone(gc)
        self.assertIsNotNone(td)

        self.assertFalse(np.isnan(gc.X).any())

    def test_shape(self):
        """ Check shape of outputs. """
        n = 3
        cols = 10
        samples = 1000
        data = self.generate_dataset(n=n, cols=cols, samples=samples)
        self.send_signal(self.widget.Inputs.data, data)

        ngenes = cols
        ncomponents = cols // 2
        self.widget.ngenes = ngenes
        self.widget.ncomponents = ncomponents
        self.widget.fit()

        gc = self.get_output(self.widget.Outputs.genes_components)
        td = self.get_output(self.widget.Outputs.transformed_data)

        self.assertEqual(td.X.shape, (samples, ncomponents))
        self.assertEqual(gc.X.shape, (ngenes, ncomponents))

from sklearn import datasets
from Orange.data import ContinuousVariable, DiscreteVariable
import numpy as np

from Orange.data import Table, Domain
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.tests.utils import simulate

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

    def test_source_id_defult(self):
        """
        Check whether Source ID is set as a default source ID attribute.
        """
        w = self.widget
        data = self.generate_dataset()
        data1 = Table("iris")
        data = Table(
            Domain(
                data.domain.attributes, data.domain.class_var,
                metas=[DiscreteVariable("Source ID", values=["a", "b"])]),
            data.X, data.Y, metas=np.random.randint(0, 2, (len(data), 1)))
        data1 = Table(
            Domain(
                data1.domain.attributes, data1.domain.class_var,
                metas=[DiscreteVariable("Source ID", values=["a", "b", ])]),
            data1.X, data1.Y, metas=np.random.randint(0, 2, (len(data1), 1)))
        w._reset_settings()
        self.send_signal(w.Inputs.data, data)
        # without manual fix source ID would be class since it is first in list
        self.assertEqual("Source ID", str(w.source_id))
        # we change setting manually (it is stored in widget settings) and it should stay even when new signal comes
        cbox = w.controls.source_id
        simulate.combobox_activate_item(cbox, "class")
        self.assertEqual("class", str(w.source_id))
        self.send_signal(w.Inputs.data, None)
        self.send_signal(w.Inputs.data, data)
        self.assertEqual("class", str(w.source_id))

        self.send_signal(w.Inputs.data, data1)
        self.assertEqual("Source ID", str(w.source_id))

        self.send_signal(w.Inputs.data, data)
        self.assertEqual("class", str(w.source_id))

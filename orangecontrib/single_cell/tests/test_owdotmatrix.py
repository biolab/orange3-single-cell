from unittest.mock import patch, Mock

import numpy as np

from Orange.data import Table, Domain, DiscreteVariable
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.tests.utils import simulate

from orangecontrib.single_cell.widgets.owdotmatrix import OWDotMatrix


class TestOWDotMatrix(WidgetTest):

    def setUp(self) -> None:
        self.iris = Table("iris")
        self.widget = self.create_widget(OWDotMatrix)

    def test_set_data(self):
        w = self.widget

        # send iris - dataset where analysis should pass
        self.send_signal(w.Inputs.data, self.iris)

        self.assertIsNotNone(w.data)
        self.assertIsNotNone(w.cluster_var)
        self.assertFalse(w.Error.no_discrete_variable.is_shown())

        # send housing - dataset where analysis should fail - no discrete var
        self.send_signal(w.Inputs.data, Table("housing"))

        self.assertIsNone(w.data)
        self.assertIsNone(w.cluster_var)
        self.assertTrue(w.Error.no_discrete_variable.is_shown())

        # reset input
        self.send_signal(w.Inputs.data, None)

        self.assertIsNone(w.data)
        self.assertIsNone(w.cluster_var)
        self.assertFalse(w.Error.no_discrete_variable.is_shown())

    def test_data_aggregation(self):
        w = self.widget

        # send iris - dataset where analysis should pass
        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(3, len(w.clusters_unordered))
        self.assertTupleEqual((3, 4), w.aggregated_data.shape)

        self.send_signal(w.Inputs.data, None)
        self.assertIsNone(w.clusters_unordered)
        self.assertIsNone(w.aggregated_data)

    def test_data_aggregation_change_attr(self):
        """
        With this test we check whether discrete attribute for aggregation changes correctly
        """
        w = self.widget

        # add additional discrete attribute to data
        iris = Table(
            Domain(
                self.iris.domain.attributes,
                self.iris.domain.class_var,
                metas=[DiscreteVariable("a", values=["a", "b"])]),
            self.iris.X,
            self.iris.Y,
            metas=[[0]] * 100 + [[1]] * 50
        )

        # send iris - dataset where analysis should pass
        self.send_signal(w.Inputs.data, iris)
        cbox = self.widget.controls.cluster_var
        simulate.combobox_activate_index(cbox, 0)

        self.assertEqual(3, len(w.clusters_unordered))
        self.assertTupleEqual((3, 4), w.aggregated_data.shape)

        cbox = self.widget.controls.cluster_var
        simulate.combobox_activate_index(cbox, 2)
        self.assertEqual(2, len(w.clusters_unordered))
        self.assertTupleEqual((2, 4), w.aggregated_data.shape)

    def test_info_string(self):
        w = self.widget
        insum = w.info.set_input_summary = Mock()
        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(insum.call_args[0][0], "150")
        self.assertEqual(insum.call_args[0][1], "4 genes\n150 cells\n3 clusters")

        self.send_signal(w.Inputs.data, None)
        self.assertEqual(insum.call_args[0][0], w.info.NoInput)

    def test_calculate_table_values_no_norm(self):
        w = self.widget
        w.biclustering = False
        w.transpose = False
        w.log_scale = False
        w.normalize = False

        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(Table, type(w.matrix))
        np.testing.assert_almost_equal(w.matrix.X, w._norm_min_max(w.matrix_before_norm.X))
        np.testing.assert_almost_equal(w.matrix.X, w._norm_min_max(w.aggregated_data))
        np.testing.assert_almost_equal(w.matrix_before_norm.X, w.aggregated_data)
        self.assertListEqual(sorted(w.clusters), sorted(self.iris.domain.class_var.values))

    def test_calculate_table_values_norm(self):
        w = self.widget
        w.biclustering = False
        w.transpose = False
        w.log_scale = False
        w.normalize = True

        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(Table, type(w.matrix))
        self.assertTupleEqual(w.matrix.X.shape, w.matrix_before_norm.X.shape)
        self.assertTupleEqual(w.matrix.X.shape, w.aggregated_data.shape)
        self.assertListEqual(sorted(w.clusters), sorted(self.iris.domain.class_var.values))

    def test_calculate_table_values_log(self):
        w = self.widget
        w.biclustering = False
        w.transpose = False
        w.log_scale = True
        w.normalize = False

        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(Table, type(w.matrix))
        self.assertTupleEqual(w.matrix.X.shape, w.matrix_before_norm.X.shape)
        self.assertTupleEqual(w.matrix.X.shape, w.aggregated_data.shape)
        self.assertListEqual(sorted(w.clusters), sorted(self.iris.domain.class_var.values))

    def test_calculate_table_values_clustering(self):
        w = self.widget
        w.biclustering = True
        w.transpose = False
        w.log_scale = False
        w.normalize = False

        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(Table, type(w.matrix))
        self.assertTupleEqual(w.matrix.X.shape, w.matrix_before_norm.X.shape)
        self.assertTupleEqual(w.matrix.X.shape, w.aggregated_data.shape)
        self.assertListEqual(sorted(w.clusters), sorted(self.iris.domain.class_var.values))
        self.assertListEqual(sorted(map(str, w.matrix.domain.attributes)),
                             sorted(map(str, self.iris.domain.attributes)))

    def test_calculate_table_values_transpose(self):
        w = self.widget
        w.biclustering = False
        w.transpose = True
        w.log_scale = False
        w.normalize = False

        self.send_signal(w.Inputs.data, self.iris)
        self.assertEqual(Table, type(w.matrix))
        self.assertTupleEqual(w.matrix.X.shape, w.matrix_before_norm.X.shape)
        self.assertTupleEqual(w.matrix.X.shape, w.aggregated_data.T.shape)
        self.assertListEqual(sorted(w.clusters), sorted(map(str, self.iris.domain.attributes)))
        self.assertListEqual(sorted(map(str, w.matrix.domain.attributes)),
                             sorted(self.iris.domain.class_var.values))

    def test_output(self):
        w = self.widget

        self.send_signal(w.Inputs.data, self.iris)
        ann_data = self.get_output(w.Outputs.annotated_data)
        sel_data = self.get_output(w.Outputs.selected_data)
        cont_data = self.get_output(w.Outputs.contingency)

        self.assertIsNotNone(ann_data)
        self.assertIsNone(sel_data)
        self.assertIsNotNone(cont_data)

        np.testing.assert_almost_equal(ann_data.X, w.data.X)
        np.testing.assert_almost_equal(ann_data.Y, w.data.Y)
        self.assertEqual(1, len(ann_data.domain.metas))
        np.testing.assert_almost_equal(cont_data.X, w.matrix.X)
        self.assertEqual(1, len(cont_data.domain.metas))
        self.assertEqual(1, len(cont_data.domain.metas))
        self.assertEqual("iris", str(cont_data.domain.metas[0]))

        self.send_signal(w.Inputs.data, None)

        ann_data = self.get_output(w.Outputs.annotated_data)
        sel_data = self.get_output(w.Outputs.selected_data)
        cont_data = self.get_output(w.Outputs.contingency)

        self.assertIsNone(ann_data)
        self.assertIsNone(sel_data)
        self.assertIsNone(cont_data)

    def test_output_transpose(self):
        w = self.widget

        self.send_signal(w.Inputs.data, self.iris)
        w.controls.transpose.click()

        ann_data = self.get_output(w.Outputs.annotated_data)
        sel_data = self.get_output(w.Outputs.selected_data)
        cont_data = self.get_output(w.Outputs.contingency)

        self.assertIsNotNone(ann_data)
        self.assertIsNone(sel_data)
        self.assertIsNotNone(cont_data)

        np.testing.assert_almost_equal(ann_data.X, w.data.X)
        np.testing.assert_almost_equal(ann_data.Y, w.data.Y)
        self.assertEqual(1, len(ann_data.domain.metas))
        np.testing.assert_almost_equal(cont_data.X, w.matrix.X)
        self.assertEqual(1, len(cont_data.domain.metas))
        self.assertEqual(1, len(cont_data.domain.metas))
        self.assertEqual("Gene", str(cont_data.domain.metas[0]))

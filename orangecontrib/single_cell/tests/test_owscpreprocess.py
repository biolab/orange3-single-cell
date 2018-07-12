import numpy as np

from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.single_cell.preprocess.scpreprocess import (
    LogarithmicScale, Binarize, Normalize,
    NormalizeGroups, NormalizeSamples, Standardize, SelectMostVariableGenes)
from orangecontrib.single_cell.widgets import owscpreprocess
from orangecontrib.single_cell.widgets.owscpreprocess import OWscPreprocess


class TestOWscPreprocess(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWscPreprocess)

    def test_available_preprocessors(self):
        self.assertEqual(self.widget.preprocessors.rowCount(), 5)

    def test_default_preprocessors(self):
        self.assertEqual(self.widget.preprocessormodel.rowCount(), 4)

    def test_preprocess(self):
        data = Table("iris")
        self.send_signal(self.widget.Inputs.data, data)
        pp_data = self.get_output(self.widget.Outputs.preprocessed_data)
        self.assertIsInstance(pp_data, Table)
        self.assertNotEqual(pp_data, data)

    def test_discrete_attributes(self):
        self.send_signal(self.widget.Inputs.data, Table("iris"))
        pp_data = self.get_output(self.widget.Outputs.preprocessed_data)
        self.assertIsNotNone(pp_data)
        self.send_signal("Data", Table("heart_disease"))
        pp_data = self.get_output(self.widget.Outputs.preprocessed_data)
        self.assertIsNone(pp_data)
        self.assertTrue(self.widget.Error.discrete_attributes.is_shown())
        self.send_signal("Data", None)
        self.assertFalse(self.widget.Error.discrete_attributes.is_shown())

    def test_missing_values(self):
        data = Table("iris")
        data[0, 3] = np.nan
        self.send_signal(self.widget.Inputs.data, data)
        pp_data = self.get_output(self.widget.Outputs.preprocessed_data)
        self.assertIsInstance(pp_data, Table)
        self.assertNotEqual(pp_data, data)
        self.assertTrue(self.widget.Warning.missing_values.is_shown())
        self.send_signal("Data", None)
        self.assertFalse(self.widget.Warning.missing_values.is_shown())
        pp_data = self.get_output(self.widget.Outputs.preprocessed_data)
        self.assertIsNone(pp_data)


class TestLogarithmicScaleEditor(WidgetTest):
    def test_editor_default(self):
        widget = owscpreprocess.LogarithmicScaleEditor()
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, LogarithmicScale)
        self.assertEqual(p.base, LogarithmicScale.BinaryLog)


class TestBinarizeEditor(WidgetTest):
    def test_editor_default(self):
        widget = owscpreprocess.BinarizeEditor()
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, Binarize)
        self.assertEqual(p.condition, Binarize.GreaterOrEqual)
        self.assertEqual(p.threshold, 1)


class TestNormalizeEditor(WidgetTest):
    def test_editor_default(self):
        pp = owscpreprocess.OWscPreprocess()
        widget = owscpreprocess.NormalizeEditor(master=pp)
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, NormalizeSamples)
        self.assertEqual(p.method, Normalize.CPM)

    def test_editor_normalize_groups(self):
        data = Table("iris")
        pp = owscpreprocess.OWscPreprocess()
        pp.set_data(data)
        widget = owscpreprocess.NormalizeEditor(master=pp)
        params = {"method": Normalize.Median,
                  "group_var": data.domain.class_var,
                  "group_by": True}
        widget.setParameters(params)
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, NormalizeGroups)
        self.assertEqual(p.group_var, data.domain.class_var)
        self.assertEqual(p.method, Normalize.Median)


class TestStandardizeEditor(WidgetTest):
    def test_editor_default(self):
        widget = owscpreprocess.StandardizeEditor()
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, Standardize)
        self.assertIsNone(p.lower_bound)
        self.assertIsNone(p.upper_bound)


class TestSelectGenesEditor(WidgetTest):
    def test_editor_default(self):
        widget = owscpreprocess.SelectGenesEditor()
        p = widget.createinstance(widget.parameters())
        self.assertIsInstance(p, SelectMostVariableGenes)
        self.assertEqual(p.method, SelectMostVariableGenes.Dispersion)
        self.assertEqual(p.n_genes, 1000)
        self.assertEqual(p.n_groups, 20)

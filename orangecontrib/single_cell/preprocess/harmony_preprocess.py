from Orange.data import Table, Domain
from Orange.preprocess.preprocess import Preprocess
from .harmony_full_numpy import harmony_numpy_fun
from Orange.data import Domain, ContinuousVariable

class HarmonyModel:
    def __init__(self, batch_vars, harmony_params):
        """
        batch_vars      = lista di Orange.Variable
        harmony_params  = parametri da passare a harmony_numpy_fun
        """
        self.batch_vars = batch_vars
        self.harmony_params = harmony_params

    def transform(self, data):
        X_raw = data.X
        batch_list = [data.get_column(v) for v in self.batch_vars]

        # PCA + Harmony
        Z_corr = harmony_numpy_fun(
            X_raw=X_raw,
            batches=batch_list,
            n_pcs=50,     
            **self.harmony_params
        )

        k = Z_corr.shape[1]     
        attrs = [ContinuousVariable(f"PC{i+1}") for i in range(k)]

        class_vars = list(data.domain.class_vars)
        metas = list(data.domain.metas)

        new_domain = Domain(attrs, class_vars, metas)
        new_Y = data.Y
        new_metas = data.metas
        new_data = Table(new_domain, Z_corr, new_Y, new_metas)

        return new_data

    def __call__(self, data):
        return self.transform(data)

class HarmonyNormalizer(Preprocess):
    def __init__(self, batch_vars=(), **harmony_params):
        """
        batch_vars      = lista di nomi delle variabili meta (stringhe)
        harmony_params  = parametri di harmony_numpy_fun
        """
        self.batch_vars = batch_vars
        self.harmony_params = harmony_params

    def __call__(self, data):

        vars_obj = [data.domain[b] for b in self.batch_vars]

        model = HarmonyModel(vars_obj, self.harmony_params)

        return model.transform(data)

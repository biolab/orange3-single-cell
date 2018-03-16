import numpy as np
from Orange.data import Domain, Table, ContinuousVariable, DiscreteVariable
from Orange.data.util import SharedComputeValue
from Orange.preprocess.preprocess import Preprocess, Continuize


# Link / inverse link functions
LINK_IDENTITY = "Identity link"
LINK_LOG = "Log link"

LINKS = {
    LINK_IDENTITY: lambda x: x,
    LINK_LOG: lambda x: np.log(x)
}

INV_LINKS = {
    LINK_IDENTITY: lambda x: x,
    LINK_LOG: lambda x: np.exp(x)
}


class ScBatchShared(SharedComputeValue):
    """Places the values of shared data within the corresponding variable column."""
    def compute(self, data, shared_data):
        assert self.variable is not None
        return shared_data.get_column_view(self.variable)[0] if self.variable in shared_data.domain else np.nan


class SCBatchNormalizer(Preprocess):
    """ Instantiate a new domain with transformations defined by the model. """
    def __init__(self, link=LINK_IDENTITY, nonzero_only=True, batch_vars=()):
        self.link = link
        self.nonzero_only = nonzero_only
        self.batch_vars = batch_vars

    def __call__(self, data):
        proj = ScBatchNormalizeModel(self.link, self.nonzero_only, self.batch_vars)
        proj.fit(data)
        attributes = [var.copy(compute_value=ScBatchShared(proj, variable=var))
                      for var in data.domain.attributes]
        for var in attributes:
            var.number_of_decimals = max(3, var.number_of_decimals)
        normalized_domain = Domain(
            attributes, data.domain.class_vars, data.domain.metas)
        return data.transform(normalized_domain)


class ScBatchNormalizeModel:

    def __init__(self, link=LINK_IDENTITY, nonzero_only=True, batch_vars=()):
        """
        :param link: Link function key.
        :param nonzero_only: Fit only on non-zero values.
        :param batch_vars: *names* of batch variables (must be meta).
        """
        if link == LINK_LOG and not nonzero_only:
            raise ValueError("Log link must be used with nonzero_only=True !")
        self.link = link
        self.nonzero_only = nonzero_only
        self.batch_vars = batch_vars
        self.models = dict()

    def _design_matrix(self, data):
        """ Create a design matrix with Continuized variables and a bias term.
            Can contain both meta and class variables. """
        assert len(self.batch_vars) > 0
        df = data[:, self.batch_vars]
        dom = df.domain.metas + df.domain.class_vars
        X = np.hstack((df.metas,
                       df.Y.reshape((len(data), len(df.domain.class_vars)))))
        Z = Continuize()(Table.from_numpy(domain=Domain(attributes=dom), X=X)).X
        return np.hstack((np.ones((len(df), 1)), Z))

    def fit(self, data):
        """ Fit one model per gene with least-squares. """
        atts = data.domain.attributes
        assert all([isinstance(a, ContinuousVariable) for a in atts])
        assert all([isinstance(data.domain[b], ContinuousVariable)
                    or isinstance(data.domain[b], DiscreteVariable) for b in self.batch_vars])
        if len(self.batch_vars) == 0:
            return
        Z = self._design_matrix(data)
        k = Z.shape[1]
        if self.nonzero_only:
            for i, a in enumerate(atts):
                nz = np.where(data.X[:, i])[0]
                Y = LINKS[self.link](data.X[:, i][nz])
                w = np.linalg.lstsq(Z[nz], Y, rcond=None)[0].reshape((k, 1)) \
                    if len(nz) else np.zeros((k, 1))
                self.models[a.name] = w
        else:
            Y = LINKS[self.link](data.X)
            if not (np.all(np.isfinite(Y)) and np.all(np.logical_not(np.isnan(Y)))):
                raise ValueError("Transformed data contains NaN/inf values. "
                                 "Use a different link function.")
            W = np.linalg.lstsq(Z, Y, rcond=None)[0]
            self.models = dict(((a.name, w.reshape((k, 1))) for a, w in zip(atts, W.T)))

    def transform(self, data):
        """ Apply transformation. Genes in new data must have had been available to fit. """
        if len(self.batch_vars) == 0:
            return Table.from_table(domain=data.domain,
                                    source=data)
        else:
            atts = data.domain.attributes
            assert all((a.name in self.models for a in atts))
            assert (data.domain.index(b) for b in self.batch_vars)
            Z = self._design_matrix(data)
            W = np.hstack((self.models[a.name].reshape((Z.shape[1], 1)) for a in atts))
            Xc = data.X.copy()
            if self.nonzero_only:
                nz = np.where(Xc)
                Xn = INV_LINKS[self.link](LINKS[self.link](Xc[nz]) - Z.dot(W)[nz])
                Xc[nz] = Xn

            else:
                Xc = INV_LINKS[self.link](LINKS[self.link](Xc) - Z.dot(W))
            return Table.from_numpy(domain=data.domain,
                                    metas=data.metas,
                                    X=Xc,
                                    Y=data.Y,
                                    W=data.W)

    def __call__(self, data):
        """ Transform data (alias). """
        return self.transform(data)

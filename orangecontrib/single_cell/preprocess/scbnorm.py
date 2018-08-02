import numpy as np
from scipy.special import betainc as betai
from sklearn.preprocessing import OneHotEncoder

from Orange.data import Domain, Table, ContinuousVariable, DiscreteVariable, Variable
from Orange.data.util import SharedComputeValue
from Orange.preprocess.preprocess import Preprocess, Continuize
from Orange.preprocess.score import Scorer

# Link / inverse link functions
LINK_IDENTITY = "Identity link"
LINK_LOG = "Log link"

LINKS = {
    LINK_IDENTITY: lambda x: x,
    LINK_LOG: np.log
}

INV_LINKS = {
    LINK_IDENTITY: lambda x: x,
    LINK_LOG: np.exp
}

# Percentage standard deviation
PERC_SD = 0.01


class ScBatchScorer(Scorer):
    """
    For Continuous batch variables,
        calculate a percentage of genes significantly (p < alpha)
        correlated with the variable.

    For Discrete batch variables,
        calculate the *relative size of union* of
        *significantly correlated genes* (p < alpha) for each value (one-hot encoding).
    """

    feature_type = Variable
    class_type = Variable
    supports_sparse_data = True
    friendly_name = "score"
    name = "ScBatchScore"

    @staticmethod
    def _std(A):
        """ Numpy implementation fails on matrices with one column. """
        return np.power(np.var(A, axis=0), 0.5)

    @staticmethod
    def correlations(A, B):
        """ Fast approximation to Pearson corr. and p-values for variables in two matrices."""
        # Correlation coefficients
        M = (A - A.mean(axis=0)).T.dot(B - B.mean(axis=0)) / A.shape[0]
        S = np.outer(ScBatchScorer._std(A), ScBatchScorer._std(B))
        S[S == 0] = 1.0
        M[S == 0] = 0.0
        rf = M / S

        # P-values
        df = A.shape[0] - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = betai(0.5 * df, 0.5, np.array(df / (df + ts), dtype=float))
        return rf, pf

    def __init__(self, alpha=0.05):
        self.alpha = alpha

    def score_data(self, data, feature):
        if feature is None:
            raise ValueError("Scorer %s computes on a per-feature basis. ", self.__class__)
        if not all((isinstance(att, ContinuousVariable) for att in data.domain.attributes)):
            raise ValueError("All variables in the data must be Continuous!")

        a = data.get_column_view(feature.name)[0].reshape((len(data), 1))
        if isinstance(feature, ContinuousVariable):
            _, p = self.correlations(a, data.X)
            w = (p < self.alpha).mean()
        else:
            B = OneHotEncoder(sparse=False).fit_transform(a)
            _, P = self.correlations(B, data.X)
            w = (P < self.alpha).sum(axis=0).astype(bool).mean()
        return w

    def __call__(self, data, feature=None):
        return self.score_data(data, feature)



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
            if np.any(data.X < 0):
                if self.link == LINK_IDENTITY:
                    raise ValueError("Data contains negative values. "
                                     "Remove 'nonzero_true' option.")
                else:
                    raise ValueError("Data contains negative values.")
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
        """ Apply transformation. Genes in new data must have had been available to fit.

        In case of linear regression (identity link), with nonzero_only=True,
        make sure zeros remains zeros and the original non-zeros are all larger than zero.

        This combination must not be used if the original data contains negative
        values.
        """
        if len(self.batch_vars) == 0:
            return data
        else:
            atts = data.domain.attributes
            assert all((a.name in self.models for a in atts))
            assert all(b in data.domain for b in self.batch_vars)
            Z = self._design_matrix(data)
            W = np.hstack((self.models[a.name].reshape((Z.shape[1], 1)) for a in atts))
            Xc = data.X.copy()
            if self.nonzero_only:
                nz = np.where(Xc)
                Xc[nz] = INV_LINKS[self.link](LINKS[self.link](Xc[nz]) - Z.dot(W)[nz])
                if self.link == LINK_IDENTITY:
                    mask = Xc != 0
                    sd = Xc.std() * PERC_SD
                    Xc = Xc - (Xc.min(axis=0) * mask) + (sd * mask)
            else:
                Xc = INV_LINKS[self.link](LINKS[self.link](Xc) - Z.dot(W))
            new_data = data.copy()
            new_data.X = Xc
            return new_data

    def __call__(self, data):
        """ Transform data (alias). """
        return self.transform(data)

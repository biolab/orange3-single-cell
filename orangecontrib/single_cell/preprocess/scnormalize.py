import numpy as np
import scipy.sparse as sp

from Orange.data.table import Table
from Orange.data.domain import Domain
from Orange.data.variable import ContinuousVariable
from Orange.projection import Projection, Projector
from Orange.data.util import SharedComputeValue
from Orange.preprocess.preprocess import Preprocess


class ScShared(SharedComputeValue):
    """Places the values of shared data within the coresponding variable column."""
    def compute(self, data, shared_data):
        assert self.variable is not None
        return np.array(shared_data[:, self.variable]).ravel()


class ScNormalizeProjection(Projection, Preprocess):
    """
    Defines a transformation of the domain.
    Initialized with an instance of ScNormalizeModel, which has __call__(self, data) method,
    and is used to compute shared_data in ScShared.
    """
    name = 'ScNormalize'
    def __init__(self, proj, domain):
        super().__init__(proj=proj)
        self.orig_domain = domain
        self.domain = Domain(
            [ContinuousVariable(name=var.name,
                                compute_value=ScShared(proj, variable=var))
             for var in self.orig_domain.attributes],
            domain.class_vars, domain.metas)


class ScNormalizeProjector(Projector):
    """
    Fits data to a model instance and returns a Projection instance based on current model instances
    Projector has a default __call__(self, data) method that calls self.fit(data.X, data.Y) and
    sets domain.
    """
    name = 'ScNormalize'
    supports_sparse = True

    def __init__(self,
                 domain,
                 equalize_var=None,
                 normalize_cells=True,
                 log_base=2):
        super().__init__()
        self.equalize_var = equalize_var
        self.normalize_cells = normalize_cells
        self.log_base = log_base
        self.domain = domain


    def fit(self, X, Y=None):
        proj = ScNormalizeModel(self.equalize_var,
                                self.normalize_cells,
                                self.log_base)
        proj.fit(X, Y)
        return ScNormalizeProjection(proj, self.domain)



class ScNormalizeModel:
    """
    Implements model fitting and transformation.
    Parameters infered from data are stored inside an instance of this class.
     A simple ad-hoc normalization to provide basic raw count pre-processing.
    """

    def __init__(self, equalize_var=None, normalize_cells=True, log_base=2):
        """
        :param equalize_var: Equalization variable.
        :param normalize_cells: Normalize cell profiles.
        :param log_base: Base for log-trasnform. Use None to skip.
        """
        self.equalize_var = equalize_var
        self.normalize_cells = normalize_cells
        self.log_base = log_base
        self.target_row_mean = 1
        self.size_factors = {}

    def fit(self, X, Y=None):
        """
        Infer row normalization parameters from the data.
        :param X: Continuous data matrix.
        :param Y: Grouping values
        :return:
        """
        # Equalize based on read depth per library / match mean read count per cell
        # Must not store indices
        if Y is not None:
            libraries = dict([(lib, np.where(Y == lib)[0]) for lib in set(Y)])
            lib_sizes = dict()
            for lib, inxs in sorted(libraries.items()):
                lib_sizes[lib] = np.median(X[inxs, :].sum(axis=1))
            self.target_row_mean = min(lib_sizes.values())
            for lib in libraries.keys():
                self.size_factors[lib] = self.target_row_mean / lib_sizes[lib]
        else:
            self.target_row_mean = np.median(np.array(X.sum(axis=1)))

    def __call__(self, data):
        """
        :param data: Data to be transformed.
        :return:
        """
        return self.transform(data)

    def transform(self, data):
        """
        Transform data based on inferred parameters.
        :param data: Data table with expression values as counts.
                    Columns are genes and rows are cells.
        :return: Data table with normalized values.
        """
        # Result in expected number of reads
        Xeq = data.X.copy()
        n = Xeq.shape[0]

        # Normalize cell profiles

        if self.normalize_cells:
            # Each cell is normalized independently by default
            rs = np.array(Xeq.sum(axis=1))
            rsm = np.ones((n, )) * self.target_row_mean
            factors = rs / rsm

            # Override with library size factor, if provided. Else, each row is
            # treated as a separate group
            if self.equalize_var is not None:
                vals = np.array(list(map(lambda lib: self.size_factors.get(lib, np.nan),
                                    data.get_column_view(self.equalize_var)[0])))
                inxs = np.logical_not(np.isnan(vals))
                factors[inxs] = vals[inxs]

            Xd = sp.dia_matrix((factors.ravel(), 0), shape=(n, n))
            Xeq = Xd.dot(Xeq)

        # Log transform log(1 + x)
        if self.log_base is not None:
            if sp.isspmatrix(Xeq):
                Xeq = Xeq.log1p() / np.log(self.log_base)
            else:
                Xeq = np.log(1 + Xeq) / np.log(self.log_base)

        # Preserve sparsity
        X_new = Xeq.tocsr() if sp.isspmatrix(Xeq) else Xeq
        data_new = Table.from_numpy(domain=data.domain,
                                    X=X_new,
                                    Y=data.Y,
                                    W=data.W,
                                    metas=data.metas)
        return data_new

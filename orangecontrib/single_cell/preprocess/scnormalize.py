import numpy as np
import scipy.sparse as sp

from Orange.data.table import Table
from Orange.preprocess.preprocess import Preprocess

class ScNormalize(Preprocess):
    """
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


    def __call__(self, data):
        """
        :param data: Data table with expression values as counts.
                    Columns are genes and rows are cells.
        :return: Data table with normalized values.
        """
        # Result in expected number of reads
        Xeq = data.X.copy()

        # Equalize based on read depth per library / match mean read count per cell
        if self.equalize_var is not None:
            y, _  = data.get_column_view(self.equalize_var)
            lib_sizes = dict()
            libraries = dict([(lib, np.where(y == lib)[0]) for lib in set(y)])
            for lib, inxs in sorted(libraries.items()):
                lib_sizes[lib] = data.X[inxs, :].sum(axis=1).mean()
            target = min(lib_sizes.values())
            size_factors = dict([(lib, target / float(size)) for lib, size in lib_sizes.items()])
            for lib, inxs in sorted(libraries.items()):
                Xeq[inxs, :] *= size_factors[lib]

        # Normalize by cells, sweep columns by means / median
        if self.normalize_cells:
            rs = np.array(Xeq.sum(axis=1).reshape((Xeq.shape[0], 1)))
            rsm = np.median(rs)
            Xd = sp.dia_matrix(((rsm / rs).ravel(), 0), shape=(len(rs), len(rs)))
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
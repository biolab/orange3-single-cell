import os
import random
from itertools import chain
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import loompy as lp
import xlrd

from Orange.data import (
    ContinuousVariable, DiscreteVariable, Domain, Table
)
from Orange.data.io import (
    Compression, open_compressed, PickleReader,
    guess_data_type, sanitize_variable
)
from Orange.widgets.utils.filedialogs import RecentPath


def separator_from_filename(file_name):
    """Get separator from file extension

    :param file_name: str
    :return: str
    """
    return "," if os.path.splitext(file_name)[1] == ".csv" else "\t"


def get_data_loader(file_name):
    """Get instance of data loader according to file extension

    :param file_name: str
    :return: Loader
    """
    base, ext = os.path.splitext(file_name)
    if ext in Compression.all:
        _, ext = os.path.splitext(base)
    if ext == ".mtx":
        return MtxLoader(file_name)
    elif ext == ".count":
        return CountLoader(file_name)
    elif ext == ".csv":
        return CsvLoader(file_name)
    elif ext in (".pkl", ".pickle"):
        return PickleLoader(file_name)
    elif ext in (".xls", ".xlsx"):
        return ExcelLoader(file_name)
    elif ext == ".loom":
        return LoomLoader(file_name)
    elif ext == ".h5ad":
        return H5ADLoader(file_name)
    else:
        return Loader(file_name)


class Loader:
    """Class loads sample file (and it's supplementary annotation
    files, if present) and creates Orange.data.Table out of its data
    """
    separator = "\t"

    def __init__(self, file_name=""):
        # file parameters
        self._file_name = file_name
        self._recent_path = RecentPath.create(file_name, [])
        self.file_size = None
        self.n_rows = None
        self.n_cols = None
        self.sparsity = None
        self._set_file_size()
        self._set_file_parameters()

        # reading parameters
        self._leading_cols = 0
        self._leading_rows = 0
        self._use_rows_mask = None
        self._use_cols_mask = None
        self._row_annot_header = 0
        self._row_annot_columns = None
        self._col_annot_header = 0
        self._col_annot_columns = None

        # GUI parameters
        self.FIXED_FORMAT = True
        self.ENABLE_ANNOTATIONS = True
        self.header_rows_count = None
        self.header_cols_count = None
        self.transposed = None
        self.sample_rows_enabled = None
        self.sample_cols_enabled = None
        self.sample_rows_p = None
        self.sample_cols_p = None
        self.row_annotations_enabled = True
        self.row_annotation_file = None  # type: RecentPath
        self.col_annotations_enabled = True
        self.col_annotation_file = None  # type: RecentPath

        # errors
        self.errors = {}  # type: Dict[str, Tuple]
        self.__reset_error_messages()

    @property
    def recent_path(self):
        return self._recent_path

    @recent_path.setter
    def recent_path(self, path):
        self._recent_path = path
        self._file_name = path.abspath

    @property
    def n_genes(self):
        return self.n_rows if self.transposed else self.n_cols

    @property
    def n_cells(self):
        return self.n_cols if self.transposed else self.n_rows

    @property
    def leading_rows(self):
        return self._leading_rows

    @leading_rows.setter
    def leading_rows(self, value):
        self._leading_rows = value

    @property
    def leading_cols(self):
        return self._leading_cols

    @leading_cols.setter
    def leading_cols(self, value):
        self._leading_cols = value

    def _set_file_size(self):
        try:
            st = os.stat(self._file_name)
        except OSError:
            pass
        else:
            self.file_size = st.st_size

    def _set_file_parameters(self):
        try:
            args = (self._file_name,)
            col_kwargs = {"sep": self.separator, "index_col": None, "nrows": 1}
            row_kwargs = {"sep": self.separator, "usecols": [1]}
            self.n_cols = self.df_read_func(*args, **col_kwargs).shape[1]
            self.n_rows = self.df_read_func(*args, **row_kwargs).shape[0]
            self._set_sparsity()
        except Exception:
            pass

    def _set_sparsity(self):
        """Get approximate sparsity if number of columns is bigger than 100."""
        if self.n_cols is not None and self.n_rows is not None:
            n_max = 100
            np.random.seed(42)
            use_cols = np.arange(1, self.n_cols - 1) if self.n_cols < n_max \
                else np.random.randint(1, self.n_cols - 1, n_max)

            data = self.df_read_func(
                self._file_name, usecols=use_cols, sep=self.separator,
                skipfooter=(max(0, self.n_rows - n_max))
            )

            n_all = data.shape[0] * data.shape[1]
            self.sparsity = (n_all - np.count_nonzero(data)) / n_all

    def _load_data(self, skip_row=None, header_rows=None, header_cols=None,
                   use_cols=None, **kwargs):
        df = self.df_read_func(
            self._file_name, sep=self.separator, index_col=header_cols,
            header=header_rows, skiprows=skip_row, usecols=use_cols
        )
        if skip_row is not None:
            self._use_rows_mask = np.array(self._use_rows_mask, dtype=bool)
        if self.transposed:
            df = df.transpose()
            self._use_rows_mask, self._use_cols_mask = \
                self._use_cols_mask, self._use_rows_mask
            self.leading_rows, self.leading_cols = \
                self.leading_cols, self.leading_rows

        attrs = [ContinuousVariable.make(str(g)) for g in df.columns]

        return attrs, df.values, df.iloc[:, :0], df.index

    @staticmethod
    def df_read_func(*args, **kwargs):
        return pd.read_csv(*args, **kwargs)

    def __call__(self):
        self.__reset_error_messages()

        rst = np.random.RandomState(0x667)
        if self.transposed:
            skip_row, skip_col = self.__skip_col(rst), self.__skip_row(rst)
        else:
            skip_col, skip_row = self.__skip_col(rst), self.__skip_row(rst)

        header_rows = self.__header_rows()
        header_cols, header_cols_indices = self.__header_cols()

        usecols, skip_col, skip_row, header_cols, header_rows = \
            self.__update_reading_parameters(
                skip_col, skip_row, header_cols_indices,
                header_cols, header_rows
            )

        try:
            attrs, X, meta_df, meta_df_index = self._load_data(
                skip_row=skip_row, skip_col=skip_col,
                header_rows=header_rows, header_cols=header_cols,
                use_cols=usecols, transpose=self.transposed
            )
        except Exception as e:
            self.errors["reading_error"] = (e, None)
            return None

        meta_parts = (meta_df,)
        if self.row_annotations_enabled and self.row_annotation_file:
            meta_df, row_annot_df = self.__update_metas(
                meta_df, meta_df_index, X)
            if row_annot_df is not None:
                meta_parts = (meta_df, row_annot_df)

        if self.col_annotations_enabled and self.col_annotation_file:
            attrs = self.__update_attributes(attrs, X)

        return self.__into_orange_table(attrs, X, meta_parts)

    def __header_rows(self):
        header_rows = self.header_rows_count
        header_rows_indices = []
        if header_rows == 0:
            header_rows = None
        elif header_rows == 1:
            header_rows = 0
            header_rows_indices = [0]
        else:
            header_rows = list(range(header_rows))
            header_rows_indices = header_rows
        self.leading_rows = len(header_rows_indices)
        return header_rows

    def __header_cols(self):
        header_cols = self.header_cols_count
        header_cols_indices = []
        if header_cols == 0:
            header_cols = None
        elif header_cols == 1:
            header_cols = 0
            header_cols_indices = [0]
        else:
            header_cols = list(range(header_cols))
            header_cols_indices = header_cols
        self.leading_cols = len(header_cols_indices)
        return header_cols, header_cols_indices

    def __skip_row(self, rstate):
        skip_row = None
        if self.sample_rows_enabled:
            p = self.sample_rows_p
            if p < 100:
                def skip_row(i, p=p):
                    return i > 3 and rstate.uniform(0, 100) > p
        return skip_row

    def __skip_col(self, rstate):
        skip_col = None
        if self.sample_cols_enabled:
            p = self.sample_cols_p
            if p < 100:
                def skip_col(i, p=p):
                    return i > 3 and rstate.uniform(0, 100) > p
        return skip_col

    def __update_reading_parameters(self, skip_col, skip_row, header_indices,
                                    header_cols, header_rows):
        self._use_rows_mask = None
        self._use_cols_mask = None
        usecols = None

        if skip_col is not None:
            self._use_cols_mask = np.array([
                not skip_col(i) or i in header_indices
                for i in range(self.n_cols)], dtype=bool)
            usecols = np.flatnonzero(self._use_cols_mask)

        if skip_row is not None:
            self._use_rows_mask = []  # record the used rows

            def skip_row(i, test=skip_row):
                r = test(i)
                self._use_rows_mask.append(r)
                return r

        # if skipping more than one leading columns/rows, use the
        # first one when creating domain
        if isinstance(header_cols, list) and self.transposed:
            if skip_col is None:
                mask = np.array([i not in header_cols[1:] for i
                                 in range(self.n_cols)], dtype=bool)
            else:
                mask = np.array([i not in header_cols[1:] and not skip_col(i)
                                 for i in range(self.n_cols)], dtype=bool)
            self._use_cols_mask = mask
            usecols = np.flatnonzero(mask)
            header_cols = 0

        if isinstance(header_rows, list) and not self.transposed:
            if skip_row is None:
                skip_row = header_rows[1:]
            else:
                self._use_rows_mask = []

                def skip_row(i, test=skip_row, header=header_rows[1:]):
                    r = test(i) or i in header
                    self._use_rows_mask.append(r)
                    return r

            header_rows = 0

        return usecols, skip_col, skip_row, header_cols, header_rows

    def __update_metas(self, meta_df, meta_df_index, X):
        file_name = self.row_annotation_file.abspath
        row_annot_df = pd.read_csv(
            file_name, sep=separator_from_filename(file_name),
            header=self._row_annot_header, names=self._row_annot_columns
        )
        if self._use_rows_mask is not None:
            # NOTE: we account for column header/ row index
            expected = len(self._use_rows_mask) - self.leading_rows
        else:
            expected = X.shape[0]
        if len(row_annot_df) != expected:
            self.errors["row_annot_mismatch"] = (expected, len(row_annot_df))
            row_annot_df = None

        if row_annot_df is not None and self._use_rows_mask is not None:
            # use the same sample indices
            indices = np.flatnonzero(self._use_rows_mask[self.leading_rows:])
            row_annot_df = row_annot_df.iloc[indices]
            # if path.endswith(".count") and row_annot.endswith('.meta'):
            #     assert np.all(row_annot_df.iloc[:, 0] == df.index)

        if row_annot_df is not None and meta_df_index is not None:
            # Try to match the leading columns with the meta_df_index.
            # If found then drop the columns (or index if the level does
            # not have a name but the annotation col does)
            drop_cols = []
            drop_index_level = []
            for i in range(meta_df_index.nlevels):
                meta_df_level = meta_df_index.get_level_values(i)
                if np.all(row_annot_df.iloc[:, i] == meta_df_level):
                    if meta_df_level.name is None:
                        drop_index_level.append(i)
                    elif meta_df_level.name == row_annot_df.columns[i].name:
                        drop_cols.append(i)

            if drop_cols:
                row_annot_df = row_annot_df.drop(columns=drop_cols)

            if drop_index_level:
                for i in reversed(drop_index_level):
                    if isinstance(meta_df.index, pd.MultiIndex):
                        meta_df_index = meta_df_index.droplevel(i)
                    else:
                        assert i == 0
                        meta_df_index = pd.RangeIndex(meta_df_index.size)
                meta_df = pd.DataFrame({}, index=meta_df_index)
        return meta_df, row_annot_df

    def __update_attributes(self, attrs, X):
        file_name = self.col_annotation_file.abspath
        col_annot_df = pd.read_csv(
            file_name, sep=separator_from_filename(file_name),
            header=self._col_annot_header, names=self._col_annot_columns
        )
        if self._use_cols_mask is not None:
            expected = len(self._use_cols_mask) - self.leading_cols
        else:
            expected = X.shape[1]
        if len(col_annot_df) != expected:
            self.errors["col_annot_mismatch"] = (expected, len(col_annot_df))
            col_annot_df = None
        if col_annot_df is not None and self._use_cols_mask is not None:
            indices = np.flatnonzero(self._use_cols_mask[self.leading_cols:])
            col_annot_df = col_annot_df.iloc[indices]

        if col_annot_df is not None:
            assert len(col_annot_df) == X.shape[1]
            if not attrs and X.shape[1]:  # No column names yet
                attrs = [ContinuousVariable.make(str(v))
                         for v in col_annot_df.iloc[:, 0]]
            names = [str(c) for c in col_annot_df.columns]
            for var, values in zip(attrs, col_annot_df.values):
                var.attributes.update(
                    {n: v for n, v in zip(names, values)})
        return attrs

    def __into_orange_table(self, attrs, X, meta_parts):
        if not attrs and X.shape[1]:
            attrs = Domain.from_numpy(X).attributes

        try:
            metas = None
            M = None
            if meta_parts:
                meta_parts = [df_.reset_index() if not df_.index.is_integer()
                              else df_ for df_ in meta_parts]
                metas, M = self.__guess_metas(meta_parts)

            domain = Domain(attrs, metas=metas)
            table = Table.from_numpy(domain, X, None, M)
        except ValueError:
            table = None
            rows = self.leading_cols if self.transposed else self.leading_rows
            cols = self.leading_rows if self.transposed else self.leading_cols
            self.errors["inadequate_headers"] = (rows, cols)
        return table

    @staticmethod
    def __guess_metas(meta_parts):
        def guessed_var(i, var_name, dtype):
            if pd.core.dtypes.common.is_numeric_dtype(dtype):
                return ContinuousVariable.make(var_name)
            orig_values = M[:, i]
            val_map, values, var_type = guess_data_type(orig_values)
            values, variable = sanitize_variable(
                val_map, values, orig_values, var_type, {}, name=var_name)
            M[:, i] = values
            return variable

        M = np.hstack(tuple(df_.values for df_ in meta_parts))
        return [guessed_var(i, name, dtype) for i, (name, dtype) in
                enumerate(zip(chain(*(_.columns for _ in meta_parts)),
                              chain(*(_.dtypes for _ in meta_parts))))], M

    def __reset_error_messages(self):
        self.errors = {"row_annot_mismatch": (),
                       "col_annot_mismatch": (),
                       "inadequate_headers": (),
                       "reading_error": ()}

    def copy(self):
        loader = self.__class__(self._file_name)
        for key in vars(loader):
            setattr(loader, key, getattr(self, key))
        return loader


class MtxLoader(Loader):
    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 0
        self.header_cols_count = 0
        self.FIXED_FORMAT = False
        self.transposed = True
        self._row_annot_header = None
        self._row_annot_columns = ["Barcodes"]
        self._col_annot_header = None
        self._col_annot_columns = ["Id", "Gene"]
        self._set_annotation_files()
        self._set_enable_annotations()

    @property
    def leading_rows(self):
        return self._leading_rows

    @leading_rows.setter
    def leading_rows(self, value):
        self._leading_rows = 0

    @property
    def leading_cols(self):
        return self._leading_cols

    @leading_cols.setter
    def leading_cols(self, value):
        self._leading_cols = 0

    def _set_annotation_files(self):
        dir_name, _ = os.path.split(self._file_name)
        genes_path = os.path.join(dir_name, "genes.tsv")
        if os.path.isfile(genes_path):
            self.col_annotation_file = RecentPath.create(genes_path, [])
        barcodes_path = os.path.join(dir_name, "barcodes.tsv")
        if os.path.isfile(barcodes_path):
            self.row_annotation_file = RecentPath.create(barcodes_path, [])

    def _set_enable_annotations(self):
        # 10x gene-barcode matrix
        # TODO: The genes/barcodes files should be unconditionally loaded
        # alongside the mtx. The row/col annotations might be used to
        # specify additional sources. For the time being they are put in
        # the corresponding comboboxes and made uneditable.
        if self.row_annotation_file is not None \
                and self.col_annotation_file is not None \
                and os.path.basename(
                self.col_annotation_file.abspath) == "genes.tsv" \
                and os.path.basename(
                self.row_annotation_file.abspath) == "barcodes.tsv":
            self.ENABLE_ANNOTATIONS = False

    def _set_file_parameters(self):
        try:
            with open_compressed(self._file_name, "rb") as f:
                self.n_rows, self.n_cols, non_zero_el = scipy.io.mminfo(f)[:3]
                all_el = self.n_rows * self.n_cols
                self.sparsity = (all_el - non_zero_el) / all_el
        except OSError:
            pass
        except ValueError:
            pass

    def _load_data(self, skip_row=None, skip_col=None, **kwargs):
        X = scipy.io.mmread(self._file_name)
        if self.transposed:
            X = X.T
        if skip_row is not None:
            self._use_rows_mask = np.array(
                [not skip_row(i) for i in range(X.shape[0])]
            )
            X = X.tocsr()[np.flatnonzero(self._use_rows_mask)]
        if skip_col is not None:
            self._use_cols_mask = np.array(
                [not skip_col(i) for i in range(X.shape[1])]
            )
            X = X.tocsc()[:, np.flatnonzero(self._use_cols_mask)]
        X = X.todense(order="F")
        if self._use_rows_mask is not None:
            meta_df = pd.DataFrame(
                {}, index=np.flatnonzero(self._use_rows_mask))
        else:
            meta_df = pd.DataFrame({}, index=pd.RangeIndex(X.shape[0]))
        return [], X, meta_df, meta_df.index


class CountLoader(Loader):
    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 1
        self.header_cols_count = 1
        self.FIXED_FORMAT = False
        self.transposed = True
        self._set_annotation_files()

    def _set_annotation_files(self):
        dir_name, basename = os.path.split(self._file_name)
        basename_no_ext, _ = os.path.splitext(basename)
        meta_path = os.path.join(dir_name, basename_no_ext + ".meta")
        if os.path.isfile(meta_path):
            self.row_annotation_file = RecentPath.create(meta_path, [])


class CsvLoader(Loader):
    separator = ","


class PickleLoader(Loader):
    separator = None

    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 0
        self.header_cols_count = 0
        self.FIXED_FORMAT = False
        self.ENABLE_ANNOTATIONS = False
        self.transposed = False
        self.row_annotations_enabled = False
        self.col_annotations_enabled = False

    def _set_file_parameters(self):
        pass

    def _load_data(self):
        random.seed(0)
        reader = PickleReader(self._file_name)
        table = reader.read()
        self.n_rows, self.n_cols = table.X.shape
        attrs = self.__attributes(table.domain.attributes)
        metas = table.domain.metas + table.domain.class_vars
        domain = Domain(attrs, metas=metas)
        return table.transform(domain)[self.__row_indices()]

    def __call__(self):
        return self._load_data()

    def __row_indices(self):
        indices = range(self.n_rows)
        if self.sample_rows_enabled:
            p = self.sample_rows_p
            if p < 100 and self.n_rows > 3:
                return random.sample(indices, int(self.n_rows * p / 100))
        return indices

    def __attributes(self, attributes):
        if self.sample_cols_enabled:
            p = self.sample_cols_p
            if p < 100 and self.n_cols > 3:
                return random.sample(attributes, int(self.n_cols * p / 100))
        return attributes


class ExcelLoader(Loader):
    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 1
        self.header_cols_count = 1

    def _set_file_parameters(self):
        try:
            sheet = xlrd.open_workbook(self._file_name).sheet_by_index(0)
            self.n_cols = sheet.ncols
            self.n_rows = sheet.nrows
            self._set_sparsity()
        except Exception:
            pass

    @staticmethod
    def df_read_func(*args, **kwargs):
        return pd.read_excel(*args, **kwargs)


class LoomLoader(Loader):
    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 0
        self.header_cols_count = 0
        self.FIXED_FORMAT = False
        self.ENABLE_ANNOTATIONS = False
        self.transposed = True
        self.row_annotations_enabled = False
        self.col_annotations_enabled = False

    def _set_file_parameters(self):
        try:
            with lp.connect(self._file_name) as ds:
                self.n_rows, self.n_cols = ds.shape
                all_el = self.n_rows * self.n_cols
                self.sparsity = (all_el - ds.sparse().count_nonzero()) / all_el
        except OSError:
            pass

    def _load_data(self, skip_row=None, skip_col=None, **kwargs):
        with lp.connect(self._file_name) as ds:
            X = ds[:, :].T
            if skip_row is not None:
                mask = np.array([not skip_row(i) for i in range(X.shape[1])])
                self._use_rows_mask = mask
            else:
                self._use_rows_mask = np.ones(X.shape[1], dtype=bool)
            if skip_col is not None:
                mask = np.array([not skip_col(i) for i in range(X.shape[0])])
                self._use_cols_mask = mask
            else:
                self._use_cols_mask = np.ones(X.shape[0], dtype=bool)

            X = X[self._use_cols_mask, :]
            X = X[:, self._use_rows_mask]
            gene_names = ds.ra.Gene[self._use_rows_mask] \
                if hasattr(ds.ra, "Gene") else []
            attrs = [ContinuousVariable.make(str(g)) for g in gene_names]
            meta_df = pd.DataFrame({key: ds.ca[key][self._use_cols_mask]
                                    for key in ds.ca.keys()})

        return attrs, X, meta_df, meta_df.index


class H5ADLoader(Loader):
    def __init__(self, file_name):
        super().__init__(file_name)
        self.header_rows_count = 0
        self.header_cols_count = 0
        self.FIXED_FORMAT = False
        self.ENABLE_ANNOTATIONS = False
        self.transposed = False
        self.row_annotations_enabled = False
        self.col_annotations_enabled = False

    def _set_file_parameters(self):
        import anndata
        try:
            adata = anndata.read_h5ad(self._file_name)
            self.n_rows, self.n_cols = adata.shape
            all_el = self.n_rows * self.n_cols
            if sp.issparse(adata.X):
                self.sparsity = (all_el - adata.X.tocsr().count_nonzero()) / all_el
            else:
                self.sparsity = (all_el - np.count_nonzero(adata.X)) / all_el
        except OSError:
            pass

    def _load_data(self, skip_row=None, skip_col=None, **kwargs):
        import anndata
        adata = anndata.read_h5ad(self._file_name)

        if skip_row is not None:
            mask = np.array([not skip_row(i) for i in range(adata.shape[0])])
            self._use_rows_mask = mask
        else:
            self._use_rows_mask = np.ones(adata.shape[0], dtype=bool)
        if skip_col is not None:
            mask = np.array([not skip_col(i) for i in range(adata.shape[1])])
            self._use_cols_mask = mask
        else:
            self._use_cols_mask = np.ones(adata.shape[1], dtype=bool)

        adata = adata[self._use_rows_mask, :]
        adata = adata[:, self._use_cols_mask]

        attrs = [ContinuousVariable.make(str(g)) for g in adata.var_names]
        meta_df = pd.DataFrame({key: adata.obs[key] for key in adata.obs.keys()})

        X = adata.X.toarray() if sp.issparse(adata.X) else adata.X

        return attrs, X, meta_df, meta_df.index


class Concatenate:
    INTERSECTION, UNION = range(2)

    @classmethod
    def concatenate(cls, concat_type, data_collection):
        if not data_collection:
            return None

        concat_data, source_var = cls.append_source_name(*data_collection[0])
        for data, source_name in data_collection[1:]:
            attrs1 = concat_data.domain.attributes
            attrs2 = data.domain.attributes
            if concat_type == cls.INTERSECTION:
                attrs = set(attrs1).intersection(attrs2)
            elif concat_type == cls.UNION:
                attrs = set(attrs1 + attrs2)
            metas = set(concat_data.domain.metas + data.domain.metas)

            def key(var):
                return var.name if isinstance(var.name, str) else ""

            domain = Domain(sorted(attrs, key=key),
                            metas=sorted(metas, key=key))
            concat_data_t = concat_data.transform(domain)
            data_t = data.transform(domain)
            source_var.values + (source_name, )
            data_t[:, source_var] = np.full(
                (len(data), 1), len(source_var.values) - 1, dtype=object
            )
            concat_data = Table.concatenate((concat_data_t, data_t), axis=0)
        return concat_data

    @staticmethod
    def append_source_name(data, name):
        source_var = DiscreteVariable("source", values=[name])
        metas = data.domain.metas + (source_var,)
        domain = Domain(data.domain.attributes, metas=metas)
        data = data.transform(domain)
        data[:, source_var] = np.full((len(data), 1), 0, dtype=object)
        return data, source_var

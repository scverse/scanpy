"""
Annotated Data
"""
from collections import Mapping, Sequence
from collections import OrderedDict
from enum import Enum
import numpy as np
from numpy import ma
from numpy.lib.recfunctions import append_fields
from scipy import sparse as sp
from scipy.sparse.sputils import IndexMixin
from ..utils import merge_dicts


SMP_NAMES = 'smp_names'
VAR_NAMES = 'var_names'


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sp.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


class BoundRecArr(np.recarray):
    """
    A np.recarray that can be constructed from a dict.

    Is bound to AnnData to allow adding fields.

    The column "smp_name"/"var_name" plays a role analogous to the `index` in
    the `pd.Series` class.

    Attributes
    ----------
    _index_name : str
        Either SMP_NAMES or VAR_NAMES.
    _parent : AnnData
        The reference to an AnnData object to which the array is bound.
    """

    def __new__(cls, source, index_name, parent, n_row=None):
        if source is None:  # empty array
            cols = [np.arange(n_row).astype(str)]
            dtype = [(index_name, cols[0].dtype)]
        elif isinstance(source, np.recarray):
            cols = [source[n] for n in source.dtype.names]
            dtype = source.dtype
        else:
            if not isinstance(source, Mapping):
                raise ValueError(
                    'meta needs to be a recarray or dictlike, not {}'
                    .format(type(source)))
            # meta is dict-like
            names = list(source.keys())
            cols = [np.asarray(col) for col in source.values()]
            if index_name not in source:
                names.append(index_name)
                cols.append(np.arange(len(cols[0]) if cols else n_row).astype(str))
            dtype = list(zip(names, [str(c.dtype) for c in cols]))
        try:
            dtype = np.dtype(dtype)
        except TypeError:
            # TODO: fix compat with Python 2
            # print(dtype, file=sys.stderr)
            raise

        arr = np.recarray.__new__(cls, (len(cols[0]),), dtype)
        arr._parent = parent  # add parent as attribute
        arr._index_name = index_name  # add index_name as attribute

        for i, name in enumerate(dtype.names):
            arr[name] = np.array(cols[i], dtype=dtype[name])

        return arr

    def flipped(self):
        old_index_name = self._index_name
        new_index_name = SMP_NAMES if old_index_name == VAR_NAMES else VAR_NAMES

        flipped = BoundRecArr(self, new_index_name, self._parent, len(self))
        flipped.dtype.names = tuple(
            new_index_name if n == old_index_name else n
            for n in self.dtype.names)

        return flipped

    def copy(self):
        new = super(BoundRecArr, self).copy()
        new._index_name = self._index_name
        new._parent = self._parent
        return new

    @property
    def columns(self):
        return [c for c in self.dtype.names if not c == self._index_name]

    def _keys_view(self, arr, keys):
        """Get a multi-column view rather than a copy.

        This is consistent with the behavior of single columns.

        Will be implemented that way from numpy 1.13 on anyways, so we do not use
        in __getitem__, instead suppress the warning raised there.

        Note
        ----
        http://stackoverflow.com/questions/15182381/how-to-return-a-view-of-several-columns-in-numpy-structured-array
        """
        dtype2 = np.dtype({name:arr.dtype.fields[name] for name in keys})
        return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

    def __getitem__(self, k):
        """Either a single one- or multi-column or mulitiple one-colum items."""
        import warnings  # ignore FutureWarning about multi-column access of structured arrays
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                return super(BoundRecArr, self).__getitem__(k)
            except:  # access multiple columns with single key k
                keys = []
                i = 0
                for name in self.dtype.names:
                    if name.startswith(k) and 'of' in name:
                        keys.append(k + str(i) + 'of' + name.split('of')[-1])
                        i += 1
                        if i == int(name.split('of')[-1]):
                            break
                if i == 0:
                    raise ValueError('Could not find item corresponding to key' + k)
                return super(BoundRecArr, self).__getitem__(keys)

    def __setitem__(self, keys, values):
        """Either a single one- or multi-column or mulitiple one-colum items."""
        if isinstance(keys, str):
            keys = [keys]
            if (not hasattr(values[0], '__len__')  # check if values is nested
                or len(values[0]) == 1 or isinstance(values[0], str)):
                values = [values]
        keys = np.array(keys)
        values = np.array(values)  # sequence of arrays or matrix with n_keys *rows*
        if not (len(keys) == len(values) or len(keys) == 1):
            raise ValueError('You passed {} column keys but {} arrays as columns. '
                             'If you passed a matrix instead of a sequence '
                             'of arrays, try transposing it.'
                             .format(len(keys), len(values)))

        # if we have several values but a single key, generate unique keys
        if len(values) > len(keys):
            keys = [keys[0] + str(i) + 'of' + str(len(values))
                    for i in range(len(values))]

        present = np.intersect1d(keys, self.dtype.names)
        absent = np.setdiff1d(keys, self.dtype.names)

        if any(present):
            for k, v in zip(present, values[np.in1d(keys, present)]):
                super(BoundRecArr, self).__setitem__(k, v)

        if any(absent):
            attr = 'smp' if self._index_name == SMP_NAMES else 'var'
            if values.shape[1] > len(self):
                raise ValueError('New column has too many entries ({} > {})'
                                 .format(values.shape[1], len(self)))
            source = append_fields(self, absent, values[np.in1d(keys, absent)],
                                   usemask=False, asrecarray=True)
            new = BoundRecArr(source, self._index_name, self._parent)
            setattr(self._parent, attr, new)

    def to_df(self):
        """Return pd.dataframe."""
        import pandas as pd
        return pd.DataFrame().from_records(self, index=self._index_name)


def _check_dimensions(data, smp, var):
    n_smp, n_var = data.shape
    if len(smp) != n_smp:
        raise ValueError('Sample metadata needs to have the same amount of '
                         'rows as data has ({}), but has {} rows'
                         .format(n_smp, smp.shape[0]))
    if len(var) != n_var:
        raise ValueError('Feature metadata needs to have the same amount of '
                         'rows as data has columns ({}), but has {} rows'
                         .format(n_var, var.shape[0]))


class AnnData(IndexMixin):

    def __init__(self, data=None, smp=None, var=None, add=None):
        """
        Annotated Data

        Stores a data matrix X of dimensions n_samples x n_variables,
        e.g. n_cells x n_genes, with the possibility to store an arbitrary
        number of annotations for both samples and variables, and
        additional arbitrary unstructured annotation via add.

        You can access additional annotation elements directly from AnnData:
        >>> adata = AnnData(np.eye(3), k=1)
        >>> assert adata['k'] == 1

        Parameters
        ----------
        data : dict, np.ndarray, np.ma.MaskedArray, sp.spmatrix
            The data matrix `X`
            X : np.ndarray, np.ma.MaskedArray, sp.spmatrix
                A n_samples x n_variables data matrix.
            or a dict containing `X` as 'X' and possibly
            'row_names' / 'smp_names' : list, np.ndarray, optional
                A n_samples array storing names for samples.
            'col_names' / 'var_names' : list, np.ndarray, optional
                A n_variables array storing names for variables.
            'row' / 'smp' : dict, optional
                A dict with row annotation.
            'col' / 'var' : dict, optional
                A dict with column annotation.
        smp : np.recarray, dict
            A n_samples x ? record array containing sample names (`smp_names`)
            and other sample annotation in the columns. A passed dict is
            converted to a record array.
        var : np.recarray, dict
            The same as `smp`, but of shape n_variables x ? for annotation of
            variables.
        add : dict
            Unstructured annotation for the whole dataset.

        Attributes
        ----------
        X, smp, var from the Parameters.
        """
        if isinstance(data, Mapping):
            if any((smp, var, add)):
                raise ValueError('If `data` is a dict no further arguments must be provided.')
            X, smp, var, add = self.from_ddata(data)
        else:
            X = data

        # check data type of X
        for s_type in StorageType:
            if isinstance(X, s_type.value):
                self.storage_type = s_type
                break
        else:
            class_names = ', '.join(c.__name__ for c in StorageType.classes())
            raise ValueError('X needs to be of one of the following types [{}] not {}'
                             .format(class_names, type(X)))

        # use lower precision, is enough for all current applications
        if X.dtype == np.float64:
            X = X.astype(np.float32)

        if len(X.shape) == 1:
            X.shape = (X.shape[0], 1)
        if X.dtype.names is None and len(X.shape) != 2:
            raise ValueError('X needs to be 2-dimensional, not '
                             '{}D'.format(len(X.shape)))

        n_smp, n_var = X.shape

        self.X = X

        self.smp = BoundRecArr(smp, SMP_NAMES, self, n_smp)
        self.var = BoundRecArr(var, VAR_NAMES, self, n_var)

        _check_dimensions(X, self.smp, self.var)

        self.add = add

    def __repr__(self):
        return ('AnnData object with attributes\n'
                '    X : array-like data matrix of shape n_samples x n_variables = '
                + str(self.X.shape[0]) + ' x ' + str(self.X.shape[1]) + '\n'
                '    smp : dict-like annotation for each sample (e.g., the sample names)\n'
                '    var : dict-like annotation for each variable (e.g., the variable names)\n'
                '    add : dict-like additional unstructured annotation')

    def from_ddata(self, ddata):
        smp, var = OrderedDict(), OrderedDict()

        add = dict(ddata.items())
        del ddata

        X = add['X']
        del add['X']

        if 'row_names' in add:
            smp['smp_names'] = add['row_names']
            del add['row_names']
        elif 'smp_names' in add:
            smp['smp_names'] = add['smp_names']
            del add['smp_names']

        if 'col_names' in add:
            var['var_names'] = add['col_names']
            del add['col_names']
        elif 'var_names' in add:
            var['var_names'] = add['var_names']
            del add['var_names']

        smp = merge_dicts(smp, add.get('row', {}), add.get('smp', {}))
        var = merge_dicts(var, add.get('col', {}), add.get('var', {}))
        for k in ['row', 'smp', 'col', 'var']:
            if k in add:
                del add[k]

        return X, smp, var, add

    def to_ddata(self):
        smp = OrderedDict([(k, self.smp[k]) for k in self.smp_keys()])
        var = OrderedDict([(k, self.var[k]) for k in self.var_keys()])
        d = {'X': self.X, 'smp': smp, 'var': var,
             'smp_names': self.smp_names, 'var_names': self.var_names}
        for k, v in self.add.items():
            d[k] = v
        return d

    def smp_keys(self):
        """Return keys of sample annotation, excluding `smp_names`."""
        return [n for n in self.smp.dtype.names if n != SMP_NAMES]

    def var_keys(self):
        """Return keys of variable annotation, excluding `var_names`."""
        return [n for n in self.var.dtype.names if n != VAR_NAMES]

    def add_keys():
        """Return keys of addtional unstructured annotation."""
        return self.add.keys()

    @property
    def smp_names(self):
        return self.smp[SMP_NAMES]

    @smp_names.setter
    def smp_names(self, keys):
        self.smp[SMP_NAMES] = keys

    @property
    def var_names(self):
        return self.var[VAR_NAMES]

    @var_names.setter
    def var_names(self, keys):
        self.var[VAR_NAMES] = keys

    def __setattr__(self, key, value):
        names_col = dict(smp=SMP_NAMES, var=VAR_NAMES).get(key)
        if names_col and not isinstance(value, BoundRecArr):  # if smp/var is set, give it the right class
            names_orig, dim = (self.smp_names, 0) if names_col == SMP_NAMES else (self.var_names, 1)
            value_orig, value = value, BoundRecArr(value, names_col, self)
            if len(value) != self.X.shape[dim]:
                raise ValueError('New value for {!r} was converted to a reacarray of length {} instead of {}'
                                 .format(key, len(value_orig), len(self)))
            if (value[names_col] == np.arange(self.X.shape[dim])).all():  # TODO: add to constructor
                value[names_col] = names_orig
        object.__setattr__(self, key, value)

    def _normalize_indices(self, packed_index):
        smp, var = super(AnnData, self)._unpack_index(packed_index)
        smp = self._normalize_index(smp, self.smp_names)
        var = self._normalize_index(var, self.var_names)
        return smp, var

    def _normalize_index(self, index, names):
        def name_idx(i):
            if isinstance(i, str):
                # `where` returns an 1-tuple (1D array) of found indices
                i = np.where(names == i)[0][0]
                if i is None:
                    raise IndexError('Index {} not in smp_names/var_names'
                                     .format(index))
            return i

        if isinstance(index, slice):
            start = name_idx(index.start)
            stop = name_idx(index.stop)
            # string slices can only be inclusive, so +1 in that case
            if isinstance(index.stop, str):
                stop = None if stop is None else stop + 1
            step = index.step
        elif isinstance(index, (int, str)):
            start = name_idx(index)
            stop = start + 1
            step = 1
        elif isinstance(index, (Sequence, np.ndarray)):
            return np.fromiter(map(name_idx, index), 'int64')
        else:
            raise IndexError('Unknown index {!r} of type {}'
                             .format(index, type(index)))

        return slice(start, stop, step)

    def __delitem__(self, index):
        smp, var = self._normalize_indices(index)
        del self.X[smp, var]
        if var == slice(None):
            del self.smp.iloc[smp, :]
        if smp == slice(None):
            del self.var.iloc[var, :]

    def __getitem__(self, index):
        # otherwise unpack index
        smp, var = self._normalize_indices(index)
        X = self.X[smp, var]
        smp_ann = self.smp[smp]
        var_ann = self.var[var]
        assert smp_ann.shape[0] == X.shape[0], (smp, smp_ann)
        assert var_ann.shape[0] == X.shape[1], (var, var_ann)
        adata = AnnData(X, smp_ann, var_ann,  self.add)
        return adata

    def __setitem__(self, index, val):
        if isinstance(index, str):
            self.add[index] = val
            return

        smp, var = self._normalize_indices(index)
        self.X[smp, var] = val

    def __contains__(self, item):
        return item in self.add

    def get(self, key, default=None):
        return self.add.get(key, default)

    def __len__(self):
        return self.X.shape[0]

    def transpose(self):
        """Return a transposed view of the object.

        Sample axis (rows) and variable axis are interchanged. No additional memory.
        """
        if sp.issparse and isinstance(self.X, sp.csr_matrix):
            return AnnData(self.X.T.tocsr(), self.var.flipped(), self.smp.flipped(), self.add)
        return AnnData(self.X.T, self.var.flipped(), self.smp.flipped(), self.add)

    T = property(transpose)

    def copy(self):
        """Full copy with memory allocated."""
        return AnnData(self.X.copy(), self.smp.copy(), self.var.copy(), self.add.copy())


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(ma.array([[1, 2], [3, 4]], add={'mask': [0, 1, 1, 0]})
    AnnData(sp.eye(2))
    AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(Smp=['A', 'B']),
        dict(Feat=['a', 'b', 'c']))

    assert AnnData(np.array([1, 2])).X.shape == (2, 1)

    from pytest import raises
    raises(ValueError, AnnData,
           np.array([[1, 2], [3, 4]]),
           dict(TooLong=[1, 2, 3, 4]))


def test_ddata():
    ddata = dict(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        row_names=['A', 'B'],
        col_names=['a', 'b', 'c'])
    AnnData(ddata)


def test_names():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata.smp_names.tolist() == 'A B'.split()
    assert adata.var_names.tolist() == 'a b c'.split()


def test_get_subset():
    mat = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    assert mat[0, 0].X.tolist() == [[1]]
    assert mat[0, :].X.tolist() == [[1, 2, 3]]
    assert mat[:, 0].X.tolist() == [[1], [4]]
    assert mat[:, [0, 1]].X.tolist() == [[1, 2], [4, 5]]
    assert mat[:, np.array([0, 2])].X.tolist() == [[1, 3], [4, 6]]
    assert mat[:, np.array([False, True, True])].X.tolist() == [[2, 3], [5, 6]]
    assert mat[:, 1:3].X.tolist() == [[2, 3], [5, 6]]


def test_get_subset_names():
    mat = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert mat['A', 'a'].X.tolist() == [[1]]
    assert mat['A', :].X.tolist() == [[1, 2, 3]]
    assert mat[:, 'a'].X.tolist() == [[1], [4]]
    assert mat[:, ['a', 'b']].X.tolist() == [[1, 2], [4, 5]]
    assert mat[:, np.array(['a', 'c'])].X.tolist() == [[1, 3], [4, 6]]
    assert mat[:, 'b':'c'].X.tolist() == [[2, 3], [5, 6]]

    from pytest import raises
    with raises(IndexError): _ = mat[:, 'X']
    with raises(IndexError): _ = mat['X', :]
    with raises(IndexError): _ = mat['A':'X', :]
    with raises(IndexError): _ = mat[:, 'a':'X']


def test_transpose():
    mat = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    mt1 = mat.T

    # make sure to not modify the original!
    assert mat.smp_names.tolist() == ['A', 'B']
    assert mat.var_names.tolist() == ['a', 'b', 'c']

    assert SMP_NAMES in mt1.smp.dtype.names
    assert VAR_NAMES in mt1.var.dtype.names
    assert mt1.smp_names.tolist() == ['a', 'b', 'c']
    assert mt1.var_names.tolist() == ['A', 'B']
    assert mt1.X.shape == mat.X.T.shape

    mt2 = mat.transpose()
    assert np.array_equal(mt1.X, mt2.X)
    assert np.array_equal(mt1.smp, mt2.smp)
    assert np.array_equal(mt1.var, mt2.var)


def test_get_subset_add():
    mat = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                  dict(Smp=['A', 'B']),
                  dict(Feat=['a', 'b', 'c']))

    assert mat[0, 0].smp['Smp'].tolist() == ['A']
    assert mat[0, 0].var['Feat'].tolist() == ['a']


def test_append_add_col():
    mat = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    mat.smp['new'] = [1, 2]
    mat.smp[['new2', 'new3']] = [['A', 'B'], ['c', 'd']]

    from pytest import raises
    with raises(ValueError):
        mat.smp['new4'] = 'far too long'.split()


def test_set_add():
    mat = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    mat.smp = dict(smp_names=[1, 2])
    assert isinstance(mat.smp, BoundRecArr)
    assert len(mat.smp.dtype) == 1

    mat.smp = dict(a=[1, 2])  # leave smp_names and a custom column
    assert isinstance(mat.smp, BoundRecArr)
    assert len(mat.smp.dtype) == 2
    assert mat.smp_names.tolist() == [1, 2]

    from pytest import raises
    with raises(ValueError):
        mat.smp = dict(a=[1, 2, 3])


def test_multi_column_getitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp[['a', 'b']] = np.array([[0, 1], [2, 3]])
    print(adata.smp[['a', 'b']])


def test_multi_column_single_key_getitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp[['c0of2', 'c1of2']] = np.array([[0, 1], [2, 3]])
    print(adata.smp)
    print(adata.smp['c'])


def test_multi_column_single_key_setitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp['c'] = np.array([[0, 1], [2, 3]])
    print(adata.smp)

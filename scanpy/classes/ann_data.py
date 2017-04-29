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


SMP_INDEX = 'smp_names'
VAR_INDEX = 'var_names'


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sp.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


def _key_belongs_to_which_key_multicol(key, keys_multicol):
    for imk, mk in enumerate(keys_multicol):
        if key.startswith(mk) and 'of' in key:
            return imk
    return -1


def _gen_keys_from_key_multicol(key_multicol, n_keys):
    """Generates single-column keys from multicolumn key."""
    keys = [key_multicol + str(i+1) + 'of' + str(n_keys) for i in range(n_keys)]
    return keys


class BoundRecArr(np.recarray):
    """
    A np.recarray that can be constructed from a dict.

    Is bound to AnnData to allow adding fields.

    The column "smp_name"/ "var_name" plays a role analogous to the `index` in
    the `pd.Series` class.

    Attributes
    ----------
    _index_name : str
        Either SMP_INDEX or VAR_INDEX.
    _parent : AnnData
        The reference to an AnnData object to which the array is bound.
    """

    def __new__(cls, source, index_name, parent, n_row=None, keys_multicol=None):
        if source is None:  # empty array
            cols = [np.arange(n_row).astype(str)]
            dtype = [(index_name, cols[0].dtype)]
        elif isinstance(source, np.recarray):
            cols = [source[n] for n in source.dtype.names]
            dtype = source.dtype
        else:
            if not isinstance(source, Mapping):
                raise ValueError('Expected recarray or dictlike type, not {}.'
                                 .format(type(source)))
            # meta is dict-like
            names = list(source.keys())
            cols = [np.asarray(col) for col in source.values()]
            if index_name not in source:
                names.append(index_name)
                cols.append(np.arange(len(cols[0]) if cols else n_row).astype(str))
            else:
                if not isinstance(source[index_name][0], str):
                    raise ValueError('Value for key "{}" for initializing index name '
                                     'needs to be of type str, not {}.'
                                     .format(index_name, type(source[index_name][0])))
            dtype = list(zip(names, [str(c.dtype) for c in cols]))
        try:
            dtype = np.dtype(dtype)
        except TypeError:
            # TODO: fix compat with Python 2
            # print(dtype, file=sys.stderr)
            raise

        arr = np.recarray.__new__(cls, (len(cols[0]),), dtype)
        arr._parent = parent  # add parent as attribute
        arr._index_name = index_name  # add index_name as attribue

        # only the multicol keys
        arr._keys_multicol = ([] if keys_multicol is None
                              else list(keys_multicol))
        # fast lookup of single-column keys
        arr._keys_multicol_lookup = ({} if keys_multicol is None
                                     else dict.fromkeys(keys_multicol))
        # all keys (multi- and single-col keys that do not belong to multicol)
        arr._keys = []  # this excludes arr._index_name
        loop_over_keys = (key for key in arr.dtype.names if key != index_name)
        for key in loop_over_keys:
            imk = _key_belongs_to_which_key_multicol(key, arr._keys_multicol)
            if imk < 0:
                arr._keys += [key]
            else:
                if arr._keys_multicol[imk] not in arr._keys:
                    arr._keys += [arr._keys_multicol[imk]]
                    arr._keys_multicol_lookup[arr._keys_multicol[imk]] = [key]
                else:
                    arr._keys_multicol_lookup[arr._keys_multicol[imk]].append(
                        key)

        # TODO: would be nicer to make not use of __setitem__ here but instead
        # call the constructor with the data in the corresponding form
        # especially because this call __setitem__ of np.recarray, but of
        # BoundRecArray
        for i, name in enumerate(dtype.names):
            arr[name] = np.array(cols[i], dtype=dtype[name])

        return arr

    def __contains__(self, k):
        return k in set(self._keys)  # should be about as fast as for dict, pd.dataframe

    def keys(self):
        """Get keys of fields excluding the index."""
        # this replaces the columns property / also pd.dataframes yield
        # the keys as an index object upon calling .keys(), dicts anyway
        # this therefore emulates pd.dataframes and dicts well enough
        return self._keys

    def flipped(self):
        old_index_name = self._index_name
        new_index_name = SMP_INDEX if old_index_name == VAR_INDEX else VAR_INDEX

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

    def _multicol_view(self, arr, keys):
        """Get a multi-column view rather than a copy.

        This is consistent with the behavior of single columns.

        Will be implemented that way from numpy 1.13 on anyways, so we do not use
        in __getitem__, instead suppress the warning raised there.

        Note
        ----
        http://stackoverflow.com/questions/15182381/how-to-return-a-view-of-several-columns-in-numpy-structured-array
        """
        dtype2 = np.dtype({name: arr.dtype.fields[name] for name in keys})
        return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

    # TODO: __delitem__ should be aware of _keys_multicol and _keys

    def __getitem__(self, k):
        """Either a single one- or multi-column or mulitiple one-colum items."""
        import warnings  # ignore FutureWarning about multi-column access of structured arrays
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                return super(BoundRecArr, self).__getitem__(k)
            except ValueError:  # access multiple columns with key k
                keys = self._keys_multicol_lookup[k]
                return super(BoundRecArr, self).__getitem__(keys).view(
                    self.dtype[keys[0]]).reshape(self.shape + (len(self._keys_multicol_lookup[k]),))

    def __setitem__(self, keys, values):
        """Either a single one- or multi-column or mulitiple one-colum items."""
        if isinstance(keys, str):
            keys = [keys]
            # check if values is nested, if not, it's not multicolumn
            if (not hasattr(values[0], '__len__')
                or len(values[0]) == 1 or isinstance(values[0], str)):
                values = [values]
            else:  # otherwise add the multi_column_key
                if keys[0] not in self._keys_multicol:
                    key_multicol = keys[0]
                    self._keys_multicol += [key_multicol]
                    self._keys += [key_multicol]
                    # generate single-column keys
                    keys = _gen_keys_from_key_multicol(key_multicol, len(values[0]))
                    self._keys_multicol_lookup[key_multicol] = keys
                values = np.array(values)
                if values.shape[0] == self.shape[0]:
                    values = values.T
                else:
                    raise ValueError('You provided an array with {} rows but it need'
                                     'to have {}.'
                                     .format(values.shape[0], self.shape[0]))
        keys = np.array(keys)
        values = np.array(values)  # sequence of arrays or matrix with n_keys *rows*
        # update keys
        for key in keys:
            if (key != self._index_name
                and key not in self._keys
                and _key_belongs_to_which_key_multicol(key, self._keys_multicol) < 0):
                self._keys += [key]

        if len(keys) != len(values):
            raise ValueError('You passed {} column keys but {} arrays as columns. '
                             'If you passed a matrix instead of a sequence '
                             'of arrays, try transposing it.'
                             .format(len(keys), len(values)))

        if values.shape[1] != self.shape[0]:
            raise ValueError('You want to add a column with {} rows '
                             'but it need to have {} rows.'
                             .format(values.shape[1], self.shape[0]))

        present = np.intersect1d(keys, self.dtype.names)
        absent = np.setdiff1d(keys, self.dtype.names)

        if any(present):
            for k, v in zip(present, values[np.in1d(keys, present)]):
                super(BoundRecArr, self).__setitem__(k, v)

        if any(absent):
            attr = 'smp' if self._index_name == SMP_INDEX else 'var'
            if values.shape[1] > len(self):
                raise ValueError('New column has too many entries ({} > {})'
                                 .format(values.shape[1], len(self)))
            source = append_fields(self, absent, values[np.in1d(keys, absent)],
                                   usemask=False, asrecarray=True)
            new = BoundRecArr(source, self._index_name, self._parent,
                              keys_multicol=self._keys_multicol)
            setattr(self._parent, attr, new)

    def to_df(self):
        """Return pd.dataframe with index filled either with smp_names or var_names."""
        import pandas as pd
        return pd.DataFrame().from_records(self, index=self._index_name)


class AnnData(IndexMixin):

    def __init__(self, data, smp=None, var=None, add=None, dtype='float32', single_col=False):
        """Annotated Data

        Stores data matrix `X` of shape n_samples x n_variables, key-based
        access to sample and variable annotations of shapes n_samples x
        n_smp_keys and n_variables x n_var_keys and additional
        unstructured annotation in a dict.

        Sample and variable annotation but are a subclass of np.ndarray with
        access to a single or multiple columns with a single key. They provide
        the basic functionality of a pd.DataFrame, if you want all of it, use
        to_df().

        Parameters
        ----------
        data : dict, np.ndarray, np.ma.MaskedArray, sp.spmatrix
            The data matrix `X`
            X : np.ndarray, np.ma.MaskedArray, sp.spmatrix
                A n_samples x n_variables data matrix. Is flattened if either
                n_samples or n_variables is 1, so that numpys behavior is
                reproduced: `adata[:, 0].X == adata.X[:, 0]`.
            or a dict containing `X` as 'X' and possibly
            'row_names' / 'smp_names' : list, np.ndarray, optional
                A n_samples array storing names for samples.
            'col_names' / 'var_names' : list, np.ndarray, optional
                A n_variables array storing names for variables.
            'row' / 'smp' : dict, optional
                A dict with row annotation.
            'col' / 'var' : dict, optional
                A dict with column annotation.
        smp : np.ndarray, dict or None (default: None)
            A n_samples x n_smp_keys array containing sample names (`index`)
            and other sample annotation in the columns. A passed dict is
            converted to a record array.
        var : np.ndarray, dict or None (default: None)
            The same as `smp`, but of shape n_variables x n_var_keys for annotation of
            variables.
        add : dict (default: None)
            Unstructured annotation for the whole dataset.
        dtype : simple np.dtype, optional (default: float32)
            Convert data matrix to this type upon initialization.
        single_col : bool, optional (default: False)
            Interpret one-dimensional input array as column.

        Attributes
        ----------
        X, smp, var, add from the Parameters.
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

        X = X.astype(dtype)  # type conversion
        if X.dtype.names is None and len(X.shape) not in {0, 1, 2}:
            raise ValueError('X needs to be 2-dimensional, not '
                             '{}D'.format(len(X.shape)))

        self.X = X
        if len(self.X.shape) == 2:
            # TODO: int immutable, copy of references to ints in self.X.shape
            # only valid until accidental change
            self.n_smps, self.n_vars = self.X.shape
            # flatten to emulate numpys behavior upon slicing
            if self.n_smps == 1 and self.n_vars == 1:
                self.X = self.X[0, 0]
            elif self.n_smps == 1 or self.n_vars == 1:
                self.X = self.X.flatten()
        elif len(self.X.shape) == 1 and single_col:
            self.n_smps = self.X.shape[0]
            self.n_vars = 1
        elif len(self.X.shape) == 1:
            self.n_vars = self.X.shape[0]
            self.n_smps = 1
        else:
            self.n_vars = 1
            self.n_smps = 1

        smp_keys_multicol = None
        if add and 'smp_keys_multicol' in add:
            smp_keys_multicol = add['smp_keys_multicol']
        var_keys_multicol = None
        if add and 'var_keys_multicol' in add:
            var_keys_multicol = add['var_keys_multicol']

        self.smp = BoundRecArr(smp, SMP_INDEX, self, self.n_smps, smp_keys_multicol)
        self.var = BoundRecArr(var, VAR_INDEX, self, self.n_vars, var_keys_multicol)
        self._check_dimensions()

        self.add = add

    def __contains__(self, k):
        raise AttributeError("AnnData has no attribute __contains__, don't check `in adata`.")

    def __repr__(self):
        return ('AnnData object with attributes\n'
                '    X : array-like data matrix of shape\n'
                '    n_smps x n_vars = {} x {}\n'
                '    smp_names : the samples index {}\n'
                '    var_names : the variables index {}\n'
                '    smp_keys : keys to samples annotation {}\n'
                '    var_keys : keys to variables annotation {}\n'
                '    smp : key-based access to sample annotation (shape n_smps x n_smp_keys)\n'
                '    var : key-based access to variable annotation (shape n_vars x n_var_keys)\n'
                '    add : key-based access to unstructured annotation'
                .format(self.n_smps, self.n_vars, self.smp_names, self.var_names,
                        self.smp_keys(), self.var_keys()))

    @property
    def smp_names(self):
        """Samples index."""
        return self.smp[SMP_INDEX]

    @smp_names.setter
    def smp_names(self, keys):
        self.smp[SMP_INDEX] = keys

    @property
    def var_names(self):
        """Variables index."""
        return self.var[VAR_INDEX]

    @var_names.setter
    def var_names(self, keys):
        self.var[VAR_INDEX] = keys

    def smp_keys(self):
        """Return keys of sample annotation, excluding the index `smp_names`."""
        return self.smp.keys()

    def var_keys(self):
        """Return keys of variable annotation, excluding the index `var_names`."""
        return self.var.keys()

    def add_keys():
        """Return keys of addtional unstructured annotation."""
        return self.add.keys()

    def __setattr__(self, key, value):
        names_col = dict(smp=SMP_INDEX, var=VAR_INDEX).get(key)
        # if smp/ var is set, give it the right class
        if names_col and not isinstance(value, BoundRecArr):
            names_orig, dim = ((self.smp_names, 0) if names_col == SMP_INDEX
                               else (self.var_names, 1))
            value_orig, value = value, BoundRecArr(value, names_col, self)
            if len(value) != self.X.shape[dim]:
                raise ValueError('New value for {!r} was converted to a '
                                 'reacarray of length {} instead of {}'
                                 .format(key, len(value_orig), len(self)))
            if (value[names_col] == np.arange(self.X.shape[dim]).astype(str)).all():  # TODO: add to constructor
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
        adata = AnnData(X, smp_ann, var_ann, self.add)
        return adata

    def __setitem__(self, index, val):
        smp, var = self._normalize_indices(index)
        self.X[smp, var] = val

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

    def _check_dimensions(self):
        if len(self.smp) != self.n_smps:
            raise ValueError('Sample annotation needs to have the same amount of '
                             'rows as data has ({}), but has {} rows'
                             .format(self.n_smps, self.smp.shape[0]))
        if len(self.var) != self.n_vars:
            raise ValueError('Feature annotation needs to have the same amount of '
                             'rows as data has columns ({}), but has {} rows'
                             .format(self.n_vars, self.var.shape[0]))

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
        smp = OrderedDict([(k, self.smp[k]) for k in self.smp.dtype.names])
        var = OrderedDict([(k, self.var[k]) for k in self.var.dtype.names])
        d = {'X': self.X, 'smp': smp, 'var': var,
             'smp_names': self.smp_names, 'var_names': self.var_names}
        for k, v in self.add.items():
            d[k] = v
        d['smp_keys_multicol'] = self.smp._keys_multicol
        d['var_keys_multicol'] = self.var._keys_multicol
        return d


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(ma.array([[1, 2], [3, 4]]), add={'mask': [0, 1, 1, 0]})
    AnnData(sp.eye(2))
    AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(Smp=['A', 'B']),
        dict(Feat=['a', 'b', 'c']))

    assert AnnData(np.array([1, 2])).X.shape == (2,)

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


def test_creation_from_vector():
    adata = AnnData(np.array([1, 2, 3]))
    adata = AnnData(np.array([[1], [2], [3]]))


def test_slicing():
    adata = AnnData(np.array([[1, 2, 3],
                              [4, 5, 6]]))

    assert np.all(adata[:, 0].X == adata.X[:, 0])

    assert adata[0, 0].X.tolist() == 1
    assert adata[0, :].X.tolist() == [1, 2, 3]
    assert adata[:, 0].X.tolist() == [1, 4]

    assert adata[:, [0, 1]].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array([0, 2])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, np.array([False, True, True])].X.tolist() == [[2, 3], [5, 6]]
    assert adata[:, 1:3].X.tolist() == [[2, 3], [5, 6]]


def test_slicing_strings():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata['A', 'a'].X.tolist() == 1
    assert adata['A', :].X.tolist() == [1, 2, 3]
    assert adata[:, 'a'].X.tolist() == [1, 4]
    assert adata[:, ['a', 'b']].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array(['a', 'c'])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, 'b':'c'].X.tolist() == [[2, 3], [5, 6]]

    from pytest import raises
    with raises(IndexError): _ = adata[:, 'X']
    with raises(IndexError): _ = adata['X', :]
    with raises(IndexError): _ = adata['A':'X', :]
    with raises(IndexError): _ = adata[:, 'a':'X']


def test_get_subset_add():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                    dict(Smp=['A', 'B']),
                    dict(Feat=['a', 'b', 'c']))

    assert adata[0, 0].smp['Smp'].tolist() == ['A']
    assert adata[0, 0].var['Feat'].tolist() == ['a']


def test_transpose():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    adata1 = adata.T

    # make sure to not modify the original!
    assert adata.smp_names.tolist() == ['A', 'B']
    assert adata.var_names.tolist() == ['a', 'b', 'c']

    assert SMP_INDEX in adata1.smp.dtype.names
    assert VAR_INDEX in adata1.var.dtype.names
    assert adata1.smp_names.tolist() == ['a', 'b', 'c']
    assert adata1.var_names.tolist() == ['A', 'B']
    assert adata1.X.shape == adata.X.T.shape

    adata2 = adata.transpose()
    assert np.array_equal(adata1.X, adata2.X)
    assert np.array_equal(adata1.smp, adata2.smp)
    assert np.array_equal(adata1.var, adata2.var)


def test_append_add_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.smp['new'] = [1, 2]
    adata.smp[['new2', 'new3']] = [['A', 'B'], ['c', 'd']]

    from pytest import raises
    with raises(ValueError):
        adata.smp['new4'] = 'far too long'.split()


def test_set_add():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.smp = dict(smp_names=['1', '2'])
    assert isinstance(adata.smp, BoundRecArr)
    assert len(adata.smp.dtype) == 1

    adata.smp = dict(a=[3, 4])  # keep smp_names and add a custom column
    assert isinstance(adata.smp, BoundRecArr)
    assert len(adata.smp.dtype) == 2
    assert adata.smp_names.tolist() == ['1', '2']  # still the same names

    from pytest import raises
    with raises(ValueError):
        adata.smp = dict(a=[1, 2, 3])


def test_print():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                    dict(foo=['A', 'B']),
                    dict(bar=['a', 'b', 'c']))
    print(adata)


def test_multicol_getitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp[['a', 'b']] = np.array([[0, 1], [2, 3]])
    print(adata.smp[['a', 'b']])


def test_multicol_single_key_setitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp['c'] = np.array([[0, 1], [2, 3]])
    print(adata.smp.dtype.names)
    print(adata.smp['c'])


def test_boundrecarray_keys():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp['foo'] = np.array([[0, 1], [2, 3]])
    assert adata.smp_keys() == ['foo']
    assert adata.smp.keys() == ['foo']
    assert adata.smp.dtype.names == ('smp_names', 'foo1of2', 'foo2of2')

    adata.smp['d'] = np.array([[0, 1], [2, 3]])
    assert adata.smp.keys() == ['foo', 'd']
    assert 'd' in adata.smp
    from pytest import raises
    with raises(KeyError):
        adata.smp['e']


def test_n_smps():
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    assert adata.n_smps == 3
    adata1 = adata[:2, ]
    assert adata1.n_smps == 2

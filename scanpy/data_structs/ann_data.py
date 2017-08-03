"""Annotated Data
"""
from collections import Mapping, Sequence
from collections import OrderedDict
from enum import Enum

import numpy as np
from numpy import ma
from numpy.lib.recfunctions import append_fields, rec_drop_fields
import pandas as pd
from scipy import sparse as sp
from scipy.sparse.sputils import IndexMixin

from .. import logging as logg
from ..utils import merge_dicts

SMP_INDEX = 'smp_names'
VAR_INDEX = 'var_names'

STRING_TYPE = 'S50'

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
    keys = [('{}{:03}of{:03}')
            .format(key_multicol, i+1, n_keys) for i in range(n_keys)]
    return keys


class SetKeyError(ValueError):
    template = '''Currently you cannot implicitly reallocate memory:
Setting the array for key {} with dtype {} requires too much memory, \
you should init AnnData with a large enough data type from the beginning.
Probably you try to assign a string of length {} \
although the array can only store strings of length {}.'''

    def __init__(self, k, dtype, dtype_expected):
        msg = SetKeyError.template.format(k, dtype, int(dtype.itemsize), int(dtype_expected.itemsize))
        super(ValueError, self).__init__(msg)


class BoundStructArray(np.ndarray):

    def __new__(cls, source, index_key, is_attr_of, n_row=None,
                keys_multicol=None, new_index_key=None):
        """Dimensionally structured dict, lowlevel alternative to pandas dataframe.

        Behaves like a dict except that
        - data has to match shape constraints
        - slicing in the row-dimension is possible

        Behaves like a numpy array except that:
        - single and multiple columns can be accessed with keys (strings)
        - the `index` column is hidden from the user and enables sclicing in AnnData
        - you can add new columns via [] (`__setitem__`)
        - it is bound to AnnData

        Can be exported to a pandas dataframe via `.to_df()`.

        Parameters
        ----------
        index_key : str
            The key or field name that shall be used as the index. If not found
            generate a dummy index np.array(['0', '1', '2' ...]).
        is_attr_of : (object, x)
            Tuple `(o, x)` containing the object `o` to which the instance is
            bound and the name `x` of the attribute.
        new_index_key : str or None
            Only needed if the internal index_key of the object should differ from
            `index_key`.

        Attributes
        ----------
        index : np.ndarray
            Access to the index column. Can also be accessed via ['index'].
        """
        old_index_key = index_key
        new_index_key = index_key if new_index_key is None else new_index_key
        # create from existing array
        if isinstance(source, np.ndarray):
            # create from existing BoundStructArray
            if isinstance(source, BoundStructArray):
                keys_multicol = source._keys_multicol if keys_multicol is None else keys_multicol
            # we need to explicitly make a deep copy of the dtype
            arr = np.array(source, dtype=[t for t in source.dtype.descr])
            # rename the index
            arr.dtype.names = tuple(new_index_key if n == old_index_key else n
                                    for n in source.dtype.names)
        # create from None or Dict
        else:
            if source is None:  # empty array
                cols = [np.arange(n_row).astype(STRING_TYPE)]
                dtype = [(new_index_key, cols[0].dtype)]
            else:
                if not isinstance(source, Mapping):
                    raise ValueError('Expected np.ndarray or dictlike type, not {}.'
                                     .format(type(source)))
                # meta is dict-like
                names = list(source.keys())
                try:  # transform to byte-strings
                    cols = [np.asarray(col) if np.array(col[0]).dtype.char not in {'U', 'S'}
                            else np.asarray(col).astype(STRING_TYPE) for col in source.values()]
                except UnicodeEncodeError:
                    raise ValueError('Currently only support ascii strings. Don\'t use "ö" etc. for sample annotation.')

                if old_index_key not in source:
                    names.append(new_index_key)
                    cols.append(np.arange(len(cols[0]) if cols else n_row).astype(STRING_TYPE))
                else:
                    names[names.index(old_index_key)] = new_index_key
                    cols[names.index(old_index_key)] = cols[names.index(old_index_key)].astype(STRING_TYPE)
                dtype = list(zip(names, [str(c.dtype) for c in cols]))
            try:
                dtype = np.dtype(dtype)
            except TypeError:
                # TODO: fix compat with Python 2
                # print(dtype, file=sys.stderr)
                raise

            arr = np.zeros((len(cols[0]),), dtype)
            # here, we do not want to call BoundStructArray.__getitem__
            # but np.ndarray.__getitem__, therefore we avoid the following line
            # arr = np.ndarray.__new__(cls, (len(cols[0]),), dtype)
            for i, name in enumerate(dtype.names):
                arr[name] = np.array(cols[i], dtype=dtype[name])

        # generate an instance of BoundStructArray
        arr = np.asarray(arr).view(cls)
        # the index_key used to look up the index column
        arr.index_key = new_index_key
        #
        arr._is_attr_of = is_attr_of
        # only the multicol keys
        arr._keys_multicol = ([] if keys_multicol is None
                              else list(keys_multicol))
        # fast lookup of single-column keys
        arr._keys_multicol_lookup = ({} if keys_multicol is None
                                     else dict.fromkeys(keys_multicol))
        # all keys (multi- and single-col keys that do not belong to multicol)
        arr._keys = []  # this excludes self.index_key

        # initialize arr._keys
        loop_over_keys = (key for key in arr.dtype.names if key != index_key)
        for key in loop_over_keys:
            imk = _key_belongs_to_which_key_multicol(key, arr._keys_multicol)
            if imk < 0:
                arr._keys += [key]
            else:
                if arr._keys_multicol[imk] not in arr._keys:
                    arr._keys += [arr._keys_multicol[imk]]
                    arr._keys_multicol_lookup[arr._keys_multicol[imk]] = [key]
                else:
                    arr._keys_multicol_lookup[arr._keys_multicol[imk]].append(key)

        # check multicol keys
        for key in arr._keys_multicol:
            if (int(arr._keys_multicol_lookup[key][-1].split('of')[-1])
                != len(arr._keys_multicol_lookup[key])):
                logg.error('Different formats in {}. Expect weird behavior! '
                           'Recompute {}!'.format(arr._keys_multicol_lookup[key], key))

        return arr

    def __array_finalize__(self, arr):
        # determines the behavior of views
        if arr is None: return
        self.index_key = getattr(arr, 'index_key', None)
        self._is_attr_of = getattr(arr, '_is_attr_of', None)
        self._keys = getattr(arr, '_keys', None)
        self._keys_multicol = getattr(arr, '_keys_multicol', None)
        self._keys_multicol_lookup = getattr(arr, '_keys_multicol_lookup', None)

    @property
    def index(self):
        return self[self.index_key]

    @index.setter
    def index(self, names):
        self[self.index_key] = names

    def __contains__(self, k):
        return k in set(self._keys)  # should be about as fast as for dict, pd.dataframe

    def keys(self):
        """Get keys of fields excluding the index (same as `columns()`)."""
        # this replaces the columns property / also pd.dataframes yield
        # the keys as an index object upon calling .keys(), dicts anyway
        # this therefore emulates pd.dataframes and dicts well enough
        return self._keys

    def columns(self):
        """Get keys of fields excluding the index (same as `keys()`)."""
        return self._keys

    def delete_field(self, name):
        """Delete field with name."""
        if name not in self.dtype.names:
            raise ValueError('Currently, can only delete single names from {}.'
                             .format(self.dtype.names))
        new_array = rec_drop_fields(self, name)
        new = BoundStructArray(new_array, self.index_key, self._is_attr_of,
                               keys_multicol=self._keys_multicol)
        setattr(self._is_attr_of[0], self._is_attr_of[1], new)

    def copy(self):
        return BoundStructArray(self, self.index_key,
                                self._is_attr_of,
                                len(self),
                                keys_multicol=self._keys_multicol.copy())

    def copy_index_exchanged(self):
        new_index_key, new_attr = (SMP_INDEX, 'smp') if self.index_key == VAR_INDEX else (VAR_INDEX, 'var')
        return BoundStructArray(self, self.index_key,
                                (self._is_attr_of[0], new_attr),
                                len(self),
                                keys_multicol=self._keys_multicol.copy(),
                                new_index_key=new_index_key)

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
        """Either a single one- or multi-column or mulitiple one-colum items.

        Column slices yield conventional *homogeneous* numpy arrays, as is
        useful for the user.

        Row slices are BoundStructArrays.
        """
        import warnings  # ignore FutureWarning about multi-column access of structured arrays
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                # determine if we are retrieving a column
                # TODO: rewrite this in a nicer way!
                col_slice = False
                if isinstance(k, str):
                    col_slice = True
                    shape = self.shape
                    dtype = self.dtype[k]
                elif isinstance(k, list) and isinstance(k[0], str):
                    col_slice = True
                    shape = self.shape + (len(k),)
                    dtype = self.dtype[k[0]]
                elif isinstance(k, np.ndarray) and len(k.shape) > 0:
                    if k.shape == (0,):  # empty array
                        print('warning: slicing with empty list returns None')
                        return None
                    if isinstance(k[0], str):
                        col_slice = True
                        shape = self.shape + (len(k),)
                        dtype = self.dtype[k[0]]
                view = super(BoundStructArray, self).__getitem__(k)
                if col_slice:
                    view = view.view(dtype, type=np.ndarray).reshape(shape)
            except (ValueError, KeyError):  # access multiple columns with key k
                keys = self._keys_multicol_lookup[k]
                view = super(BoundStructArray, self).__getitem__(keys).view(
                    self.dtype[keys[0]]).reshape(
                        self.shape + (len(self._keys_multicol_lookup[k]),))
            if view.dtype.char == 'S':
                return view.astype('U')
            return view

    def __setitem__(self, keys, values):
        """Either a single one- or multi-column or mulitiple one-colum items."""
        if isinstance(values, pd.Series):
            values = values.values
            if values.dtype.char == 'O' or 'int' in values.dtype.name:
                values = values.astype(str)

        names_to_remove = []
        if isinstance(keys, str):
            # TODO: check that no-one accidentally overwrites the index?
            # quite unlikely though as self.index_key is not common
            # if keys == self.index_key:
            #     raise ValueError('The key {} is reserved for the index in BoundStructArray. '
            #                      .format(self.index_key))
            keys = [keys]
            # check if values is nested, if not, it's not multicolumn
            if (not hasattr(values[0], '__len__')
                or len(values[0]) == 1  # seems that we do not need this, as the previous line matches already
                or np.array(values[0]).dtype.char in {'S', 'U'}):  # a string is passed
                values = [values]
            else:  # otherwise it's a multicolumn key
                key_multicol = keys[0]
                if keys[0] not in self._keys_multicol:
                    self._keys_multicol += [key_multicol]
                    self._keys += [key_multicol]
                # generate single-column keys
                keys = _gen_keys_from_key_multicol(key_multicol, len(values[0]))
                self._keys_multicol_lookup[key_multicol] = keys
                # remove all fields from the array that are not among keys
                keys_set = set(keys)
                for name in self.dtype.names:
                    if name.startswith(key_multicol) and name not in keys_set:
                        names_to_remove.append(name)
                values = np.array(values)
                if values.shape[0] == self.shape[0]:
                    values = values.T
                else:
                    raise ValueError('You provided an array with {} rows but it need'
                                     'to have {}.'
                                     .format(values.shape[0], self.shape[0]))
        else:
            values = np.array(values)  # sequence of arrays or matrix with n_keys *rows*
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
            if (key != self.index_key
                and key not in self._keys
                and _key_belongs_to_which_key_multicol(key, self._keys_multicol) < 0):
                self._keys += [key]

        if len(keys) != len(values):
            print(keys, values)
            raise ValueError('You passed {} column keys but {} arrays as columns. '
                             'If you passed a matrix instead of a sequence '
                             'of arrays, try transposing it.'
                             .format(len(keys), len(values)))

        if values.shape[1] != self.shape[0]:
            raise ValueError('You want to add a column with {} rows '
                             'but it need to have {} rows.'
                             .format(values.shape[1], self.shape[0]))

        if values.dtype.char in {'U', 'S'}:
            try:
                itemsize = values.dtype.itemsize
                if values.dtype.char == 'U': itemsize /= 4
                if itemsize > np.dtype(STRING_TYPE).itemsize:
                    logg.m('WARNING: truncating strings to length {}'
                           .format(np.dtype(STRING_TYPE).itemsize))
                values = values.astype(STRING_TYPE)
            except UnicodeEncodeError:
                raise ValueError('Currently only support ascii strings. Don\'t use "ö" etc. for sample annotation.')

        present = np.intersect1d(keys, self.dtype.names)
        absent = np.setdiff1d(keys, self.dtype.names)

        if any(present):
            for k, v in zip(present, values[np.in1d(keys, present)]):
                if (v.dtype != self.dtype[k]
                        and v.dtype.itemsize > self.dtype[k].itemsize):
                    # TODO: need to reallocate memory
                    # or allow storing objects, or use pd.dataframes
                    raise SetKeyError(k, v.dtype, self.dtype[k])
                super(BoundStructArray, self).__setitem__(k, v)

        if any(absent):
            if values.shape[1] > len(self):
                raise ValueError('New column has too many entries ({} > {})'
                                 .format(values.shape[1], len(self)))
            source = append_fields(self, absent, values[np.in1d(keys, absent)],
                                   usemask=False, asrecarray=True)
            if names_to_remove:
                source = rec_drop_fields(source, names_to_remove)
            new = BoundStructArray(source, self.index_key, self._is_attr_of,
                                   keys_multicol=self._keys_multicol)
            setattr(self._is_attr_of[0], self._is_attr_of[1], new)

    def __str__(self):
        return self.to_df().__str__()

    def to_df(self):
        """Return pd.dataframe with index filled either with smp_names or var_names."""
        import pandas as pd
        return pd.DataFrame().from_records(self, index=self.index_key)


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
        X : array-like data matrix
        n_smps : number of samples
        n_vars : number of variables
        smp_names : sample names/index
        var_names : variable names/index
        smp_keys : keys to samples annotation
        var_keys : keys to variables annotation
        add_keys : keys to unstructured annotation
        smp : sample annotation (shape n_smps x n_smp_keys)
        var : variable annotation (shape n_vars x n_var_keys)
        add : unstructured annotation
        """
        if isinstance(data, Mapping):
            if any((smp, var, add)):
                raise ValueError('If `data` is a dict no further arguments must be provided.')
            X, smp, var, add = self.from_dict(data)
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

        # type conversion: if type doesn't match, a copy is made
        if sp.issparse(X) or isinstance(X, ma.MaskedArray):
            # TODO: maybe use view on data attribute of sparse matrix
            #       as in readwrite.read_10x_h5
            if X.dtype != np.dtype(dtype): X = X.astype(dtype)
        else:  # is np.ndarray
            X = X.astype(dtype, copy=False)

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
                if sp.issparse(self.X): self.X = self.X.toarray()
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

        self.smp = BoundStructArray(smp, SMP_INDEX, (self, 'smp'), self.n_smps, smp_keys_multicol)
        self.var = BoundStructArray(var, VAR_INDEX, (self, 'var'), self.n_vars, var_keys_multicol)
        self._check_dimensions()

        self.add = add or {}

    def __contains__(self, k):
        raise AttributeError("AnnData has no attribute __contains__, don't check `in adata`.")

    def __repr__(self):
        return ('AnnData object with n_smps x n_vars = {} x {}\n'
                '    smp_keys = {}\n'
                '    var_keys = {}\n'
                '    add_keys = {}\n'
                '    smp_names = {}\n'
                '    var_names = {}'
                .format(self.n_smps, self.n_vars,
                        self.smp.keys(), self.var.keys(), list(self.add.keys()),
                        self.smp_names, self.var_names))

    @property
    def smp_names(self):
        """Samples index."""
        return self.smp.index

    @smp_names.setter
    def smp_names(self, names):
        self.smp.index = names

    @property
    def var_names(self):
        """Variables index."""
        return self.var.index

    @var_names.setter
    def var_names(self, names):
        self.var.index = names

    def smp_keys(self):
        """Return keys of sample annotation, excluding the index `smp_names`."""
        return self.smp.keys()

    def var_keys(self):
        """Return keys of variable annotation, excluding the index `var_names`."""
        return self.var.keys()

    def add_keys(self):
        """Return keys of addtional unstructured annotation."""
        return self.add.keys()

    def __setattr__(self, key, value):
        # if smp/ var is set, make it a BoundStructArray
        if key in {'smp', 'var'} and not isinstance(value, BoundStructArray):
            index_key, names_orig, dim = ((SMP_INDEX, self.smp_names, 0) if key == 'smp'
                                          else (VAR_INDEX, self.var_names, 1))
            value_orig, value = value, BoundStructArray(value, index_key, self)
            if len(value) != self.X.shape[dim]:
                raise ValueError('New value for {!r} was converted to a '
                                 'reacarray of length {} instead of {}'
                                 .format(key, len(value_orig), len(self)))
            if (value[index_key] == np.arange(self.X.shape[dim]).astype(str)).all():  # TODO: add to constructor
                value[index_key] = names_orig
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
                i = np.where(names == i)[0]
                if len(i) == 0:  # returns array of length 0 if nothing is found
                    raise IndexError('Name "{}" is not valid variable or sample index.'
                                     .format(index))
                i = i[0]
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
        # Note: this cannot be made inplace
        # http://stackoverflow.com/questions/31916617/using-keyword-arguments-in-getitem-method-in-python
        smp, var = self._normalize_indices(index)
        X = self.X[smp, var]
        smp_ann = self.smp[smp]
        var_ann = self.var[var]
        assert smp_ann.shape[0] == X.shape[0], (smp, smp_ann)
        assert var_ann.shape[0] == X.shape[1], (var, var_ann)
        add_ann = self.add.copy()
        # slice sparse spatrices of n_smps x n_smps in self.add
        if not (isinstance(smp, slice) and
                smp.start is None and smp.step is None and smp.stop is None):
            raised_warning = False
            for k, v in self.add.items():  # TODO: make sure this really works as expected
                if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                    add_ann[k] = v.tocsc()[:, smp].tocsr()[smp, :]
                    if not raised_warning:
                        logg.warn('Slicing adjacency matrices can be dangerous. '
                                  'Consider recomputing the data graph.')
                        raised_warning = True
        adata = AnnData(X, smp_ann, var_ann, add_ann)
        return adata

    def inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        self.X = self.X[:, index]
        self.var = self.var[index]
        self.n_vars = self.X.shape[1]
        return None

    def inplace_subset_smp(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        self.X = self.X[index, :]
        self.smp = self.smp[index]
        raised_warning = False
        for k, v in self.add.items():
            if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                self.add[k] = v.tocsc()[:, index].tocsr()[index, :]
                if not raised_warning:
                    logg.warn('Slicing adjacency matrices can be dangerous. '
                              'Consider recomputing the data graph.')
                    raised_warning = True
        self.n_smps = self.X.shape[0]
        return None

    def get_smp_array(self, k):
        """Get an array along the sample dimension by first looking up
        smp_keys and then var_names."""
        x = (self.smp[k] if k in self.smp_keys()
             else self[:, k].X if k in set(self.var_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in smp_keys or var_names.'
                             .format(k))
        return x

    def get_var_array(self, k):
        """Get an array along the variables dimension by first looking up
        var_keys and then smp_names."""
        x = (self.var[k] if k in self.var_keys()
             else self[k] if k in set(self.smp_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in var_keys or smp_names.'
                             .format(k))
        return x

    def __setitem__(self, index, val):
        smp, var = self._normalize_indices(index)
        self.X[smp, var] = val

    def __len__(self):
        return self.X.shape[0]

    def transpose(self):
        """Return a transposed view of the object.

        Sample axis (rows) and variable axis are interchanged. No additional memory.
        """
        if sp.isspmatrix_csr(self.X):
            return AnnData(self.X.T.tocsr(), self.var.copy_index_exchanged(), self.smp.copy_index_exchanged(), self.add)
        return AnnData(self.X.T, self.var.copy_index_exchanged(), self.smp.copy_index_exchanged(), self.add)

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

    def to_dict(self):
        """A dict that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        d = {'X': self.X, 'smp': self.smp, 'var': self.var}
        for k, v in self.add.items():
            d[k] = v
        d['smp_keys_multicol'] = self.smp._keys_multicol
        d['var_keys_multicol'] = self.var._keys_multicol
        return d

    def from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.
        """
        add = dict(ddata.items())
        del ddata
        X = add['X']
        del add['X']
        if ('smp' in add and isinstance(add['smp'], np.ndarray)
            and 'var' in add and isinstance(add['var'], np.ndarray)):
            smp = add['smp']
            del add['smp']
            var = add['var']
            del add['var']
        else:
            smp, var = OrderedDict(), OrderedDict()
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

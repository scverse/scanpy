from functools import singledispatch
from typing import (
    Optional,
    Iterable,
    AbstractSet,
    Sequence,
    Tuple,
    Union,
    Literal,
    List,
    get_args,
)

from anndata import AnnData, utils
import numpy as np
import pandas as pd
import collections.abc as cabc
from scipy.sparse import coo_matrix, dia_matrix, spmatrix

Array = Union[np.ndarray, spmatrix]
AggType = Literal['count', 'mean', 'sum', 'var']


class Aggregate:
    """\
    Functionality for generic grouping and aggregating.

    There is currently support for count, sum, mean, and variance.

    Set `weight` for weighted sum, mean, and variance.

    Set `key_set` to a list of keys to most efficiently compute results for a subset of groups.

    **Implementation**

    Moments are computed using weighted sum aggregation of data by some feature
    via multiplication by a sparse coordinate matrix A, exposed by
    `_sparse_aggregator`. The approach works with data in ndarray or scipy sparse formats, with
    no view or copy overhead on runtime or memory, even when filtering keys.

    Runtime is effectively computation of the product A * X, i.e. the count of (non-zero)
    entries in X with multiplicity the number of group memberships for that entry. This is
    O(data) for partitions (each observation belonging to exactly one group), independent of
    the number of groups.

    Params
    ------
    groupby
        `Series` containing values for grouping by.
    data
        Data matrix for aggregation.
    weight
        Weights to be used for aggergation.
    key_set
        Subset of keys to which to filter.
    """

    _groupby: pd.Series
    _data: Array
    _weight: Union[pd.Series, Array]
    _key_set: AbstractSet[str]
    _key_index: Optional[np.ndarray]  # caution, may be stale if attributes are updated

    def __init__(
        self,
        groupby: pd.Series,
        data: Array,
        weight: Union[pd.Series, Array] = None,
        key_set: Optional[Iterable[str]] = None,
    ):
        self._groupby = groupby
        self._data = data
        self._weight = weight
        self._key_set = None if key_set is None else dict.fromkeys(key_set).keys()
        self._key_index = None

    def count(self) -> np.ndarray:
        """
        Count the number of observations in each group.

        Returns
        -------
        Array of counts.
        """
        _, key_index, _, _ = self._extract_indices()
        count_ = np.bincount(key_index)
        return count_

    def sum(self) -> Array:
        """
        Compute the sum per feature per group of observations.

        Returns
        -------
        Array of sum.
        """
        A, _ = self._sparse_aggregator(normalize=False)
        return utils.asarray(A * self._data)

    def mean(self) -> Array:
        """\
        Compute the mean per feature per group of observations.

        Returns
        -------
        Array of mean.
        """
        A, _ = self._sparse_aggregator(normalize=True)
        return utils.asarray(A * self._data)

    def count_mean_var(self, dof: int = 1) -> dict:
        """\
        Compute the count, as well as mean and variance per feature, per group of observations.

        The formula `Var(X) = E(X^2) - E(X)^2` suffers loss of precision when the variance is a
        very small fraction of the squared mean. In particular, when X is constant, the formula may
        nonetheless be non-zero. By default, our implementation resets the variance to exactly zero
        when the computed variance, relative to the squared mean, nears limit of precision of the
        floating-point significand.

        Params
        ------
        dof
            Degrees of freedom for variance.

        Returns
        -------
        dict with mean, count, and var keys.
        """
        assert dof >= 0
        A, _ = self._sparse_aggregator(normalize=True)
        count_ = np.bincount(self._key_index)
        mean_ = utils.asarray(A @ self._data)
        # sparse matrices do not support ** for elementwise power.
        mean_sq = utils.asarray(A @ _power(self._data, 2))
        if self._weight is None:
            sq_mean = mean_**2
        else:
            A_unweighted, _ = Aggregate(
                groupby=self._groupby,
                data=self._data,
                weight=self._weight,
                key_set=self._key_set,
            )._sparse_aggregator()
            mean_unweighted = utils.asarray(A_unweighted * self._data)
            sq_mean = 2 * mean_ * mean_unweighted + mean_unweighted**2
        var_ = mean_sq - sq_mean
        # TODO: Why these values exactly? Because they are high relative to the datatype?
        # (unchanged from original code: https://github.com/scverse/anndata/pull/564)
        precision = 2 << (42 if self._data.dtype == np.float64 else 20)
        # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
        var_[precision * var_ < sq_mean] = 0
        if dof != 0:
            var_ *= (count_ / (count_ - dof))[:, np.newaxis]
        return {'mean': mean_, 'var': var_, 'count': count_}

    def _sparse_aggregator(
        self, normalize: bool = False
    ) -> Tuple[coo_matrix, np.ndarray]:
        """
        Form a coordinate-sparse matrix A such that rows of A * X
        are weighted sums of groups of rows of X.

        A[i, j] = w includes X[j,:] in group i with weight w.

        Params
        ------
        normalize
            If true, weights for each group are normalized to sum to 1.0,
            corresponding to (weighted) mean.

        Returns
        -------
        A
            weighted sums of groups of rows of X.
        keys
            An ndarray with keys[i] the group key corresponding to row i of A.
        """
        keys, key_index, df_index, weight_value = self._extract_indices()
        if df_index is None:
            df_index = np.arange(len(key_index))
        if self._weight is None:
            weight_value = np.ones(len(key_index))
        # TODO: why a coo matrix here and a dia matrix below? (unchanged from original code: https://github.com/scverse/anndata/pull/564)
        A = coo_matrix(
            (weight_value, (key_index, df_index)),
            shape=(len(keys), self._data.shape[0]),
        )
        if normalize:
            n_row = A.shape[0]
            row_sums = np.asarray(A.sum(axis=1))
            D = dia_matrix(((row_sums.T**-1), [0]), shape=(n_row, n_row))
            A = D * A
        return A, keys

    def _filter_indices(
        self,
        keys: np.ndarray,
        key_index: np.ndarray,
        df_index: np.ndarray,
        weight_value: Optional[Union[pd.Series, Array]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, Union[pd.Series, Array, None]]:
        """Filter the values of keys, key_index, df_index, and optionally weight_value based on self._key_set.

        Parameters
        ----------
        keys
            Unique key values to be filtered.
        key_index
            Non-unique integer indices mapping keys to the df_index to be filtered.
        df_index
            An Index that the keys + key_index constitute to be filtered.
        weight_value, optional
            Weight values to be filtered., by default None

        Returns
        -------
            Filtered versions of all arguments.

        Raises
        ------
        ValueError
           If no keys in key_set found in keys.
        """
        keep = [i for i, k in enumerate(keys) if k in set(self._key_set)]
        if len(keep) == 0:
            raise ValueError("No keys in key_set found in keys.")
        elif len(keep) < len(keys):
            mask = np.in1d(key_index, keep)
            remap = np.zeros(len(keys), dtype=np.int64)
            for i, j in enumerate(keep):
                remap[j] = i
            keys = [keys[j] for j in keep]
            key_index = np.array([remap[i] for i in key_index[mask]], dtype=np.int64)
            df_index = df_index[mask]
            if weight_value is not None:
                weight_value = weight_value[mask]
        return keys, key_index, df_index, weight_value

    def _extract_indices(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, Union[pd.Series, Array, None]]:
        """Extract indices from self._groupby with the goal of building a matrix that can be multiplied with the data to produce an aggregation statistics e.g., mean or variance.
        These are filtered if a self._key_set is present.

        Returns
        -------
            Unique keys, an array mapping those unique keys to an index, said index, and a weight if present.
        """
        key_value = self._groupby
        keys, key_index = np.unique(_ndarray_from_seq(key_value), return_inverse=True)
        df_index = np.arange(len(key_index))
        if self._weight is None:
            weight_value = None
        else:
            weight_value = self._weight.values[df_index]
        if self._key_set is not None:
            keys, key_index, df_index, weight_value = self._filter_indices(
                keys, key_index, df_index, weight_value
            )
        self._key_index = key_index  # passed to count and count_mean_var to avoid re-extracting in the latter
        return keys, key_index, df_index, weight_value


def _power(X: Array, power: Union[float, int]) -> Array:
    """Generate elementwise power of a matrix.  Needed for sparse matrices because they do not support ** so the `.power` function is used.

    Parameters
    ----------
    X
        Matrix whose power is to be raised.
    power
        Integer power value

    Returns
    -------
        Matrix whose power has been raised.
    """
    return X ** power if isinstance(X, np.ndarray) else X.power(power)


def _ndarray_from_seq(lst: Sequence):
    # prevents expansion of iterables as axis
    n = len(lst)
    if n > 0 and isinstance(lst[0], cabc.Iterable):
        arr = np.empty(n, dtype=object)
        arr[:] = lst
    else:
        arr = np.array(lst)
    return arr


def _superset_columns(df: pd.DataFrame, groupby_key: str) -> List[str]:
    """\
    Find all columns which are a superset of the key column.

    Params
    ------
    df
        DataFrame which contains candidate columns.
    groupby_key
        Key for column of which to find superset of columns.

    Returns
    -------
    Superset columns.
    """
    columns = []
    groupy_key_codes = df[groupby_key].astype('category')
    for key in df:
        if key != groupby_key:
            key_codes = df[key].astype('category')
            if all(
                key_codes[groupy_key_codes == group_key_code].nunique() == 1
                for group_key_code in groupy_key_codes
            ):
                columns += [key]
    return columns


def _df_grouped(df: pd.DataFrame, key: str, key_set: List[str]) -> pd.DataFrame:
    """\
    Generate a grouped-by dataframe (no aggregation) by
    a key with columns that are supersets of the key column.

    Params
    ------
    df
        DataFrame to be grouped.
    key
        Column to be grouped on.
    key_set
        Values in the `key` column to keep before groupby.

    Returns
    -------
    pd.DataFrame: Grouped-by Dataframe.
    """
    df = df.copy()
    if key_set is not None:
        df = df[df[key].isin(key_set)]
    if pd.api.types.is_categorical_dtype(df[key]):
        df[key] = df[key].cat.remove_unused_categories()
    return df.groupby(key).first()[_superset_columns(df, key)]


@singledispatch
def aggregated(
    adata: AnnData,
    by: str,
    func: Union[AggType, Iterable[AggType]],
    *,
    dim: Literal['obs', 'var'] = 'obs',
    weight_key: Optional[str] = None,
    key_set: Optional[Iterable[str]] = None,
    dof: int = 1,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
    varm: Optional[str] = None,
) -> AnnData:
    """\
    Aggregate data based on one of the columns of one of the axes (`obs` or `var`).
    If none of `layer`, `obsm`, or `varm` are passed in, `X` will be used for aggregation data.
    If `func` only has length 1 or is just an `AggType`, then aggregation data is written to `X`.
    Otherwise, it is written to `layers` or `xxxm` as appropriate for the dimensions of the aggregation data.

    Params
    ------
    adata
        :class:`~anndata.AnnData` to be aggregated.
    by
        Key of the column to be grouped-by.
    func
        How to aggregate.
    dim
        Axis on which to find group by column.
    weight_key
        Key of the `dim` containing weights for a weighted sum aggregation.
    key_set
        Subset of dim on which to filter.
    dof
        Degrees of freedom for variance. Defaults to 1.
    layer
        If not None, key for aggregation data.
    obsm
        If not None, key for aggregation data.
    varm
        If not None, key for aggregation data.

    Returns
    -------
    Aggregated :class:`~anndata.AnnData`.
    """
    data = adata.X
    write_to_xxxm = None
    if sum([varm is None, obsm is None, layer is None]) > 1:
        raise TypeError("Please only provide one (or none) of varm, obsm, or layer")
    if not varm is None:
        data = adata.varm[varm]
        write_to_xxxm = True  # the data will have to be transposed so this is accurate
    elif not obsm is None:
        data = adata.obsm[obsm]
        write_to_xxxm = True
    elif not layer is None:
        data = adata.layers[layer]
        if dim == 'var':
            data = data.T
    elif dim == 'var':
        # i.e., all of `varm`, `obsm`, `layers` are None so we use `X` which must be transposed
        data = data.T
    return aggregated(
        data,
        groupby_df=getattr(adata, dim),
        dim=dim,
        by=by,
        write_to_xxxm=write_to_xxxm,
        no_groupby_df=getattr(adata, 'var' if dim == 'obs' else 'obs'),
        weight_key=weight_key,
        key_set=key_set,
        func=func,
        dof=dof,
    )


@aggregated.register(np.ndarray)
@aggregated.register(spmatrix)
def aggregated_from_array(
    data,
    groupby_df: pd.DataFrame,
    func: Union[AggType, Iterable[AggType]],
    dim: str,
    by: str,
    write_to_xxxm: bool,
    no_groupby_df: pd.DataFrame,
    weight_key: Optional[str] = None,
    key_set: Optional[Iterable[str]] = None,
    dof: int = 1,
) -> AnnData:
    """Aggregate data based on one of the columns of one of a `~pd.DataFrame`."""
    groupby = Aggregate(
        groupby=groupby_df[by],
        data=data,
        weight=groupby_df[weight_key] if weight_key is not None else None,
        key_set=key_set,
    )
    # groupby df is put in `obs`, nongroupby in `var` to be transposed later as appropriate
    obs_var_dict = {'obs': _df_grouped(groupby_df, by, key_set), 'var': no_groupby_df}
    data_dict = {
        'layers': {},
        'X': None,
        'obsm': {},
    }
    write_key = 'obsm' if write_to_xxxm else 'layers'
    funcs = set([func] if isinstance(func, str) else func)
    if unknown := funcs - set(get_args(AggType)):
        raise ValueError(f'… {unknown} …')
    if 'sum' in funcs:  # sum is calculated separately from the rest
        agg = groupby.sum()
        # put aggregation in X if it is the only one and the aggregation data is not coming from `xxxm`
        if len(funcs) == 1 and not write_to_xxxm:  
            data_dict['X'] = agg
        else:
            data_dict[write_key]['sum'] = agg
    # here and below for count, if var is present, these can be calculate alongside var
    if 'mean' in funcs and 'var' not in funcs:  
        agg = groupby.mean()
        if len(funcs) == 1 and not write_to_xxxm:
            data_dict['X'] = agg
        else:
            data_dict[write_key]['mean'] = agg
    if 'count' in funcs and 'var' not in funcs:
        obs_var_dict['obs']['count'] = groupby.count()  # count goes in dim df
    if 'var' in funcs:
        agg = groupby.count_mean_var(dof)
        if len(funcs) == 1 and not write_to_xxxm:
            data_dict['X'] = agg['var']
        else:
            data_dict[write_key]['var'] = agg['var']
            if 'mean' in funcs:
                data_dict[write_key]['mean'] = agg['mean']
            if 'count' in funcs:
                obs_var_dict['obs']['count'] = agg['count']
    adata_agg = AnnData(**{**data_dict, **obs_var_dict})
    if dim == 'var':
        return adata_agg.T
    return adata_agg

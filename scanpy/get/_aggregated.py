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
)

from anndata import AnnData, utils
import numpy as np
import pandas as pd
import collections.abc as cabc
from scipy.sparse import coo_matrix, dia_matrix, spmatrix

Array = Union[np.ndarray, spmatrix]


class Aggregate:
    """
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
    _groupby
        `Series` containing values for grouping by.
    _data
        Data matrix for aggregation.
    _weight
        Weights to be used for aggergation.
    _key_set
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
        """
        Compute the mean per feature per group of observations.

        Returns
        -------
            Array of mean.
        """
        A, _ = self._sparse_aggregator(normalize=True)
        return utils.asarray(A * self._data)

    def count_mean_var(self, dof: int = 1) -> dict:
        """
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

    def _extract_indices(self):
        def _filter_indices(key_set, keys, key_index, df_index, weight_value=None):
            keep = [i for i, k in enumerate(keys) if k in set(key_set)]
            if len(keep) == 0:
                raise ValueError("No keys in key_set found in adata.obs[key].")
            elif len(keep) < len(keys):
                mask = np.in1d(key_index, keep)
                remap = np.zeros(len(keys), dtype=np.int64)
                for i, j in enumerate(keep):
                    remap[j] = i
                keys = [keys[j] for j in keep]
                key_index = np.array(
                    [remap[i] for i in key_index[mask]], dtype=np.int64
                )
                df_index = df_index[mask]
                if weight_value is not None:
                    weight_value = weight_value[mask]
            return keys, key_index, df_index, weight_value

        key_value = self._groupby
        keys, key_index = np.unique(_ndarray_from_seq(key_value), return_inverse=True)
        df_index = np.arange(len(key_index))
        if self._weight is None:
            weight_value = None
        else:
            weight_value = self._weight.values[df_index]
        if self._key_set is not None:
            keys, key_index, df_index, weight_value = _filter_indices(
                self._key_set, keys, key_index, df_index, weight_value
            )
        self._key_index = key_index  # passed to count and count_mean_var to avoid re-extracting in the latter
        return keys, key_index, df_index, weight_value


def _power(X, power):
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
    """Find all columns which are a superset of the key column.

    Args:
        df (pd.DataFrame): DataFrame which contains candidate columns.
        groupby_key (str): Key for column of which to find superset of columns.

    Returns:
        List[str]: Superset columns.
    """
    columns = []
    groupy_key_codes = df[groupby_key].astype('category')
    for key in df:
        if key != groupby_key:
            key_codes = df[key].astype('category')
            if all(
                [
                    key_codes[groupy_key_codes == group_key_code].nunique() == 1
                    for group_key_code in groupy_key_codes
                ]
            ):
                columns += [key]
    return columns


def _df_grouped(df: pd.DataFrame, key: str, key_set: List[str]) -> pd.DataFrame:
    """Generate a grouped-by dataframe (no aggregation) by a key with columns that are supersets of the key column

    Args:
        df (pd.DataFrame): DataFrame to be grouped.
        key (str): Column to be grouped on.
        key_set (List[str]): values in the `key` column to keep before groupby.

    Returns:
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
    how: Literal['count', 'mean', 'sum', 'count_mean_var'] = 'count_mean_var',
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

    Parameters
    ----------
        adata:
            :class:`~anndata.AnnData` to be aggregated.
        by:
            Key of the column to be grouped-by.
        how: 
            How to aggregate. Defaults to 'count_mean_var'.
        dim:
            Axis on which to find group by column. Defaults to 'obs'.
        weight_key:
            Key of the `dim` containing weights for a weighted sum aggregation. Defaults to None.
        key_set:
            Subset of dim on which to filter. Defaults to None.
        dof:
            Degrees of freedom for variance. Defaults to 1.
        layer:
            If not None, key for aggregation data. Defaults to None.
        obsm:
            If not None, key for aggregation data. Defaults to None.
        varm:
            If not None, key for aggregation data. Defaults to None.

    Returns
    -------
        AnnData:
            Aggregated :class:`~anndata.AnnData`.
    """
    data = adata.X
    write_to_xxxm = None
    if varm is not None:
        data = adata.varm[varm]
        write_to_xxxm = True  # the data will have to be transposed so this is accurate
    elif obsm is not None:
        data = adata.obsm[obsm]
        write_to_xxxm = True
    elif layer is not None:
        data = adata.layers[layer]
        if dim == 'var':
            data = data.T
    elif (
        dim == 'var'
    ):  # i.e., all of `varm`, `obsm`, `layers` are None so we use `X` which must be transposed
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
        how=how,
        dof=dof,
    )


@aggregated.register(np.ndarray)
@aggregated.register(spmatrix)
def aggregated_from_array(
    data,
    groupby_df: pd.DataFrame,
    dim: str,
    by: str,
    write_to_xxxm: bool,
    no_groupby_df: pd.DataFrame,
    weight_key: Optional[str] = None,
    key_set: Optional[Iterable[str]] = None,
    how: Literal['count', 'mean', 'sum', 'count_mean_var'] = 'count_mean_var',
    dof: int = 1,
) -> AnnData:
    """\
    Aggregate data based on one of the columns of one of a `~pd.DataFrame`.

    Parameters
    ----------
        data:
            Data for aggregation.
        groupby_df:
            `~pd.DataFrame` with column to be grouped on.
        dim:
            Key of AnnData corresponding to the dim on which the grouped by data belongs.
        by:
            Key of the groupby `~pd.DataFrame` for grouping.
        write_to_xxxm:
            Whether or not to write aggregation data to `varm` or `obsm` (based on `dim`)
        no_groupby_df:
            `~pd.DataFrame` on the opposite dim of dim.
        weight_key:
            Key of the `dim` containing weights for a weighted sum aggregation. Defaults to None.
        key_set:
            Defaults to None. Subset of dim on which to filter.
        how:
            How to aggregate. Defaults to 'count_mean_var'.
        dof: 
            Degrees of freedom for variance. Defaults to 1.

    Returns
    -------
        AnnData:
            Aggregated :class:`~anndata.AnnData`.
    """
    groupby = Aggregate(
        groupby=groupby_df[by],
        data=data,
        weight=groupby_df[weight_key] if weight_key is not None else None,
        key_set=key_set,
    )
    # groupby df is put in `obs`, nongroupby in `var` to be transposed later as appropriate
    obs_var_dict = {'obs': _df_grouped(groupby_df, by, key_set), 'var': no_groupby_df}
    data_dict = {}
    if how == 'count':
        obs_var_dict['obs']['count'] = groupby.count()  # count goes in df
    elif how == 'mean':
        agg = groupby.mean()
        data_dict = {'obsm': {'mean': agg}} if write_to_xxxm else {'X': agg}
    elif how == 'sum':
        agg = groupby.sum()
        data_dict = {'obsm': {'sum': agg}} if write_to_xxxm else {'X': agg}
    else:
        agg = groupby.count_mean_var(dof)
        write_key = 'obsm' if write_to_xxxm else 'layers'
        obs_var_dict['obs']['count'] = agg['count']  # count in df
        data_dict = {
            write_key: {'mean': agg['mean'], 'var': agg['var']}
        }  # others in layers/obsm
    adata_agg = AnnData(**{**data_dict, **obs_var_dict})
    if dim == 'var':
        return adata_agg.T
    return adata_agg

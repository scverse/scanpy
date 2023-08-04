from functools import cached_property
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


class GroupBy:
    """
    Functionality for grouping and aggregating AnnData observations by key on `obs`
     and if you want to do `var`, transpose the AnnData object first.

    There is currently support for count, sum, mean, and varience per group, and for scores
    derived from these per pair of groups.

    Set `weight` for weighted sum, mean, and variance.

    Set `key_set` to a list of keys to most efficiently compute results for a subset of groups.

    **Implementation**

    Moments are computed using weighted sum aggregation of AnnData obsevations per variable
    (i.e., feature) via multiplication by a sparse coordinate matrix A, exposed by
    `_sparse_aggregator`. The approach works with data in ndarray or scipy sparse formats, with
    no view or copy overhead on runtime or memory, even when filtering keys.

    Runtime is effectively computation of the product A * X, i.e. the count of (non-zero)
    entries in X with multiplicity the number of group memberships for that entry. This is
    O(data) for partitions (each observation belonging to exactly one group), independent of
    the number of groups.

    To compute scores, first statistics are computed for each group in at least one pair, and
    then scores are computed for each pair using the statistics. Runtime is dominated by the
    former, so is effectively independent of the number of pairs.

    Params
    ------
    adata
    key
        Group key field in adata.obs.
    data
        Element of the AnnData to aggregate (default None yields adata.X).  Should have the same dimensions as the AnnData object.
    weight
        Weight field in adata.obs of type float.
    key_set
        Subset of keys to which to filter.
    write_to_obsm
        Whether to write to `obsm` or not
    """

    _adata: AnnData
    _key: str
    _data: Union[np.ndarray, spmatrix]
    _weight: Optional[str]
    _key_set: AbstractSet[str]
    _key_index: Optional[np.ndarray]  # caution, may be stale if attributes are updated
    _write_to_obsm: bool

    def __init__(
        self,
        adata: AnnData,
        key: str,
        data: Union[np.ndarray, spmatrix],
        *,
        weight: Optional[str] = None,
        key_set: Optional[Iterable[str]] = None,
        write_to_obsm: bool = False,
    ):
        self._adata = adata
        self._data = data
        self._key = key
        self._weight = weight
        self._key_set = None if key_set is None else dict.fromkeys(key_set).keys()
        self._key_index = None
        self._write_to_obsm = write_to_obsm

    @cached_property
    def _superset_columns(self) -> List[str]:
        """Find all columns which are a superset of the key column.

        Returns:
            List[str]: Superset columns.
        """
        columns = []
        groupy_key_codes = self._adata.obs[self._key].astype('category')
        for key in self._adata.obs:
            if key != self._key:
                key_codes = self._adata.obs[key].astype('category')
                if all(
                    [
                        key_codes[groupy_key_codes == group_key_code].nunique() == 1
                        for group_key_code in groupy_key_codes
                    ]
                ):
                    columns += [key]
        return columns

    @cached_property
    def _df_grouped(self) -> pd.DataFrame:
        df = self._adata.obs.copy()
        if self._key_set is not None:
            df = df[df[self._key].isin(self._key_set)]
        if df[self._key].dtype.name == 'category':
            df[self._key] = df[self._key].cat.remove_unused_categories()
        return df.groupby(self._key).first()[self._superset_columns]

    @cached_property
    def _base_axis_indices(self) -> pd.Index:
        return pd.DataFrame(index=pd.Index(self._adata.var_names).copy())

    @cached_property
    def _obs_var_dict(self) -> dict:
        return {
            'obs': self._df_grouped,
            'var': self._base_axis_indices,
        }

    def count(self) -> AnnData:
        """
        Count the number of observations in each group.

        Returns
        -------
            Series of counts indexed by key.
        """
        _, key_index, _, _ = self._extract_indices()
        count_ = np.bincount(key_index)
        obs_var_dict = self._obs_var_dict
        obs_var_dict['obs']['count'] = count_
        return AnnData(**obs_var_dict)

    def sum(self) -> AnnData:
        """
        Compute the sum per feature per group of observations.

        Returns
        -------
            AnnData with sum in X indexed on obs by key with var from adata.
        """
        A, _ = self._sparse_aggregator(normalize=False)
        X = utils.asarray(A * self._data)
        data_dict = (
            {'obsm': {'sum': X}} if self._write_to_obsm else { 'X': X }
        )
        return AnnData(**{**self._obs_var_dict, **data_dict})

    def mean(self) -> AnnData:
        """
        Compute the mean per feature per group of observations.

        Returns
        -------
            AnnData with means in X indexed on obs by key with var from adata.
        """
        A, _ = self._sparse_aggregator(normalize=True)
        X = utils.asarray(A * self._data)
        data_dict = (
            {'obsm': {'mean': X}} if self._write_to_obsm else { 'X': X }
        )
        return AnnData(**{**self._obs_var_dict, **data_dict})

    def count_mean_var(self, dof: int = 1) -> AnnData:
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
            AnnData with mean and var in layers indexed on obs by key with var from adata.  Counts are in obs under counts.
        """
        assert dof >= 0
        A, _ = self._sparse_aggregator(normalize=True)
        count_ = np.bincount(self._key_index)
        mean_ = utils.asarray(A @ self._data)
        mean_sq = utils.asarray(A @ _power(self._data, 2))
        if self._weight is None:
            sq_mean = mean_**2
        else:
            A_unweighted, _ = GroupBy(
                adata=self._adata,
                data=self._data,
                key=self._key,
                key_set=self._key_set,
                write_to_obsm=self._write_to_obsm,
            )._sparse_aggregator()
            mean_unweighted = utils.asarray(A_unweighted * self._data)
            sq_mean = 2 * mean_ * mean_unweighted + mean_unweighted**2
        var_ = mean_sq - sq_mean
        precision = 2 << (42 if self._data.dtype == np.float64 else 20)
        # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
        var_[precision * var_ < sq_mean] = 0
        if dof != 0:
            var_ *= (count_ / (count_ - dof))[:, np.newaxis]
        obs_var_dict = self._obs_var_dict
        obs_var_dict['obs']['count'] = count_
        write_to_obsm = 'obsm' if self._write_to_obsm else 'layers'
        return AnnData(
            **{
                **obs_var_dict,
                write_to_obsm: {
                    'mean': mean_,
                    'var': var_,
                },
            }
        )

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

        key_value = self._adata.obs[self._key]
        keys, key_index = np.unique(_ndarray_from_seq(key_value), return_inverse=True)
        df_index = np.arange(len(key_index))
        if self._weight is None:
            weight_value = None
        else:
            weight_value = self._adata.obs[self._weight].values[df_index]
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


def aggregated(
    adata: AnnData,
    by: str,
    how: Literal['count', 'mean', 'sum', 'count_mean_var'],
    df_key: Literal['obs', 'var'] = 'obs',
    weight: Optional[str] = None,
    key_set: Optional[Iterable[str]] = None,
    dof: int = 1,
    layer=None,
    obsm=None,
    varm=None,
):
    data = adata.X
    write_to_obsm = None
    if varm is not None:
        data = adata.varm[varm]
        write_to_obsm = True  # the data will have to be transposed so this is accurate
    elif obsm is not None:
        data = adata.obsm[obsm]
        write_to_obsm = True
    elif layer is not None:
        data = adata.layers[layer]
        if df_key == 'var':
            data = data.T
    elif df_key == 'var':
        data = data.T
    groupby = GroupBy(
        adata=adata if df_key == 'obs' else adata.T,
        data=data,
        key=by,
        weight=weight,
        key_set=key_set,
        write_to_obsm=write_to_obsm,
    )
    if how == 'count':
        data = groupby.count()
    elif how == 'mean':
        data = groupby.mean()
    elif how == 'sum':
        data = groupby.sum()
    else:
        data = groupby.count_mean_var(dof)
    if df_key == 'var':
        return data.T
    return data

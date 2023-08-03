from collections import defaultdict
from typing import (
    Optional,
    Iterable,
    AbstractSet,
    Sequence,
    Tuple,
    Union,
    NamedTuple,
    Literal,
)

from anndata import AnnData, utils
import numpy as np
import pandas as pd
import collections.abc as cabc
from scipy.sparse import coo_matrix, dia_matrix, spmatrix


class CountMeanVar(NamedTuple):
    count: Optional[pd.Series] = None
    mean: pd.DataFrame = None
    var: Optional[pd.DataFrame] = None

    @classmethod
    def from_df(cls, df: pd.DataFrame) -> "CountMeanVar":
        return CountMeanVar(
            count=df["count"] if "count" in df.columns else None,
            mean=df["mean"],
            var=df["var"] if "var" in df.columns else None,
        )

    def map(self, f=lambda v: v, **fs) -> "CountMeanVar":
        fs = defaultdict(lambda: f, fs)
        return CountMeanVar(
            **{
                stat: fs[stat](val)
                for stat, val in self._asdict().items()
                if val is not None
            }
        )


Score = Literal[
    "diff-score", "fold-score", "t-score", "v-score", "t-score-pooled", "v-score-pooled"
]


class GroupBy:
    """
    Functionality for grouping and aggregating AnnData observations by key, per variable.

    There is currently support for count, sum, mean, and varience per group, and for scores
    derived from these per pair of groups.

    Set `weight` for weighted sum, mean, and variance.

    Set `explode` to True and use a key of type tuple to assign observations to multiple groups.
    In this case, repetition of a key confers multiplicity of the observation in the group.

    Set `key_set` to a list of keys to most efficiently compute results for a subset of groups.

    NaN values propagate, with the exception that `score_pairs` sets non-finite scores to 0 by
    default. Use the pd_* methods to instead mask NaN values. These slower methods convert data
    to dense format and do not currently support weight, explode, or key_set.

    **Implementation**

    Moments are computed using weighted sum aggregation of AnnData obsevations per variable
    (i.e., feature) via multiplication by a sparse coordinate matrix A, exposed by
    `sparse_aggregator`. The approach works with data in ndarray or scipy sparse formats, with
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
        Element of the AnnData to aggregate (default None yields adata.X)
    weight
        Weight field in adata.obs of type float.
    explode
        If False, each observation is assigned to the group keyed by adata.obs[key].
        If True, each observation is assigned to all groups in tuple adata.obs[key].
    key_set
        Subset of keys to which to filter.
    """

    adata: AnnData
    key: str
    data: Union[np.ndarray, spmatrix]
    weight: Optional[str]
    explode: bool
    key_set: AbstractSet[str]
    _key_index: Optional[np.ndarray]  # caution, may be stale if attributes are updated

    def __init__(
        self,
        adata: AnnData,
        key: str,
        *,
        data: Union[np.ndarray, spmatrix, None] = None,
        weight: Optional[str] = None,
        explode: bool = False,
        key_set: Optional[Iterable[str]] = None,
    ):
        self.adata = adata
        self.data = adata.X if data is None else data
        self.key = key
        self.weight = weight
        self.explode = explode
        self.key_set = None if key_set is None else dict.fromkeys(key_set).keys()
        self._key_index = None

    def count(self) -> pd.Series:
        """
        Count the number of observations in each group.

        Returns
        -------
            Series of counts indexed by key.
        """
        keys, key_index, _, _ = self._extract_indices()
        count_ = np.bincount(key_index)
        return pd.Series(
            data=count_,
            index=pd.Index(keys, name=self.key, tupleize_cols=False),
            name="count",
        )

    def sum(self) -> pd.DataFrame:
        """
        Compute the sum per feature per group of observations.

        Returns
        -------
        DataFrame of sums indexed by key with columns from adata.
        """
        A, keys = self.sparse_aggregator(normalize=False)
        return pd.DataFrame(
            index=pd.Index(keys, name=self.key, tupleize_cols=False),
            columns=self.adata.var_names.copy(),
            data=utils.asarray(A * self.data),
        )

    def mean(self) -> pd.DataFrame:
        """
        Compute the mean per feature per group of observations.

        Returns
        -------
            DataFrame of means indexed by key with columns from adata.
        """
        A, keys = self.sparse_aggregator(normalize=True)
        return pd.DataFrame(
            index=pd.Index(keys, name=self.key, tupleize_cols=False),
            columns=self.adata.var_names.copy(),
            data=utils.asarray(A * self.data),
        )

    def var(self, dof: int = 1) -> pd.DataFrame:
        """
        Compute the variance per feature per group of observations.

        See also count_mean_var, which is more efficient when the mean is also desired.

        Params
        ------
        dof
            Degrees of freedom for variance.

        Returns
        -------
            DataFrame of variances indexed by key with columns from adata.
        """
        return self.count_mean_var(dof).var

    def count_mean_var(self, dof: int = 1) -> CountMeanVar:
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
            Dictionary with keys (count, mean, var) and values the corresponding Series and DataFrames.
        """
        assert dof >= 0
        A, keys = self.sparse_aggregator(normalize=True)
        count_ = np.bincount(self._key_index)
        mean_ = utils.asarray(A @ self.data)
        mean_sq = utils.asarray(A @ _power(self.data, 2))
        if self.weight is None:
            sq_mean = mean_ ** 2
        else:
            A_unweighted, _ = GroupBy(
                self.adata, self.key, explode=self.explode, key_set=self.key_set
            ).sparse_aggregator()
            mean_unweighted = utils.asarray(A_unweighted * self.data)
            sq_mean = 2 * mean_ * mean_unweighted + mean_unweighted ** 2
        var_ = mean_sq - sq_mean
        precision = 2 << (42 if self.data.dtype == np.float64 else 20)
        # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
        var_[precision * var_ < sq_mean] = 0
        if dof != 0:
            var_ *= (count_ / (count_ - dof))[:, np.newaxis]

        index = pd.Index(keys, name=self.key, tupleize_cols=False)
        count_sr = pd.Series(index=index, data=count_, name="count")
        mean_df = pd.DataFrame(
            index=index.copy(), columns=self.adata.var_names.copy(), data=mean_
        )
        var_df = pd.DataFrame(
            index=index.copy(), columns=self.adata.var_names.copy(), data=var_
        )
        return CountMeanVar(count=count_sr, mean=mean_df, var=var_df)

    def sparse_aggregator(
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
        keys, key_index, obs_index, weight_value = self._extract_indices()
        if obs_index is None:
            obs_index = np.arange(len(key_index))
        if self.weight is None:
            weight_value = np.ones(len(key_index))
        A = coo_matrix(
            (weight_value, (key_index, obs_index)),
            shape=(len(keys), self.data.shape[0]),
        )
        if normalize:
            n_row = A.shape[0]
            row_sums = np.asarray(A.sum(axis=1))
            D = dia_matrix(((row_sums.T ** -1), [0]), shape=(n_row, n_row))
            A = D * A
        return A, keys

    def _extract_indices(self):
        def _filter_indices(key_set, keys, key_index, obs_index, weight_value=None):
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
                obs_index = obs_index[mask]
                if weight_value is not None:
                    weight_value = weight_value[mask]
            return keys, key_index, obs_index, weight_value

        key_value = self.adata.obs[self.key]
        if self.explode:
            assert isinstance(
                key_value.iloc[0], tuple
            ), "key type must be tuple to explode"
            keys, key_index = np.unique(
                _ndarray_from_seq([k for ks in key_value for k in ks]),
                return_inverse=True,
            )
            obs_index = np.array([i for i, ks in enumerate(key_value) for _ in ks])
        else:
            keys, key_index = np.unique(
                _ndarray_from_seq(key_value), return_inverse=True
            )
            obs_index = np.arange(len(key_index))
        if self.weight is None:
            weight_value = None
        else:
            weight_value = self.adata.obs[self.weight].values[obs_index]
        if self.key_set is not None:
            keys, key_index, obs_index, weight_value = _filter_indices(
                self.key_set, keys, key_index, obs_index, weight_value
            )
        self._key_index = key_index  # passed to count and count_mean_var to avoid re-extracting in the latter
        return keys, key_index, obs_index, weight_value

    def pd_mean(self) -> pd.DataFrame:
        """
        Slower implementation of mean that masks NaN values.
        """
        assert (
            (self.weight is None) and (self.explode is False) and (self.key_set is None)
        )
        df = pd.DataFrame(
            index=self.adata.obs[self.key],
            columns=self.adata.var_names,
            data=utils.asarray(self.data),
        )
        return df.groupby(self.key).mean()

    def _pd_count_mean_var_df(self) -> pd.DataFrame:
        assert (
            (self.weight is None) and (self.explode is False) and (self.key_set is None)
        )
        aggs = ["count", "mean", "var"]
        df = pd.DataFrame(
            index=self.adata.obs[self.key],
            columns=self.adata.var_names,
            data=utils.asarray(self.data),
        )
        return df.groupby(self.key).agg(aggs).swaplevel(axis=1).sort_index(axis=1)

    def pd_count_mean_var(self) -> CountMeanVar:
        """
        Slower implementation of count_mean_var that masks NaN values.
        """
        return CountMeanVar.from_df(self._pd_count_mean_var_df())


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

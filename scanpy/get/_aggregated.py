from __future__ import annotations

from collections.abc import Iterable, Sequence, Set
from functools import singledispatch
from typing import TYPE_CHECKING, Literal, NamedTuple, get_args
from typing import Union as _U

import numpy as np
import pandas as pd
from anndata import AnnData, utils
from scipy import sparse

if TYPE_CHECKING:
    from numpy.typing import NDArray

Array = _U[np.ndarray, sparse.spmatrix]
AggType = Literal["count", "mean", "sum", "var"]


class CMV(NamedTuple):
    count: NDArray[np.integer]
    mean: NDArray[np.floating]
    var: NDArray[np.floating]


class Indices(NamedTuple):
    keys: np.ndarray
    key_index: np.ndarray
    df_index: np.ndarray
    weight_value: pd.Series | Array | None


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
        Weights to be used for aggregation.
    """

    def __init__(self, groupby, data, weight=None):
        self.groupby = groupby
        self.indicator_matrix = sparse_indicator(groupby)
        self.data = data
        self.weight = weight

    groupby: pd.Series
    data: Array
    weight: pd.Series | Array
    key_set: Set[str] | None

    def count(self) -> np.ndarray:
        """\
        Count the number of observations in each group.

        Returns
        -------
        Array of counts.
        """
        # pattern = self.data._with_data(np.broadcast_to(1, len(self.data.data)))
        # return self.indicator_matrix @ pattern
        return self.indicator_matrix @ (self.data != 0)

    def sum(self) -> Array:
        """\
        Compute the sum per feature per group of observations.

        Returns
        -------
        Array of sum.
        """
        return utils.asarray(self.indicator_matrix @ self.data)

    def mean(self) -> Array:
        """\
        Compute the mean per feature per group of observations.

        Returns
        -------
        Array of mean.
        """
        return (
            utils.asarray(self.indicator_matrix @ self.data)
            / np.bincount(self.groupby.codes)[:, None]
        )

    def count_mean_var(self, dof: int = 1, *, _indices: Indices | None = None) -> CMV:
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
        Object with `count`, `mean`, and `var` attributes.
        """
        assert dof >= 0
        count_ = self.count()
        group_counts = np.bincount(self.groupby.codes)
        mean_ = self.mean()
        # sparse matrices do not support ** for elementwise power.
        mean_sq = (
            utils.asarray(self.indicator_matrix @ _power(self.data, 2))
            / group_counts[:, None]
        )
        if self.weight is None:
            sq_mean = mean_**2
        else:
            A_unweighted = sparse_indicator(self.groupby)
            # , _ = Aggregate(
            #     groupby=self.groupby,
            #     data=self.data,
            #     weight=self.weight,  # TODO: why pass weights when creating unweighted A?
            #     key_set=self.key_set,
            # )._sparse_aggregator()
            mean_unweighted = utils.asarray(A_unweighted @ self.data)
            sq_mean = 2 * mean_ * mean_unweighted + mean_unweighted**2
        var_ = mean_sq - sq_mean
        # TODO: Why these values exactly? Because they are high relative to the datatype?
        # (unchanged from original code: https://github.com/scverse/anndata/pull/564)
        precision = 2 << (42 if self.data.dtype == np.float64 else 20)
        # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
        var_[precision * var_ < sq_mean] = 0
        if dof != 0:
            var_ *= (group_counts / (group_counts - dof))[:, np.newaxis]
        return CMV(count=count_, mean=mean_, var=var_)


# def count_mean_var_spd(by, data):
#     sums = np.zeros((by.shape[0],data.shape[1]))
#     counts = np.zeros((by.shape[0],data.shape[1]))
#     sums = by.toarray() @ data
#     counts = by.toarray() @ data._with_data(np.ones(len(data.data),dtype=data.data.dtype))
#     n_cells = np.array(by.sum(axis= 1).astype(data.dtype))
#     means = sums/n_cells
#     sq_mean = by.toarray() @ data.multiply(data)/n_cells
#     var = sq_mean - np.power(means, 2)
#     var *= n_cells / (n_cells - 1)
#     return sums, counts, means, var


def _power(X: Array, power: float | int) -> Array:
    """\
    Generate elementwise power of a matrix.

    Needed for non-square sparse matrices because they do not support ** so the `.power` function is used.

    Params
    ------
    X
        Matrix whose power is to be raised.
    power
        Integer power value

    Returns
    -------
    Matrix whose power has been raised.
    """
    return X**power if isinstance(X, np.ndarray) else X.power(power)


def _ndarray_from_seq(lst: Sequence):
    # prevents expansion of iterables as axis
    n = len(lst)
    if n > 0 and isinstance(lst[0], Iterable):
        arr = np.empty(n, dtype=object)
        arr[:] = lst
    else:
        arr = np.array(lst)
    return arr


@singledispatch
def aggregated(
    adata: AnnData,
    by: str,
    func: AggType | Iterable[AggType],
    *,
    dim: Literal["obs", "var"] = "obs",
    weight_key: str | None = None,
    dof: int = 1,
    layer: str | None = None,
    obsm: str | None = None,
    varm: str | None = None,
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
    # TODO replace with get helper
    if sum(p is not None for p in [varm, obsm, layer]) > 1:
        raise TypeError("Please only provide one (or none) of varm, obsm, or layer")
    if varm is not None:
        data = adata.varm[varm]
    elif obsm is not None:
        data = adata.obsm[obsm]
    elif layer is not None:
        data = adata.layers[layer]
        if dim == "var":
            data = data.T
    elif dim == "var":
        # i.e., all of `varm`, `obsm`, `layers` are None so we use `X` which must be transposed
        data = data.T
    return aggregated(
        data,
        groupby_df=getattr(adata, dim),
        dim=dim,
        by=by,
        # write_to_xxxm=write_to_xxxm,
        no_groupby_df=getattr(adata, "var" if dim == "obs" else "obs"),
        weight_key=weight_key,
        # key_set=key_set,
        func=func,
        dof=dof,
    )


@aggregated.register(np.ndarray)
@aggregated.register(sparse.spmatrix)
def aggregated_from_array(
    data,
    groupby_df: pd.DataFrame,
    func: AggType | Iterable[AggType],
    dim: str,
    by: str,
    no_groupby_df: pd.DataFrame,
    weight_key: str | None = None,
    dof: int = 1,
) -> AnnData:
    """Aggregate data based on one of the columns of one of a `~pd.DataFrame`."""
    categorical = _combine_categories(groupby_df, by)
    groupby = Aggregate(
        groupby=categorical,
        data=data,
        weight=groupby_df[weight_key] if weight_key is not None else None,
    )
    # groupby df is put in `obs`, nongroupby in `var` to be transposed later as appropriate
    adata_kw = dict(
        X=None,
        layers={},
        obs=pd.DataFrame(index=categorical.categories),
        var=no_groupby_df,
        obsm={},
    )
    funcs = set([func] if isinstance(func, str) else func)
    if unknown := funcs - set(get_args(AggType)):
        raise ValueError(f"func {unknown} is not one of {get_args(AggType)}")
    if "sum" in funcs:  # sum is calculated separately from the rest
        agg = groupby.sum()
        adata_kw["layers"]["sum"] = agg
    # here and below for count, if var is present, these can be calculate alongside var
    if "mean" in funcs and "var" not in funcs:
        agg = groupby.mean()
        adata_kw["layers"]["mean"] = agg
    if "count" in funcs and "var" not in funcs:
        adata_kw["layers"]["count"] = groupby.count()  # count goes in dim df
    if "var" in funcs:
        aggs = groupby.count_mean_var(dof)
        adata_kw["layers"]["var"] = aggs.var
        if "mean" in funcs:
            adata_kw["layers"]["mean"] = aggs.mean
        if "count" in funcs:
            adata_kw["layers"]["count"] = aggs.count

    adata_agg = AnnData(**adata_kw)
    if dim == "var":
        return adata_agg.T
    return adata_agg


def _combine_categories(label_df: pd.DataFrame, cols: list[str]) -> pd.Categorical:
    from itertools import product

    if isinstance(cols, str):
        cols = [cols]

    df = pd.DataFrame(
        {c: pd.Categorical(label_df[c]).remove_unused_categories() for c in cols},
    )
    result_categories = [
        "_".join(map(str, x)) for x in product(*[df[c].cat.categories for c in cols])
    ]
    n_categories = [len(df[c].cat.categories) for c in cols]

    factors = np.ones(len(cols) + 1, dtype=np.int32)  # First factor needs to be 1
    np.cumsum(n_categories[::-1], out=factors[1:])
    factors = factors[:-1][::-1]

    # TODO: pick a more optimal bit width
    final_codes = np.zeros(df.shape[0], dtype=np.int32)
    for factor, c in zip(factors, cols):
        final_codes += df[c].cat.codes * factor

    return pd.Categorical.from_codes(
        final_codes, categories=result_categories
    ).remove_unused_categories()


def sparse_indicator(
    categorical, weights: None | np.ndarray = None
) -> sparse.coo_matrix:
    if weights is None:
        weights = np.broadcast_to(1.0, len(categorical))
    A = sparse.coo_matrix(
        (weights, (categorical.codes, np.arange(len(categorical)))),
        shape=(len(categorical.categories), len(categorical)),
    )
    return A

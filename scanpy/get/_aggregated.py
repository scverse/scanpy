from __future__ import annotations

from collections.abc import Iterable, Sequence, Set
from functools import singledispatch
from typing import Literal, NamedTuple, get_args
from typing import Union as _U

import numpy as np
import pandas as pd
from anndata import AnnData, utils
from scipy import sparse

Array = _U[np.ndarray, sparse.spmatrix]
AggType = Literal["count_nonzero", "mean", "sum", "var"]


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
        self.indicator_matrix = sparse_indicator(groupby, weight=weight)
        self.data = data
        self.weight = weight

    groupby: pd.Series
    data: Array
    weight: pd.Series | Array
    key_set: Set[str] | None

    def count_nonzero(self) -> np.ndarray:
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

    def mean_var(self, dof: int = 1) -> tuple[np.ndarray, np.ndarray]:
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
        return mean_, var_


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
def aggregate(
    adata: AnnData,
    by: str | list[str],
    func: AggType | Iterable[AggType],
    *,
    dim: Literal["obs", "var"] = "obs",
    weight: str | None = None,
    dof: int = 1,
    layer: str | None = None,
    obsm: str | None = None,
    varm: str | None = None,
) -> AnnData:
    """\
    Aggregate data matrix based on some categorical grouping.

    This function is useful for pseudobulking as well as plotting.

    Aggregation to perform is specified by `func`, which can be a single metric or a
    list of metrics. Each metric is computed over the group and results in a new layer
    in the output `AnnData` object.

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
    weight
        Key of the `dim` containing weights for a weighted sum aggregation.
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

    Examples
    --------

    Calculating mean expression and number of nonzero entries per cluster:

    >>> import scanpy as sc, pandas as pd
    >>> pbmc = sc.datasets.pbmc3k_processed().raw.to_adata()
    >>> pbmc.shape
    (2638, 13714)
    >>> aggregated = sc.get.aggregate(pbmc, by="louvain", func=["mean", "count_nonzero"])
    >>> aggregated
    AnnData object with n_obs × n_vars = 8 × 13714
        obs: 'louvain'
        var: 'n_cells'
        layers: 'mean', 'count_nonzero'

    We can group over multiple columns:

    >>> pbmc.obs["percent_mito_binned"] = pd.cut(pbmc.obs["percent_mito"], bins=5)
    >>> sc.get.aggregate(pbmc, by=["louvain", "percent_mito_binned"], func=["mean", "count_nonzero"])
    AnnData object with n_obs × n_vars = 40 × 13714
        obs: 'louvain', 'percent_mito_binned'
        var: 'n_cells'
        layers: 'mean', 'count_nonzero'

    Note that this filters out any combination of groups that wasn't present in the original data.
    """
    if dim not in ["obs", "var"]:
        raise ValueError(f"dim must be one of 'obs' or 'var', was '{dim}'")
    # TODO replace with get helper
    data = adata.X
    if sum(p is not None for p in [varm, obsm, layer]) > 1:
        raise TypeError("Please only provide one (or none) of varm, obsm, or layer")
    if varm is not None:
        if dim != "var":
            raise ValueError("varm can only be used when dim is 'var'")
        data = adata.varm[varm]
    elif obsm is not None:
        if dim != "obs":
            raise ValueError("obsm can only be used when dim is 'obs'")
        data = adata.obsm[obsm]
    elif layer is not None:
        data = adata.layers[layer]
        if dim == "var":
            data = data.T
    elif dim == "var":
        # i.e., all of `varm`, `obsm`, `layers` are None so we use `X` which must be transposed
        data = data.T

    dim_df = getattr(adata, dim)
    categorical, new_label_df = _combine_categories(dim_df, by)
    if isinstance(weight, str):
        weight = dim_df[weight]
    # Actual computation
    layers = aggregate(
        data,
        by=categorical,
        func=func,
        dof=dof,
        weight=weight,
    )
    result = AnnData(
        layers=layers,
        obs=new_label_df,
        var=getattr(adata, "var" if dim == "obs" else "obs"),
    )

    # result = aggregate(
    #     data,
    #     groupby_df=getattr(adata, dim),
    #     by=by,
    #     # write_to_xxxm=write_to_xxxm,
    #     no_groupby_df=getattr(adata, "var" if dim == "obs" else "obs"),
    #     weight_key=weight_key,
    #     # key_set=key_set,
    #     func=func,
    #     dof=dof,
    # )

    if dim == "var":
        return result.T
    else:
        return result


@aggregate.register(np.ndarray)
@aggregate.register(sparse.spmatrix)
def aggregate_array(
    data,
    by: pd.Categorical,
    func: AggType | Iterable[AggType],
    *,
    dof: int = 1,
    weight: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    groupby = Aggregate(groupby=by, data=data, weight=weight)
    result = {}

    funcs = set([func] if isinstance(func, str) else func)
    if unknown := funcs - set(get_args(AggType)):
        raise ValueError(f"func {unknown} is not one of {get_args(AggType)}")

    if "sum" in funcs:  # sum is calculated separately from the rest
        agg = groupby.sum()
        result["sum"] = agg
    # here and below for count, if var is present, these can be calculate alongside var
    if "mean" in funcs and "var" not in funcs:
        agg = groupby.mean()
        result["mean"] = agg
    if "count_nonzero" in funcs:
        result["count_nonzero"] = groupby.count_nonzero()
    if "var" in funcs:
        mean_, var_ = groupby.mean_var(dof)
        result["var"] = var_
        if "mean" in funcs:
            result["mean"] = mean_

    return result


# @aggregate.register(np.ndarray)
# @aggregate.register(sparse.spmatrix)
# def aggregate_from_array(
#     data,
#     groupby_df: pd.DataFrame,
#     func: AggType | Iterable[AggType],
#     by: str,
#     no_groupby_df: pd.DataFrame,
#     weight_key: str | None = None,
#     dof: int = 1,
# ) -> AnnData:
#     """Aggregate data based on one of the columns of one of a `~pd.DataFrame`."""
#     categorical, new_label_df = _combine_categories(groupby_df, by)
#     groupby = Aggregate(
#         groupby=categorical,
#         data=data,
#         weight=groupby_df[weight_key] if weight_key is not None else None,
#     )
#     # groupby df is put in `obs`, nongroupby in `var` to be transposed later as appropriate
#     adata_kw = dict(
#         X=None,
#         layers={},
#         obs=new_label_df,
#         var=no_groupby_df,
#         obsm={},
#     )
#     funcs = set([func] if isinstance(func, str) else func)
#     if unknown := funcs - set(get_args(AggType)):
#         raise ValueError(f"func {unknown} is not one of {get_args(AggType)}")
#     if "sum" in funcs:  # sum is calculated separately from the rest
#         agg = groupby.sum()
#         adata_kw["layers"]["sum"] = agg
#     # here and below for count, if var is present, these can be calculate alongside var
#     if "mean" in funcs and "var" not in funcs:
#         agg = groupby.mean()
#         adata_kw["layers"]["mean"] = agg
#     if "count_nonzero" in funcs:
#         adata_kw["layers"]["count_nonzero"] = groupby.count_nonzero()
#     if "var" in funcs:
#         mean_, var_ = groupby.mean_var(dof)
#         adata_kw["layers"]["var"] = var_
#         if "mean" in funcs:
#             adata_kw["layers"]["mean"] = mean_

#     adata_agg = AnnData(**adata_kw)
#     return adata_agg


def _combine_categories(
    label_df: pd.DataFrame, cols: list[str]
) -> tuple[pd.Categorical, pd.DataFrame]:
    """
    Returns both the result categories and a dataframe labelling each row
    """
    from itertools import product

    if isinstance(cols, str):
        cols = [cols]

    df = pd.DataFrame(
        {c: pd.Categorical(label_df[c]).remove_unused_categories() for c in cols},
    )
    n_categories = [len(df[c].cat.categories) for c in cols]

    # It's like np.concatenate([x for x in product(*[range(n) for n in n_categories])])
    code_combinations = np.indices(n_categories).reshape(len(n_categories), -1)
    result_categories = [
        "_".join(map(str, x)) for x in product(*[df[c].cat.categories for c in cols])
    ]

    # Dataframe with unique combination of categories for each row
    new_label_df = pd.DataFrame(
        {
            c: pd.Categorical.from_codes(code_combinations[i], df[c].cat.categories)
            for i, c in enumerate(cols)
        },
        index=result_categories,
    )

    # Calculating result codes
    factors = np.ones(len(cols) + 1, dtype=np.int32)  # First factor needs to be 1
    np.cumsum(n_categories[::-1], out=factors[1:])
    factors = factors[:-1][::-1]

    code_array = np.zeros((len(cols), df.shape[0]), dtype=np.int32)
    for i, c in enumerate(cols):
        code_array[i] = df[c].cat.codes
    code_array *= factors[:, None]

    result_categorical = pd.Categorical.from_codes(
        code_array.sum(axis=0), categories=result_categories
    )

    # Filter unused categories
    result_categorical = result_categorical.remove_unused_categories()
    new_label_df = new_label_df.loc[result_categorical.categories]

    return result_categorical, new_label_df


def sparse_indicator(
    categorical, weight: None | np.ndarray = None
) -> sparse.coo_matrix:
    if weight is None:
        weight = np.broadcast_to(1, len(categorical))
    A = sparse.coo_matrix(
        (weight, (categorical.codes, np.arange(len(categorical)))),
        shape=(len(categorical.categories), len(categorical)),
    )
    return A

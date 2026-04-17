from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING, Literal, TypedDict

import numpy as np
import pandas as pd
from anndata import AnnData
from fast_array_utils.stats._power import power as fau_power  # TODO: upstream
from scipy import sparse
from sklearn.utils.sparsefuncs import csc_median_axis_0

from scanpy._compat import CSBase, CSRBase, DaskArray

from .._utils import _resolve_axis, get_literal_vals
from ._kernels import agg_sum_csc, agg_sum_csr, mean_var_csc, mean_var_csr
from .get import _check_mask

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable

    from numpy.typing import NDArray

type Array = np.ndarray | CSBase | DaskArray
type ConstantDtypeAgg = Literal["count_nonzero", "sum", "median"]
type AggType = ConstantDtypeAgg | Literal["mean", "var"]


class Aggregate[ArrayT: np.ndarray | CSBase]:
    """Functionality for generic grouping and aggregating.

    There is currently support for count_nonzero, sum, mean, and variance.

    **Implementation**

    Moments are computed using weighted sum aggregation of data by some feature
    via multiplication by a sparse coordinate matrix A.

    Runtime is effectively computation of the product `A @ X`, i.e. the count of (non-zero)
    entries in X with multiplicity the number of group memberships for that entry.
    This is `O(data)` for partitions (each observation belonging to exactly one group),
    independent of the number of groups.

    Params
    ------
    groupby
        :class:`~pandas.Categorical` containing values for grouping by.
    data
        Data matrix for aggregation.
    mask
        Mask to be used for aggregation.
    """

    def __init__(
        self,
        groupby: pd.Categorical,
        data: ArrayT,
        *,
        mask: NDArray[np.bool] | None = None,
    ) -> None:
        self.groupby = groupby
        if (missing := groupby.isna()).any():
            mask = mask & ~missing if mask is not None else ~missing
        self.indicator_matrix = sparse_indicator(groupby, mask=mask)
        if isinstance(data, CSBase):
            # TODO: Look into if this can be CSR and fast for dense
            self.indicator_matrix = self.indicator_matrix.tocsr()
        self.data = data

    groupby: pd.Categorical
    indicator_matrix: CSRBase | sparse.coo_array
    data: ArrayT

    def count_nonzero(self) -> NDArray[np.integer]:
        """Count the number of observations in each group.

        Returns
        -------
        Array of counts.

        """
        data = self.data
        if isinstance(data, CSBase):
            data = type(data)(
                (np.ones(data.nnz, dtype="uint8"), data.indices, data.indptr),
                shape=data.shape,
            )
        else:
            data = (data != 0).astype("uint8")
        return self._sum(data=data)

    def _sum(self, data: ArrayT):
        if isinstance(data, np.ndarray):
            res = self.indicator_matrix @ data
            if isinstance(res, CSBase):
                return res.toarray()
            return res
        dtype = np.int64 if np.issubdtype(data.dtype, np.integer) else np.float64
        out = np.zeros((self.indicator_matrix.shape[0], data.shape[1]), dtype=dtype)
        (agg_sum_csr if isinstance(data, CSRBase) else agg_sum_csc)(
            self.indicator_matrix, data, out
        )
        if isinstance(data, CSBase):
            nnz = np.count_nonzero(out)
            if nnz / out.size < 0.5:  # heuristic for when to return sparse vs dense
                return type(data)(out)  # convert to sparse type of input
        return out

    def sum(self) -> np.ndarray:
        """Compute the sum per feature per group of observations.

        Returns
        -------
        Array of sum.

        """
        return self._sum(self.data)

    def mean(self) -> Array:
        """Compute the mean per feature per group of observations.

        Returns
        -------
        Array of mean.

        """
        return self.sum() / np.bincount(self.groupby.codes)[:, None]

    def mean_var(self, dof: int = 1) -> tuple[np.ndarray, np.ndarray]:
        """Compute the count, as well as mean and variance per feature, per group of observations.

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
        if isinstance(self.data, np.ndarray):
            mean_ = self.mean()
            # sparse matrices do not support ** for elementwise power.
            mean_sq = self._sum(_power(self.data, 2)) / group_counts[:, None]
            sq_mean = mean_**2
            var_ = mean_sq - sq_mean
        else:
            mean_, var_ = (
                mean_var_csr if isinstance(self.data, CSRBase) else mean_var_csc
            )(self.indicator_matrix, self.data)
            sq_mean = mean_**2
        # TODO: Why these values exactly? Because they are high relative to the datatype?
        # (unchanged from original code: https://github.com/scverse/anndata/pull/564)
        precision = 2 << (42 if self.data.dtype == np.float64 else 20)
        # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
        var_[precision * var_ < sq_mean] = 0
        if dof != 0:
            var_ *= (group_counts / (group_counts - dof))[:, np.newaxis]
        return mean_, var_

    def median(self) -> Array:
        """Compute the median per feature per group of observations.

        Returns
        -------
        Array of median.

        """
        medians = []
        for group in np.unique(self.groupby.codes):
            group_mask = self.groupby.codes == group
            group_data = self.data[group_mask]
            if isinstance(group_data, CSBase):
                if group_data.format != "csc":
                    group_data = group_data.tocsc()
                medians.append(csc_median_axis_0(group_data))
            else:
                medians.append(np.median(group_data, axis=0))
        return np.array(medians)


def _power(x: Array, power: float) -> Array:
    """Generate elementwise power of a matrix.

    Needed for non-square sparse matrices because they do not support `**` so the `.power` function is used.

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
    return x**power if isinstance(x, np.ndarray) else x.power(power)


def aggregate(  # noqa: PLR0912
    adata: AnnData,
    by: str | Collection[str],
    func: AggType | Iterable[AggType],
    *,
    axis: Literal["obs", 0, "var", 1] | None = None,
    mask: NDArray[np.bool] | str | None = None,
    dof: int = 1,
    layer: str | None = None,
    obsm: str | None = None,
    varm: str | None = None,
) -> AnnData:
    """Aggregate data matrix based on some categorical grouping.

    This function is useful for pseudobulking as well as plotting.

    Aggregation to perform is specified by `func`, which can be a single metric or a
    list of metrics. Each metric is computed over the group and results in a new layer
    in the output `AnnData` object.

    If none of `layer`, `obsm`, or `varm` are passed in, `X` will be used for aggregation data.

    .. array-support:: get.aggregate

    Params
    ------
    adata
        :class:`~anndata.AnnData` to be aggregated.
    by
        Key of the column to be grouped-by.
    func
        How to aggregate.
    axis
        Axis on which to find group by column.
    mask
        Boolean mask (or key to column containing mask) to apply along the axis.
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
    >>> aggregated = sc.get.aggregate(
    ...     pbmc, by="louvain", func=["mean", "count_nonzero"]
    ... )
    >>> aggregated
    AnnData object with n_obs × n_vars = 8 × 13714
        obs: 'louvain', 'n_obs_aggregated'
        var: 'n_cells'
        layers: 'mean', 'count_nonzero'

    We can group over multiple columns:

    >>> pbmc.obs["percent_mito_binned"] = pd.cut(pbmc.obs["percent_mito"], bins=5)
    >>> sc.get.aggregate(
    ...     pbmc, by=["louvain", "percent_mito_binned"], func=["mean", "count_nonzero"]
    ... )
    AnnData object with n_obs × n_vars = 40 × 13714
        obs: 'louvain', 'percent_mito_binned', 'n_obs_aggregated'
        var: 'n_cells'
        layers: 'mean', 'count_nonzero'

    Note that this filters out any combination of groups that wasn't present in the original data.

    """
    if not isinstance(adata, AnnData):
        msg = (
            "sc.get.aggregate is currently only implemented for AnnData input, "
            f"was passed {type(adata)}."
        )
        raise NotImplementedError(msg)
    if axis is None:
        axis = 1 if varm else 0
    axis, axis_name = _resolve_axis(axis)
    mask = _check_mask(adata, mask, axis_name)
    data = adata.X
    if sum(p is not None for p in [varm, obsm, layer]) > 1:
        msg = "Please only provide one (or none) of varm, obsm, or layer"
        raise TypeError(msg)

    if varm is not None:
        if axis != 1:
            msg = "varm can only be used when axis is 1"
            raise ValueError(msg)
        data = adata.varm[varm]
    elif obsm is not None:
        if axis != 0:
            msg = "obsm can only be used when axis is 0"
            raise ValueError(msg)
        data = adata.obsm[obsm]
    elif layer is not None:
        data = adata.layers[layer]
        if axis == 1:
            data = data.T
    elif axis == 1:
        # i.e., all of `varm`, `obsm`, `layers` are None so we use `X` which must be transposed
        data = data.T

    dim_df = getattr(adata, axis_name)
    categorical, new_label_df = _combine_categories(dim_df, by)

    # Add number of obs aggregated into each group
    new_label_df["n_obs_aggregated"] = (
        pd.Series(categorical).value_counts().reindex(new_label_df.index)
    )
    # Actual computation
    layers = _aggregate(
        data,
        by=categorical,
        func=func,
        mask=mask,
        dof=dof,
    )

    # Define new var dataframe
    if obsm or varm:
        if isinstance(data, pd.DataFrame):
            # Check if there could be labels
            var = pd.DataFrame(index=data.columns)
        else:
            # Create them otherwise
            var = pd.DataFrame(index=pd.RangeIndex(data.shape[1]).astype(str))
    else:
        var = getattr(adata, "var" if axis == 0 else "obs")

    # It's all coming together
    result = AnnData(layers=layers, obs=new_label_df, var=var)

    if axis == 1:
        return result.T
    else:
        return result


@singledispatch
def _aggregate(
    data,
    by: pd.Categorical,
    func: AggType | Iterable[AggType],
    *,
    mask: NDArray[np.bool] | None = None,
    dof: int = 1,
) -> dict[AggType, np.ndarray | DaskArray]:
    msg = f"Data type {type(data)} not supported for aggregation"
    raise NotImplementedError(msg)


class MeanVarDict(TypedDict):
    mean: DaskArray
    var: DaskArray


def aggregate_dask_mean_var(
    data: DaskArray,
    by: pd.Categorical,
    *,
    mask: NDArray[np.bool] | None = None,
    dof: int = 1,
) -> MeanVarDict:
    mean = aggregate_dask(data, by, "mean", mask=mask, dof=dof)["mean"]
    sq_mean = aggregate_dask(fau_power(data, 2), by, "mean", mask=mask, dof=dof)["mean"]
    # TODO: If we don't compute here, the results are not deterministic under the process cluster for sparse.
    if isinstance(data._meta, CSRBase):
        sq_mean = sq_mean.compute()
    var = sq_mean - fau_power(mean, 2)
    if dof != 0:
        group_counts = np.bincount(by.codes)
        var *= (group_counts / (group_counts - dof))[:, np.newaxis]
    return MeanVarDict(mean=mean, var=var)


@_aggregate.register(DaskArray)
def aggregate_dask(
    data: DaskArray,
    by: pd.Categorical,
    func: AggType | Iterable[AggType],
    *,
    mask: NDArray[np.bool] | None = None,
    dof: int = 1,
) -> dict[AggType, DaskArray]:
    import dask

    if not isinstance(data._meta, CSBase | np.ndarray):
        msg = f"Got {type(data._meta)} meta in DaskArray but only csr_matrix/csr_array and ndarray are supported."
        raise ValueError(msg)
    chunked_axis, _ = (0, 1) if isinstance(data._meta, CSRBase | np.ndarray) else (1, 0)

    funcs = set([func] if isinstance(func, str) else func)
    if "median" in funcs:
        msg = "Dask median calculation not supported.  If you want a median-of-medians calculation, please open an issue."
        raise NotImplementedError(msg)
    has_mean, has_var = (v in funcs for v in ["mean", "var"])
    funcs_no_var_or_mean = funcs - {"var", "mean"}

    if funcs_no_var_or_mean:
        funcs_list = list(funcs_no_var_or_mean)
        by_codes = np.asarray(by.codes)
        mask_arr = np.asarray(mask) if mask is not None else None

        @dask.delayed
        def aggregate_chunk(block, block_idx):
            subset = slice(block_idx[0], block_idx[1])
            by_subsetted = (
                pd.Categorical.from_codes(by_codes[subset], categories=by.categories)
                if chunked_axis == 0
                else by
            )
            mask_subsetted = (
                mask_arr[subset]
                if (mask_arr is not None and chunked_axis == 0)
                else mask_arr
            )
            return {
                f: _aggregate(block, by_subsetted, f, mask=mask_subsetted, dof=dof)[f]
                for f in funcs_list
            }

        @dask.delayed
        def add_aggs(a, b):
            return {f: a[f] + b[f] for f in funcs_list}

        blocks = data.to_delayed().ravel()

        offset = 0
        delayed_chunks = []
        for i, block in enumerate(blocks):
            block_idx = (offset, offset + data.chunks[chunked_axis][i])
            delayed_chunks.append(aggregate_chunk(block, block_idx))
            offset += data.chunks[chunked_axis][i]

        while len(delayed_chunks) > 1:
            delayed_chunks = [
                add_aggs(delayed_chunks[i], delayed_chunks[i + 1])
                if i + 1 < len(delayed_chunks)
                else delayed_chunks[i]
                for i in range(0, len(delayed_chunks), 2)
            ]

        aggregated = {
            f: dask.array.from_delayed(
                dask.delayed(lambda r, f=f: r[f])(delayed_chunks[0]),
                shape=(len(by.categories), data.shape[1]),
                dtype=np.float64,
            )
            for f in funcs_list
        }

    if has_var:
        aggredated_mean_var = aggregate_dask_mean_var(data, by, mask=mask, dof=dof)
        aggregated["var"] = aggredated_mean_var["var"]
        if has_mean:
            aggregated["mean"] = aggredated_mean_var["mean"]
    # division must come after, not before, the summation for numerical precision
    # i.e., we can't just call map blocks over the mean function.
    elif has_mean:
        group_counts = np.bincount(by.codes)
        aggregated["mean"] = (
            aggregate_dask(data, by, "sum", mask=mask, dof=dof)["sum"]
            / group_counts[:, None]
        )
    return aggregated


@_aggregate.register(pd.DataFrame)
def aggregate_df(data, by, func, *, mask=None, dof=1) -> dict[AggType, np.ndarray]:
    return _aggregate(data.values, by, func, mask=mask, dof=dof)


@_aggregate.register(np.ndarray)
@_aggregate.register(CSBase)
def aggregate_array(
    data: Array,
    by: pd.Categorical,
    func: AggType | Iterable[AggType],
    *,
    mask: NDArray[np.bool] | None = None,
    dof: int = 1,
) -> dict[AggType, np.ndarray]:
    groupby = Aggregate(groupby=by, data=data, mask=mask)
    result = {}

    funcs = set([func] if isinstance(func, str) else func)
    if unknown := funcs - get_literal_vals(AggType):
        msg = f"func {unknown} is not one of {get_literal_vals(AggType)}"
        raise ValueError(msg)

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
    if "median" in funcs:
        agg = groupby.median()
        result["median"] = agg
    return result


def _combine_categories(
    label_df: pd.DataFrame, cols: Collection[str] | str
) -> tuple[pd.Categorical, pd.DataFrame]:
    """Return both the result categories and a dataframe labelling each row."""
    from itertools import product

    if isinstance(cols, str):
        cols = [cols]

    df = pd.DataFrame(
        {c: pd.Categorical(label_df[c]).remove_unused_categories() for c in cols},
    )
    n_categories = [len(df[c].cat.categories) for c in cols]

    # It's like np.concatenate([x for x in product(*[range(n) for n in n_categories])])
    code_combinations = np.indices(n_categories).reshape(len(n_categories), -1)
    result_categories = pd.Index([
        "_".join(map(str, x)) for x in product(*[df[c].cat.categories for c in cols])
    ])

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
    np.cumprod(n_categories[::-1], out=factors[1:])
    factors = factors[:-1][::-1]

    code_array = np.zeros((len(cols), df.shape[0]), dtype=np.int32)
    for i, c in enumerate(cols):
        code_array[i] = df[c].cat.codes
    code_array *= factors[:, None]

    codes = code_array.sum(axis=0)
    codes = np.where(np.any(code_array < 0, axis=0), -1, codes)
    result_categorical = pd.Categorical.from_codes(codes, categories=result_categories)

    # Filter unused categories
    result_categorical = result_categorical.remove_unused_categories()
    new_label_df = new_label_df.loc[result_categorical.categories]

    return result_categorical, new_label_df


def sparse_indicator(
    categorical: pd.Categorical,
    *,
    mask: NDArray[np.bool] | None = None,
) -> sparse.coo_array:
    # TODO: why is this float64.  This is a scanpy 2.0 problem maybe?
    mask = (
        np.broadcast_to(1.0, len(categorical)) if mask is None else mask.astype("uint8")
    )
    # can’t have -1s in the codes, but (as long as it’s valid), the value is ignored, so set to 0 where masked
    codes = np.where(mask, categorical.codes, 0)
    a = sparse.coo_array(
        (mask, (codes, np.arange(len(categorical)))),
        shape=(len(categorical.categories), len(categorical)),
    )
    return a

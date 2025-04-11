from __future__ import annotations

from functools import singledispatch, wraps
from typing import TYPE_CHECKING
from warnings import warn

import numba
import numpy as np
import pandas as pd
from scipy import sparse

from scanpy.preprocessing._distributed import materialize_as_ndarray
from scanpy.preprocessing._utils import _get_mean_var

from .._compat import CSBase, CSRBase, DaskArray, _register_union, njit
from .._utils import _doc_params, axis_nnz, axis_sum
from ._docs import (
    doc_adata_basic,
    doc_expr_reps,
    doc_obs_qc_args,
    doc_obs_qc_returns,
    doc_qc_metric_naming,
    doc_var_qc_returns,
)

if TYPE_CHECKING:
    from collections.abc import Collection

    from anndata import AnnData


def _choose_mtx_rep(adata, *, use_raw: bool = False, layer: str | None = None):
    is_layer = layer is not None
    if use_raw and is_layer:
        msg = (
            "Cannot use expression from both layer and raw. You provided:"
            f"{use_raw=!r} and {layer=!r}"
        )
        raise ValueError(msg)
    if is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    else:
        return adata.X


@_doc_params(
    doc_adata_basic=doc_adata_basic,
    doc_expr_reps=doc_expr_reps,
    doc_obs_qc_args=doc_obs_qc_args,
    doc_qc_metric_naming=doc_qc_metric_naming,
    doc_obs_qc_returns=doc_obs_qc_returns,
)
def describe_obs(  # noqa: PLR0913
    adata: AnnData,
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    qc_vars: Collection[str] = (),
    percent_top: Collection[int] | None = (50, 100, 200, 500),
    layer: str | None = None,
    use_raw: bool = False,
    log1p: bool | None = True,
    inplace: bool = False,
    X=None,
    parallel=None,
) -> pd.DataFrame | None:
    """Describe observations of anndata.

    Calculates a number of qc metrics for observations in AnnData object. See
    section `Returns` for a description of those metrics.

    Note that this method can take a while to compile on the first call. That
    result is then cached to disk to be used later.

    Params
    ------
    {doc_adata_basic}
    {doc_qc_metric_naming}
    {doc_obs_qc_args}
    {doc_expr_reps}
    log1p
        Add `log1p` transformed metrics.
    inplace
        Whether to place calculated metrics in `adata.obs`.
    X
        Matrix to calculate values on. Meant for internal usage.

    Returns
    -------
    QC metrics for observations in adata. If inplace, values are placed into
    the AnnData's `.obs` dataframe.

    {doc_obs_qc_returns}

    """
    if parallel is not None:
        warn(
            "Argument `parallel` is deprecated, and currently has no effect.",
            FutureWarning,
            stacklevel=2,
        )
    # Handle whether X is passed
    if X is None:
        X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
        if isinstance(X, sparse.coo_matrix):
            X = sparse.csr_matrix(X)  # COO not subscriptable  # noqa: TID251
        if isinstance(X, CSBase):
            X.eliminate_zeros()
    obs_metrics = pd.DataFrame(index=adata.obs_names)
    obs_metrics[f"n_{var_type}_by_{expr_type}"] = materialize_as_ndarray(
        axis_nnz(X, axis=1)
    )
    if log1p:
        obs_metrics[f"log1p_n_{var_type}_by_{expr_type}"] = np.log1p(
            obs_metrics[f"n_{var_type}_by_{expr_type}"]
        )
    obs_metrics[f"total_{expr_type}"] = np.ravel(axis_sum(X, axis=1))
    if log1p:
        obs_metrics[f"log1p_total_{expr_type}"] = np.log1p(
            obs_metrics[f"total_{expr_type}"]
        )
    if percent_top:
        percent_top = sorted(percent_top)
        proportions = top_segment_proportions(X, percent_top)
        for i, n in enumerate(percent_top):
            obs_metrics[f"pct_{expr_type}_in_top_{n}_{var_type}"] = (
                proportions[:, i] * 100
            )
    for qc_var in qc_vars:
        obs_metrics[f"total_{expr_type}_{qc_var}"] = np.ravel(
            axis_sum(X[:, adata.var[qc_var].values], axis=1)
        )
        if log1p:
            obs_metrics[f"log1p_total_{expr_type}_{qc_var}"] = np.log1p(
                obs_metrics[f"total_{expr_type}_{qc_var}"]
            )
        obs_metrics[f"pct_{expr_type}_{qc_var}"] = (
            obs_metrics[f"total_{expr_type}_{qc_var}"]
            / obs_metrics[f"total_{expr_type}"]
            * 100
        )
    if inplace:
        adata.obs[obs_metrics.columns] = obs_metrics
    else:
        return obs_metrics
    return None


@_doc_params(
    doc_adata_basic=doc_adata_basic,
    doc_expr_reps=doc_expr_reps,
    doc_qc_metric_naming=doc_qc_metric_naming,
    doc_var_qc_returns=doc_var_qc_returns,
)
def describe_var(
    adata: AnnData,
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    layer: str | None = None,
    use_raw: bool = False,
    inplace: bool = False,
    log1p: bool = True,
    X: CSBase | sparse.coo_matrix | np.ndarray | None = None,
) -> pd.DataFrame | None:
    """Describe variables of anndata.

    Calculates a number of qc metrics for variables in AnnData object. See
    section `Returns` for a description of those metrics.

    Params
    ------
    {doc_adata_basic}
    {doc_qc_metric_naming}
    {doc_expr_reps}
    inplace
        Whether to place calculated metrics in `adata.var`.
    X
        Matrix to calculate values on. Meant for internal usage.

    Returns
    -------
    QC metrics for variables in adata. If inplace, values are placed into the
    AnnData's `.var` dataframe.

    {doc_var_qc_returns}

    """
    # Handle whether X is passed
    if X is None:
        X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
        if isinstance(X, sparse.coo_matrix):
            X = sparse.csr_matrix(X)  # COO not subscriptable  # noqa: TID251
        if isinstance(X, CSBase):
            X.eliminate_zeros()
    var_metrics = pd.DataFrame(index=adata.var_names)
    var_metrics[f"n_cells_by_{expr_type}"], var_metrics[f"mean_{expr_type}"] = (
        materialize_as_ndarray((axis_nnz(X, axis=0), _get_mean_var(X, axis=0)[0]))
    )
    if log1p:
        var_metrics[f"log1p_mean_{expr_type}"] = np.log1p(
            var_metrics[f"mean_{expr_type}"]
        )
    var_metrics[f"pct_dropout_by_{expr_type}"] = (
        1 - var_metrics[f"n_cells_by_{expr_type}"] / X.shape[0]
    ) * 100
    var_metrics[f"total_{expr_type}"] = np.ravel(axis_sum(X, axis=0))
    if log1p:
        var_metrics[f"log1p_total_{expr_type}"] = np.log1p(
            var_metrics[f"total_{expr_type}"]
        )
    if inplace:
        adata.var[var_metrics.columns] = var_metrics
        return None
    return var_metrics


@_doc_params(
    doc_adata_basic=doc_adata_basic,
    doc_expr_reps=doc_expr_reps,
    doc_obs_qc_args=doc_obs_qc_args,
    doc_qc_metric_naming=doc_qc_metric_naming,
    doc_obs_qc_returns=doc_obs_qc_returns,
    doc_var_qc_returns=doc_var_qc_returns,
)
def calculate_qc_metrics(
    adata: AnnData,
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    qc_vars: Collection[str] | str = (),
    percent_top: Collection[int] | None = (50, 100, 200, 500),
    layer: str | None = None,
    use_raw: bool = False,
    inplace: bool = False,
    log1p: bool = True,
    parallel: bool | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Calculate quality control metrics.

    Calculates a number of qc metrics for an AnnData object, see section
    `Returns` for specifics. Largely based on `calculateQCMetrics` from scater
    :cite:p:`McCarthy2017`. Currently is most efficient on a sparse CSR or dense matrix.

    Note that this method can take a while to compile on the first call. That
    result is then cached to disk to be used later.

    Parameters
    ----------
    {doc_adata_basic}
    {doc_qc_metric_naming}
    {doc_obs_qc_args}
    {doc_expr_reps}
    inplace
        Whether to place calculated metrics in `adata`'s `.obs` and `.var`.
    log1p
        Set to `False` to skip computing `log1p` transformed annotations.

    Returns
    -------
    Depending on `inplace` returns calculated metrics
    (as :class:`~pandas.DataFrame`) or updates `adata`'s `obs` and `var`.

    {doc_obs_qc_returns}

    {doc_var_qc_returns}

    Example
    -------
    Calculate qc metrics for visualization.

    .. plot::
        :context: close-figs

        import scanpy as sc
        import seaborn as sns

        pbmc = sc.datasets.pbmc3k()
        pbmc.var["mito"] = pbmc.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(pbmc, qc_vars=["mito"], inplace=True)
        sns.jointplot(
            data=pbmc.obs,
            x="log1p_total_counts",
            y="log1p_n_genes_by_counts",
            kind="hex",
        )

    .. plot::
        :context: close-figs

        sns.histplot(pbmc.obs["pct_counts_mito"])

    """
    if parallel is not None:
        warn(
            "Argument `parallel` is deprecated, and currently has no effect.",
            FutureWarning,
            stacklevel=2,
        )
    # Pass X so I only have to do it once
    X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
    if isinstance(X, sparse.coo_matrix):
        X = sparse.csr_matrix(X)  # COO not subscriptable  # noqa: TID251
    if isinstance(X, CSBase):
        X.eliminate_zeros()

    # Convert qc_vars to list if str
    if isinstance(qc_vars, str):
        qc_vars = [qc_vars]

    obs_metrics = describe_obs(
        adata,
        expr_type=expr_type,
        var_type=var_type,
        qc_vars=qc_vars,
        percent_top=percent_top,
        inplace=inplace,
        X=X,
        log1p=log1p,
    )
    var_metrics = describe_var(
        adata,
        expr_type=expr_type,
        var_type=var_type,
        inplace=inplace,
        X=X,
        log1p=log1p,
    )

    if not inplace:
        return obs_metrics, var_metrics


def top_proportions(mtx: np.ndarray | CSBase | sparse.coo_matrix, n: int):
    """Calculate cumulative proportions of top expressed genes.

    Parameters
    ----------
    mtx
        Matrix, where each row is a sample, each column a feature.
    n
        Rank to calculate proportions up to. Value is treated as 1-indexed,
        `n=50` will calculate cumulative proportions up to the 50th most
        expressed gene.

    """
    if isinstance(mtx, CSBase | sparse.coo_matrix):
        if not isinstance(mtx, CSRBase):
            mtx = sparse.csr_matrix(mtx)  # noqa: TID251
        # Allowing numba to do more
        return top_proportions_sparse_csr(mtx.data, mtx.indptr, np.array(n))
    else:
        return top_proportions_dense(mtx, n)


def top_proportions_dense(mtx, n):
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.argpartition, 1, -mtx, n - 1)
    partitioned = partitioned[:, :n]
    values = np.zeros_like(partitioned, dtype=np.float64)
    for i in range(partitioned.shape[0]):
        vec = mtx[i, partitioned[i, :]]  # Not a view
        vec[::-1].sort()  # Sorting on a reversed view (e.g. a descending sort)
        vec = np.cumsum(vec) / sums[i]
        values[i, :] = vec
    return values


def top_proportions_sparse_csr(data, indptr, n):
    values = np.zeros((indptr.size - 1, n), dtype=np.float64)
    for i in numba.prange(indptr.size - 1):
        start, end = indptr[i], indptr[i + 1]
        vec = np.zeros(n, dtype=np.float64)
        if end - start <= n:
            vec[: end - start] = data[start:end]
            total = vec.sum()
        else:
            vec[:] = -(np.partition(-data[start:end], n - 1)[:n])
            total = (data[start:end]).sum()  # Is this not just vec.sum()?
        vec[::-1].sort()
        values[i, :] = vec.cumsum() / total
    return values


def check_ns(func):
    @wraps(func)
    def check_ns_inner(
        mtx: np.ndarray | CSBase | sparse.coo_matrix | DaskArray, ns: Collection[int]
    ):
        if not (max(ns) <= mtx.shape[1] and min(ns) > 0):
            msg = "Positions outside range of features."
            raise IndexError(msg)
        return func(mtx, ns)

    return check_ns_inner


@singledispatch
@check_ns
def top_segment_proportions(mtx: np.ndarray, ns: Collection[int]) -> np.ndarray:
    """Calculate total percentage of counts in top ns genes.

    Parameters
    ----------
    mtx
        Matrix, where each row is a sample, each column a feature.
    ns
        Positions to calculate cumulative proportion at. Values are considered
        1-indexed, e.g. `ns=[50]` will calculate cumulative proportion up to
        the 50th most expressed gene.

    """
    # Currently ns is considered to be 1 indexed
    ns = np.sort(ns)
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.partition, 1, mtx, mtx.shape[1] - ns)[:, ::-1][
        :, : ns[-1]
    ]
    values = np.zeros((mtx.shape[0], len(ns)))
    acc = np.zeros(mtx.shape[0])
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]


@top_segment_proportions.register(DaskArray)
@check_ns
def _(mtx: DaskArray, ns: Collection[int]) -> DaskArray:
    if not isinstance(mtx._meta, CSRBase | np.ndarray):
        msg = f"DaskArray must have csr matrix or ndarray meta, got {mtx._meta}."
        raise ValueError(msg)
    return mtx.map_blocks(
        lambda x: top_segment_proportions(x, ns), meta=np.array([])
    ).compute()


@_register_union(top_segment_proportions, CSBase)
@top_segment_proportions.register(sparse.coo_matrix)
@check_ns
def _(mtx: CSBase | sparse.coo_matrix, ns: Collection[int]) -> DaskArray:
    if not isinstance(mtx, CSRBase):
        mtx = sparse.csr_matrix(mtx)  # noqa: TID251
    return top_segment_proportions_sparse_csr(mtx.data, mtx.indptr, np.array(ns))


@njit
def top_segment_proportions_sparse_csr(data, indptr, ns):
    # work around https://github.com/numba/numba/issues/5056
    indptr = indptr.astype(np.int64)
    ns = ns.astype(np.int64)
    ns = np.sort(ns)
    maxidx = ns[-1]
    sums = np.zeros((indptr.size - 1), dtype=data.dtype)
    values = np.zeros((indptr.size - 1, len(ns)), dtype=np.float64)
    # Just to keep it simple, as a dense matrix
    partitioned = np.zeros((indptr.size - 1, maxidx), dtype=data.dtype)
    for i in numba.prange(indptr.size - 1):
        start, end = indptr[i], indptr[i + 1]
        sums[i] = np.sum(data[start:end])
        if end - start <= maxidx:
            partitioned[i, : end - start] = data[start:end]
        elif (end - start) > maxidx:
            partitioned[i, :] = -(np.partition(-data[start:end], maxidx))[:maxidx]
        partitioned[i, :] = np.partition(partitioned[i, :], maxidx - ns)
    partitioned = partitioned[:, ::-1][:, : ns[-1]]
    acc = np.zeros((indptr.size - 1), dtype=data.dtype)
    prev = 0
    # can’t use enumerate due to https://github.com/numba/numba/issues/2625
    for j in range(ns.size):
        acc += partitioned[:, prev : ns[j]].sum(axis=1)
        values[:, j] = acc
        prev = ns[j]
    return values / sums.reshape((indptr.size - 1, 1))

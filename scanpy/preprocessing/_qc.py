<<<<<<< HEAD
import numba
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse, isspmatrix_csr, isspmatrix_coo
from sklearn.utils.sparsefuncs import mean_variance_axis


def calculate_qc_metrics(adata, expr_type="counts", var_type="genes", qc_vars=(),
                         percent_top=(50, 100, 200, 500), inplace=False):
    """
    Calculate quality control metrics.

    Calculates a number of qc metrics for an AnnData object, see section
    `Returns` for specifics. Largely based on `calculateQCMetrics` from scater
    [McCarthy17]_. Currently is most efficient on a sparse CSR or dense matrix.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    expr_type : `str`, optional (default: `"counts"`)
        Name of kind of values in X.
    var_type : `str`, optional (default: `"genes"`)
        The kind of thing the variables are.
    qc_vars : `Container`, optional (default: `()`)
        Keys for boolean columns of `.var` which identify variables you could 
        want to control for (e.g. "ERCC" or "mito").
    percent_top : `Container[int]`, optional (default: `(50, 100, 200, 500)`)
        Which proportions of top genes to cover. If empty or `None` don't
        calculate. Values are considered 1-indexed, `percent_top=[50]` finds
        cumulative proportion to the 50th most expressed gene.
    inplace : bool, optional (default: `False`)
        Whether to place calculated metrics in `.obs` and `.var`

    Returns
    -------
    Union[NoneType, Tuple[pd.DataFrame, pd.DataFrame]]
        Depending on `inplace` returns calculated metrics (`pd.DataFrame`) or
        updates `adata`'s `obs` and `var`.

        Observation level metrics include:

        * `total_{var_type}_by_{expr_type}`
            E.g. "total_genes_by_counts". Number of genes with positive counts
            in a cell.
        * `total_{expr_type}`
            E.g. "total_counts". Total number of counts for a cell.
        * `pct_{expr_type}_in_top_{n}_{var_type}` - for `n` in `percent_top`
            E.g. "pct_counts_in_top_50_genes". Cumulative percentage of counts
            for 50 most expressed genes in a cell.
        * `total_{expr_type}_{qc_var}` - for `qc_var` in `qc_vars`
            E.g. "total_counts_mito". Total number of counts for variabes in
            `qc_vars`.
        * `pct_{expr_type}_{qc_var}` - for `qc_var` in `qc_vars`
            E.g. "pct_counts_mito". Proportion of total counts for a cell which
            are mitochondrial.

        Variable level metrics include:

        * `total_{expr_type}`
            E.g. "total_counts". Sum of counts for a gene.
        * `mean_{expr_type}`
            E.g. "mean counts". Mean expression over all cells.
        * `n_cells_by_{expr_type}`
            E.g. "n_cells_by_counts". Number of cells this expression is
            measured in.
        * `pct_dropout_by_{expr_type}`
            E.g. "pct_dropout_by_counts". Percentage of cells this feature does
            not appear in.


    Example
    -------
    Calculate qc metrics for visualization.

    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sns.jointplot(adata.obs, "log1p_total_counts", "log1p_n_genes_by_counts", kind="hex")
    """
    X = adata.X
    obs_metrics = pd.DataFrame(index=adata.obs_names)
    var_metrics = pd.DataFrame(index=adata.var_names)
    if isspmatrix_coo(X):
        X = csr_matrix(X)  # COO not subscriptable
    if issparse(X):
        X.eliminate_zeros()
    # Calculate obs metrics
    if issparse(X):
        obs_metrics["n_{var_type}_by_{expr_type}"] = X.getnnz(axis=1)
    else:
        obs_metrics["n_{var_type}_by_{expr_type}"] = np.count_nonzero(X, axis=1)
    obs_metrics["log1p_n_{var_type}_by_{expr_type}"] = np.log1p(
        obs_metrics["n_{var_type}_by_{expr_type}"])
    obs_metrics["total_{expr_type}"] = X.sum(axis=1)
    obs_metrics["log1p_total_{expr_type}"] = np.log1p(
        obs_metrics["total_{expr_type}"])
    if percent_top:
        percent_top = sorted(percent_top)
        proportions = top_segment_proportions(X, percent_top)
        # Since there are local loop variables, formatting must occur in their scope
        # Probably worth looking into a python3.5 compatable way to make this better
        for i, n in enumerate(percent_top):
            obs_metrics["pct_{expr_type}_in_top_{n}_{var_type}".format(**locals())] = \
                proportions[:, i] * 100
    for qc_var in qc_vars:
        obs_metrics["total_{expr_type}_{qc_var}".format(**locals())] = \
            X[:, adata.var[qc_var].values].sum(axis=1)
        obs_metrics["log1p_total_{expr_type}_{qc_var}".format(**locals())] = \
            np.log1p(
                obs_metrics["total_{expr_type}_{qc_var}".format(**locals())])
        # "total_{expr_type}" not formatted yet
        obs_metrics["pct_{expr_type}_{qc_var}".format(**locals())] = \
            obs_metrics["total_{expr_type}_{qc_var}".format(**locals())] / \
            obs_metrics["total_{expr_type}"] * 100
    # Calculate var metrics
=======
from typing import Optional, Tuple, Collection, Union
from warnings import warn

import numba
import numpy as np
import pandas as pd
from scipy.sparse import issparse, isspmatrix_csr, isspmatrix_coo
from scipy.sparse import spmatrix, csr_matrix
from sklearn.utils.sparsefuncs import mean_variance_axis

from anndata import AnnData
from ._docs import (
    doc_expr_reps,
    doc_obs_qc_args,
    doc_qc_metric_naming,
    doc_obs_qc_returns,
    doc_var_qc_returns,
    doc_adata_basic,
)
from .._utils import _doc_params


def _choose_mtx_rep(adata, use_raw=False, layer=None):
    is_layer = layer is not None
    if use_raw and is_layer:
        raise ValueError(
            "Cannot use expression from both layer and raw. You provided:"
            f"'use_raw={use_raw}' and 'layer={layer}'"
        )
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
def describe_obs(
    adata: AnnData,
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    qc_vars: Collection[str] = (),
    percent_top: Collection[int] = (50, 100, 200, 500),
    layer: Optional[str] = None,
    use_raw: bool = False,
    inplace: bool = False,
    X=None,
    parallel=None,
) -> Optional[pd.DataFrame]:
    """\
    Describe observations of anndata.

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
            FutureWarning
        )
    # Handle whether X is passed
    if X is None:
        X = _choose_mtx_rep(adata, use_raw, layer)
        if isspmatrix_coo(X):
            X = csr_matrix(X)  # COO not subscriptable
        if issparse(X):
            X.eliminate_zeros()
    obs_metrics = pd.DataFrame(index=adata.obs_names)
    if issparse(X):
        obs_metrics[f"n_{var_type}_by_{expr_type}"] = X.getnnz(axis=1)
    else:
        obs_metrics[f"n_{var_type}_by_{expr_type}"] = np.count_nonzero(X, axis=1)
    obs_metrics[f"log1p_n_{var_type}_by_{expr_type}"] = np.log1p(
        obs_metrics[f"n_{var_type}_by_{expr_type}"]
    )
    obs_metrics[f"total_{expr_type}"] = X.sum(axis=1)
    obs_metrics[f"log1p_total_{expr_type}"] = np.log1p(obs_metrics[f"total_{expr_type}"])
    if percent_top:
        percent_top = sorted(percent_top)
        proportions = top_segment_proportions(X, percent_top)
        for i, n in enumerate(percent_top):
            obs_metrics[f"pct_{expr_type}_in_top_{n}_{var_type}"] = (
                proportions[:, i] * 100
            )
    for qc_var in qc_vars:
        obs_metrics[f"total_{expr_type}_{qc_var}"] = (
            X[:, adata.var[qc_var].values].sum(axis=1)
        )
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
    layer: Optional[str] = None,
    use_raw: bool = False,
    inplace=False,
    X=None,
) -> Optional[pd.DataFrame]:
    """\
    Describe variables of anndata.

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
        X = _choose_mtx_rep(adata, use_raw, layer)
        if isspmatrix_coo(X):
            X = csr_matrix(X)  # COO not subscriptable
        if issparse(X):
            X.eliminate_zeros()
    var_metrics = pd.DataFrame(index=adata.var_names)
>>>>>>> upstream/master
    if issparse(X):
        # Current memory bottleneck for csr matrices:
        var_metrics["n_cells_by_{expr_type}"] = X.getnnz(axis=0)
        var_metrics["mean_{expr_type}"] = mean_variance_axis(X, axis=0)[0]
    else:
        var_metrics["n_cells_by_{expr_type}"] = np.count_nonzero(X, axis=0)
        var_metrics["mean_{expr_type}"] = X.mean(axis=0)
<<<<<<< HEAD
    var_metrics["log1p_mean_{expr_type}"] = np.log1p(
        var_metrics["mean_{expr_type}"])
    var_metrics["pct_dropout_by_{expr_type}"] = \
        (1 - var_metrics["n_cells_by_{expr_type}"] / X.shape[0]) * 100
    var_metrics["total_{expr_type}"] = np.ravel(X.sum(axis=0))
    var_metrics["log1p_total_{expr_type}"] = np.log1p(
        var_metrics["total_{expr_type}"])
    # Format strings
    for df in obs_metrics, var_metrics:
        new_colnames = []
        for col in df.columns:
            new_colnames.append(col.format(**locals()))
        df.columns = new_colnames
    # Return
    if inplace:
        adata.obs[obs_metrics.columns] = obs_metrics
        adata.var[var_metrics.columns] = var_metrics
    else:
        return obs_metrics, var_metrics

def top_proportions(mtx, n):
    """
=======
    var_metrics["log1p_mean_{expr_type}"] = np.log1p(var_metrics["mean_{expr_type}"])
    var_metrics["pct_dropout_by_{expr_type}"] = (
        1 - var_metrics["n_cells_by_{expr_type}"] / X.shape[0]
    ) * 100
    var_metrics["total_{expr_type}"] = np.ravel(X.sum(axis=0))
    var_metrics["log1p_total_{expr_type}"] = np.log1p(var_metrics["total_{expr_type}"])
    # Relabel
    new_colnames = []
    for col in var_metrics.columns:
        new_colnames.append(col.format(**locals()))
    var_metrics.columns = new_colnames
    if inplace:
        adata.var[var_metrics.columns] = var_metrics
    else:
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
    qc_vars: Collection[str] = (),
    percent_top: Collection[int] = (50, 100, 200, 500),
    layer: Optional[str] = None,
    use_raw: bool = False,
    inplace: bool = False,
    parallel: Optional[bool] = None,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """\
    Calculate quality control metrics.

    Calculates a number of qc metrics for an AnnData object, see section
    `Returns` for specifics. Largely based on `calculateQCMetrics` from scater
    [McCarthy17]_. Currently is most efficient on a sparse CSR or dense matrix.

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

    Returns
    -------
    Depending on `inplace` returns calculated metrics
    (as :class:`~pandas.DataFrame`) or updates `adata`'s `obs` and `var`.

    {doc_obs_qc_returns}

    {doc_var_qc_returns}

    Example
    -------
    Calculate qc metrics for visualization.

    >>> import scanpy as sc
    >>> import seaborn as sns
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sns.jointplot(
            "log1p_total_counts", "log1p_n_genes_by_counts",
            data=adata.obs, kind="hex"
        )
    """
    if parallel is not None:
        warn(
            "Argument `parallel` is deprecated, and currently has no effect.",
            FutureWarning
        )
    # Pass X so I only have to do it once
    X = _choose_mtx_rep(adata, use_raw, layer)
    if isspmatrix_coo(X):
        X = csr_matrix(X)  # COO not subscriptable
    if issparse(X):
        X.eliminate_zeros()

    obs_metrics = describe_obs(
        adata,
        expr_type=expr_type,
        var_type=var_type,
        qc_vars=qc_vars,
        percent_top=percent_top,
        inplace=inplace,
        X=X,
    )
    var_metrics = describe_var(
        adata, expr_type=expr_type, var_type=var_type, inplace=inplace, X=X
    )

    if not inplace:
        return obs_metrics, var_metrics


def top_proportions(mtx: Union[np.array, spmatrix], n: int):
    """\
>>>>>>> upstream/master
    Calculates cumulative proportions of top expressed genes

    Parameters
    ----------
<<<<<<< HEAD
    mtx : `Union[np.array, sparse.spmatrix]`
        Matrix, where each row is a sample, each column a feature.
    n : `int`
=======
    mtx
        Matrix, where each row is a sample, each column a feature.
    n
>>>>>>> upstream/master
        Rank to calculate proportions up to. Value is treated as 1-indexed,
        `n=50` will calculate cumulative proportions up to the 50th most
        expressed gene.
    """
    if issparse(mtx):
        if not isspmatrix_csr(mtx):
            mtx = csr_matrix(mtx)
        # Allowing numba to do more
<<<<<<< HEAD
        return top_proportions_sparse_csr(mtx.data, mtx.indptr, n)
    else:
        return top_proportions_dense(mtx, n)

@numba.jit
def top_proportions_dense(mtx, n):
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.argpartition, 1, -mtx, n-1)
=======
        return top_proportions_sparse_csr(mtx.data, mtx.indptr, np.array(n))
    else:
        return top_proportions_dense(mtx, n)


def top_proportions_dense(mtx, n):
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.argpartition, 1, -mtx, n - 1)
>>>>>>> upstream/master
    partitioned = partitioned[:, :n]
    values = np.zeros_like(partitioned, dtype=np.float64)
    for i in range(partitioned.shape[0]):
        vec = mtx[i, partitioned[i, :]]  # Not a view
        vec[::-1].sort()  # Sorting on a reversed view (e.g. a descending sort)
        vec = np.cumsum(vec) / sums[i]
        values[i, :] = vec
    return values

<<<<<<< HEAD
@numba.jit(parallel=True)
def top_proportions_sparse_csr(data, indptr, n):
    values = np.zeros((indptr.size-1, n), dtype=np.float64)
    for i in numba.prange(indptr.size-1):
        start, end = indptr[i], indptr[i+1]
        vec = np.zeros(n, dtype=np.float64)
        if end - start <= n:
            vec[:end-start] = data[start:end]
            total = vec.sum()
        else:
            vec[:] = -(np.partition(-data[start:end], n-1)[:n])
=======

def top_proportions_sparse_csr(data, indptr, n):
    values = np.zeros((indptr.size - 1, n), dtype=np.float64)
    for i in numba.prange(indptr.size - 1):
        start, end = indptr[i], indptr[i + 1]
        vec = np.zeros(n, dtype=np.float64)
        if end - start <= n:
            vec[:end - start] = data[start:end]
            total = vec.sum()
        else:
            vec[:] = -(np.partition(-data[start:end], n - 1)[:n])
>>>>>>> upstream/master
            total = (data[start:end]).sum()  # Is this not just vec.sum()?
        vec[::-1].sort()
        values[i, :] = vec.cumsum() / total
    return values


<<<<<<< HEAD
def top_segment_proportions(mtx, ns):
=======
def top_segment_proportions(
    mtx: Union[np.array, spmatrix], ns: Collection[int]
) -> np.ndarray:
>>>>>>> upstream/master
    """
    Calculates total percentage of counts in top ns genes.

    Parameters
    ----------
<<<<<<< HEAD
    mtx : `Union[np.array, sparse.spmatrix]`
        Matrix, where each row is a sample, each column a feature.
    ns : `Container[Int]`
=======
    mtx
        Matrix, where each row is a sample, each column a feature.
    ns
>>>>>>> upstream/master
        Positions to calculate cumulative proportion at. Values are considered
        1-indexed, e.g. `ns=[50]` will calculate cumulative proportion up to
        the 50th most expressed gene.
    """
    # Pretty much just does dispatch
    if not (max(ns) <= mtx.shape[1] and min(ns) > 0):
        raise IndexError("Positions outside range of features.")
    if issparse(mtx):
        if not isspmatrix_csr(mtx):
            mtx = csr_matrix(mtx)
<<<<<<< HEAD
        return top_segment_proportions_sparse_csr(mtx.data, mtx.indptr,
                                                  np.array(ns, dtype=np.int))
    else:
        return top_segment_proportions_dense(mtx, ns)

def top_segment_proportions_dense(mtx, ns):
    # Currently ns is considered to be 1 indexed
    ns = np.sort(ns)
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(
        np.partition, 1, mtx, mtx.shape[1] - ns)[:, ::-1][:, :ns[-1]]
=======
        return top_segment_proportions_sparse_csr(
            mtx.data, mtx.indptr, np.array(ns, dtype=np.int)
        )
    else:
        return top_segment_proportions_dense(mtx, ns)


def top_segment_proportions_dense(
    mtx: Union[np.array, spmatrix], ns: Collection[int]
) -> np.ndarray:
    # Currently ns is considered to be 1 indexed
    ns = np.sort(ns)
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.partition, 1, mtx, mtx.shape[1] - ns)[:, ::-1][:, :ns[-1]]
>>>>>>> upstream/master
    values = np.zeros((mtx.shape[0], len(ns)))
    acc = np.zeros((mtx.shape[0]))
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]

<<<<<<< HEAD
@numba.jit(parallel=True)
=======

@numba.njit(cache=True, parallel=True)
>>>>>>> upstream/master
def top_segment_proportions_sparse_csr(data, indptr, ns):
    ns = np.sort(ns)
    maxidx = ns[-1]
    sums = np.zeros((indptr.size - 1), dtype=data.dtype)
<<<<<<< HEAD
    values = np.zeros((indptr.size-1, len(ns)), dtype=np.float64)
    # Just to keep it simple, as a dense matrix
    partitioned = np.zeros((indptr.size-1, maxidx), dtype=data.dtype)
    for i in numba.prange(indptr.size - 1):
        start, end = indptr[i], indptr[i+1]
        sums[i] = np.sum(data[start:end])
        if end - start <= maxidx:
            partitioned[i, :end-start] = data[start:end]
        elif (end - start) > maxidx:
            partitioned[i, :] = - \
                (np.partition(-data[start:end], maxidx))[:maxidx]
    partitioned = np.apply_along_axis(
        np.partition, 1, partitioned, maxidx - ns)[:, ::-1][:, :ns[-1]]
    acc = np.zeros((indptr.size-1), dtype=data.dtype)
=======
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
    partitioned = partitioned[:, ::-1][:, :ns[-1]]
    acc = np.zeros((indptr.size - 1), dtype=data.dtype)
>>>>>>> upstream/master
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
<<<<<<< HEAD
    return values / sums[:, None]
=======
    return values / sums.reshape((indptr.size - 1, 1))
>>>>>>> upstream/master

import numba
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse, isspmatrix_csr, isspmatrix_coo


def calculate_qc_metrics(adata, exprs_values="counts", feature_controls=(),
                         percent_top=(50, 100, 200, 500), inplace=False):
    """
    Calculate quality control metrics.

    Calculates a number of qc metrics for an AnnData object, largely based on
    `calculateQCMetrics` from scater [McCarthy17]_. Currently is most efficient
    on a sparse CSR or dense matrix.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    exprs_values : `str`, optional (default: `"counts"`)
        Name of kind of values in X.
    feature_controls : `Container`, optional (default: `()`)
        Keys for boolean columns of `.var` which identify feature controls,
        e.g. "ERCC" or "mito".
    percent_top : `Container[int]`, optional (default: `(50, 100, 200, 500)`)
        Which proportions of top genes to cover. If empty or `None` don't
        calculate.
    inplace : bool, optional (default: `False`)
        Whether to place calculated metrics in `.obs` and `.var`

    Returns
    -------
    Union[NoneType, Tuple[pd.DataFrame, pd.DataFrame]]
        Depending on `inplace` returns calculated metrics (`pd.DataFrame`) or
        updates `adata`'s `obs` and `var`.

        Observation level metrics include:

        * `total_features_by_{exprs_values}`
        * `total_{expr_values}`
        * `pct_{expr_values}_in_top_{n}_features` - for `n` in `percent_top`
        * `total_{exprs_values}_{feature_control}` - for each `feature_control`
        * `pct_{exprs_values}_{feature_control}` - for each `feature_control`

        Variable level metrics include:

        * `total_{expr_values}`
        * `mean_{expr_values}`
        * `n_cells_by_{expr_values}`
        * `pct_dropout_by_{expr_values}`
    """
    if isspmatrix_coo(adata.X):
        X = csr_matrix(adata.X)  # COO not subscriptable
    else:
        X = adata.X
    obs_metrics = pd.DataFrame(index=adata.obs_names)
    var_metrics = pd.DataFrame(index=adata.var_names)
    # Calculate obs metrics
    obs_metrics["total_features_by_{exprs_values}"] = (
        X != 0).sum(axis=1)
    obs_metrics["log1p_total_features_by_{exprs_values}"] = np.log1p(
        obs_metrics["total_features_by_{exprs_values}"])
    obs_metrics["total_{exprs_values}"] = X.sum(axis=1)
    obs_metrics["log1p_total_{exprs_values}"] = np.log1p(
        obs_metrics["total_{exprs_values}"])
    proportions = top_segment_proportions(X, percent_top)
    # Since there are local loop variables, formatting must occur in their scope
    # Probably worth looking into a python3.5 compatable way to make this better
    for i, n in enumerate(percent_top):
        obs_metrics["pct_{exprs_values}_in_top_{n}_features".format(**locals())] = \
            proportions[:, i] * 100
    for feature_control in feature_controls:
        obs_metrics["total_{exprs_values}_{feature_control}".format(**locals())] = \
            X[:, adata.var[feature_control].values].sum(axis=1)
        obs_metrics["log1p_total_{exprs_values}_{feature_control}".format(**locals())] = \
            np.log1p(
                obs_metrics["total_{exprs_values}_{feature_control}".format(**locals())])
        # "total_{exprs_values}" not formatted yet
        obs_metrics["pct_{exprs_values}_{feature_control}".format(**locals())] = \
            obs_metrics["total_{exprs_values}_{feature_control}".format(**locals())] / \
            obs_metrics["total_{exprs_values}"] * 100
    # Calculate var metrics
    var_metrics["mean_{exprs_values}"] = np.ravel(X.mean(axis=0))
    var_metrics["log1p_mean_{exprs_values}"] = np.log1p(
        var_metrics["mean_{exprs_values}"])
    var_metrics["n_cells_by_{exprs_values}"] = np.ravel((X != 0).sum(axis=0))
    var_metrics["pct_dropout_by_{exprs_values}"] = \
        (1 - var_metrics["n_cells_by_{exprs_values}"] / X.shape[0]) * 100
    var_metrics["total_{exprs_values}"] = np.ravel(X.sum(axis=0))
    var_metrics["log1p_total_{exprs_values}"] = np.log1p(
        var_metrics["total_{exprs_values}"])
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
    Calculates cumulative proportions of top expressed genes

    Parameters
    ----------
    mtx : `Union[np.array, sparse.spmatrix]`
        Matrix, where each row is a sample, each column a feature.
    n : `int`
        Rank to calculate proportions up to. Value is treated as 1-indexed
        `n=50` will calculate cumulative proportions up to the 50th most
        expressed gene.
    """
    if issparse(mtx):
        if not isspmatrix_csr(mtx):
            mtx = csr_matrix(mtx)
        # Allowing numba to do more
        return top_proportions_sparse_csr(mtx.data, mtx.indptr, n)
    else:
        return top_proportions_dense(mtx, n)

@numba.jit
def top_proportions_dense(mtx, n):
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.argpartition, 1, -mtx, n-1)
    partitioned = partitioned[:, :n]
    values = np.zeros_like(partitioned, dtype=np.float64)
    for i in range(partitioned.shape[0]):
        vec = mtx[i, partitioned[i, :]]  # Not a view
        vec[::-1].sort()  # Sorting on a reversed view (e.g. a descending sort)
        vec = np.cumsum(vec) / sums[i]
        values[i, :] = vec
    return values

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
            total = (data[start:end]).sum()  # Is this not just vec.sum()?
        vec[::-1].sort()
        values[i, :] = vec.cumsum() / total
    return values


def top_segment_proportions(mtx, ns):
    """
    Calculates total percentage of counts in top ns genes.

    Parameters
    ----------
    mtx : `Union[np.array, sparse.spmatrix]`
        Matrix, where each row is a sample, each column a feature.
    ns : `Container[Int]`
        Positions to calculate cumulative proportion at. Values are considered
        1-indexed, e.g. `ns=[50]` will calculate cumulative proportion up to
        the 50th most expressed gene.
    """
    # Pretty much just does dispatch
    if issparse(mtx):
        if not isspmatrix_csr(mtx):
            mtx = csr_matrix(mtx)
        return top_segment_proportions_sparse_csr(mtx.data, mtx.indptr, ns)
    else:
        return top_segment_proportions_dense(mtx, ns)

def top_segment_proportions_dense(mtx, ns):
    # Currently ns is considered to be 1 indexed
    ns = np.sort(ns)
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(
        np.partition, 1, mtx, mtx.shape[1] - ns)[:, ::-1][:, :ns[-1]]
    values = np.zeros((mtx.shape[0], len(ns)))
    acc = np.zeros((mtx.shape[0]))
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]

@numba.jit(parallel=True)
def top_segment_proportions_sparse_csr(data, indptr, ns):
    ns = np.sort(ns)
    maxidx = ns[-1]
    sums = np.zeros((indptr.size - 1), dtype=data.dtype)
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
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]

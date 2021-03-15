"""Moran's I global spatial autocorrelation."""
from typing import Any, Union, Optional
from functools import singledispatch
from anndata import AnnData

import numpy as np
import pandas as pd
from scipy import sparse

import numba.types as nt
from numba import njit

from scanpy.get import _get_obs_rep

it = nt.int64
ft = nt.float64
ip = np.int64
fp = np.float64
tt = nt.UniTuple


@singledispatch
def morans_i(
    adata: AnnData,
    *,
    vals: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    use_graph: Optional[str] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
    obsp: Optional[str] = None,
    use_raw: bool = False,
) -> Union[np.ndarray, float]:
    r"""
    Calculate Moran’s I Global Autocorrelation Statistic.

    Moran’s I is a global autocorrelation statistic for some measure on a graph. It is commonly used in
    spatial data analysis to assess autocorrelation on a 2D grid. It is closely related to Geary's C,
    but not identical. More info can be found `here`<https://en.wikipedia.org/wiki/Moran%27s_I> .

    .. math::

        I=\frac{n}{S_{0}} \frac{\sum_{i=1}^{n} \sum_{j=1}^{n} w_{i, j} z_{i} z_{j}}{\sum_{i=1}^{n} z_{i}^{2}}

    Params
    ------
    adata
    vals
        Values to calculate Moran's I for. If this is two dimensional, should
        be of shape `(n_features, n_cells)`. Otherwise should be of shape
        `(n_cells,)`. This matrix can be selected from elements of the anndata
        object by using key word arguments: `layer`, `obsm`, `obsp`, or
        `use_raw`.
    use_graph
        Key to use for graph in anndata object. If not provided, default
        neighbors connectivities will be used instead.
    layer
        Key for `adata.layers` to choose `vals`.
    obsm
        Key for `adata.obsm` to choose `vals`.
    obsp
        Key for `adata.obsp` to choose `vals`.
    use_raw
        Whether to use `adata.raw.X` for `vals`.


    This function can also be called on the graph and values directly. In this case
    the signature looks like:

    Params
    ------
    g
        The graph
    vals
        The values


    See the examples for more info.

    Returns
    -------
    If vals is two dimensional, returns a 1 dimensional ndarray array. Returns
    a scalar if `vals` is 1d.


    Examples
    --------

    Calculate Morans I for each components of a dimensionality reduction:

    .. code:: python

        import scanpy as sc, numpy as np

        pbmc = sc.datasets.pbmc68k_processed()
        pc_c = sc.metrics.morans_i(pbmc, obsm="X_pca")


    It's equivalent to call the function directly on the underlying arrays:

    .. code:: python

        alt = sc.metrics.morans_i(pbmc.obsp["connectivities"], pbmc.obsm["X_pca"].T)
        np.testing.assert_array_equal(pc_c, alt)
    """
    if use_graph is None:
        # Fix for anndata<0.7
        if hasattr(adata, "obsp") and "connectivities" in adata.obsp:
            g = adata.obsp["connectivities"]
        elif "neighbors" in adata.uns:
            g = adata.uns["neighbors"]["connectivities"]
        else:
            raise ValueError("Must run neighbors first.")
    else:
        raise NotImplementedError()
    if vals is None:
        vals = _get_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return morans_i(g, vals)


@njit(
    ft(ft[:], it[:], it[:], ft[:], ft[:], ft, ft),
    cache=True,
    parallel=True,
    fastmath=True,
)
def _morans_i_vec_W(
    x: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    z: np.ndarray,
    W: fp,
    z2ss: fp,
) -> Any:

    zl = np.empty(len(z))
    N = len(indptr) - 1

    for i in range(N):
        s = slice(indptr[i], indptr[i + 1])
        i_indices = indices[s]
        i_data = data[s]
        zl[i] = np.sum(i_data * z[i_indices])
    inum = (z * zl).sum()

    return len(x) / W * inum / z2ss


@njit(ft(ft[:], it[:], it[:], ft[:]), cache=True, parallel=True, fastmath=True)
def _morans_i_vec(
    x: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
) -> Any:
    W = data.sum()
    z = x - x.mean()
    z2ss = (z * z).sum()
    return _morans_i_vec_W(x, indptr, indices, data, z, W, z2ss)


@njit(ft[:](ft[:, :], it[:], it[:], ft[:]), cache=True, parallel=True, fastmath=True)
def _morans_i_mtx(
    X: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
) -> Any:
    M, N = X.shape
    assert N == len(indptr) - 1
    W = data.sum()
    out = np.zeros(M, dtype=ft)
    for k in range(M):
        x = X[k, :].astype(ft)
        z = x - x.mean()
        z2ss = (z * z).sum()
        out[k] = _morans_i_vec_W(x, indptr, indices, data, z, W, z2ss)
    return out


@njit(
    ft[:](ft[:], it[:], it[:], it[:], it[:], ft[:], tt(it, 2)),
    cache=True,
    parallel=True,
    fastmath=True,
)
def _morans_i_mtx_csr(
    X_data: np.ndarray,
    X_indices: np.ndarray,
    X_indptr: np.ndarray,
    indices: np.ndarray,
    indptr: np.ndarray,
    data: np.ndarray,
    X_shape: tuple,
):
    M, N = X_shape
    W = data.sum()
    out = np.zeros(M, dtype=ft)
    x_data_list = np.split(X_data, X_indptr[1:-1])
    x_indices_list = np.split(X_indices, X_indptr[1:-1])
    for k in range(M):
        x = np.zeros(N, dtype=ft)
        x_index = x_indices_list[k]
        x[x_index] = x_data_list[k]
        z = x - x.mean()
        z2ss = (z * z).sum()
        out[k] = _morans_i_vec_W(x, indptr, indices, data, z, W, z2ss)
    return out


###############################################################################
# Interface (taken from gearys C)
###############################################################################
@singledispatch
def _resolve_vals(val):
    return np.asarray(val)


@_resolve_vals.register(np.ndarray)
@_resolve_vals.register(sparse.csr_matrix)
def _(val):
    return val


@_resolve_vals.register(sparse.spmatrix)
def _(val):
    return sparse.csr_matrix(val)


@_resolve_vals.register(pd.DataFrame)
@_resolve_vals.register(pd.Series)
def _(val):
    return val.to_numpy()


@morans_i.register(sparse.csr_matrix)
def _morans_i(g, vals) -> np.ndarray:
    assert g.shape[0] == g.shape[1], "`g` should be a square adjacency matrix"
    vals = _resolve_vals(vals)
    g_data = g.data.astype(fp, copy=False)
    if isinstance(vals, sparse.csr_matrix):
        assert g.shape[0] == vals.shape[1]
        return _morans_i_mtx_csr(
            vals.data.astype(fp, copy=False),
            vals.indices.astype(ip),
            vals.indptr.astype(ip),
            g.indices.astype(ip),
            g.indptr.astype(ip),
            g_data,
            vals.shape,
        )
    elif isinstance(vals, np.ndarray) and vals.ndim == 1:
        assert g.shape[0] == vals.shape[0]
        return _morans_i_vec(vals, g.indptr, g.indices, g_data)
    elif isinstance(vals, np.ndarray) and vals.ndim == 2:
        assert g.shape[0] == vals.shape[1]
        return _morans_i_mtx(
            vals.astype(fp), g.indptr.astype(ip), g.indices.astype(ip), g_data
        )
    else:
        raise NotImplementedError()


# @njit(
#     ft[:](ft[:], it[:], it[:], ft[:], it),
#     parallel=False,
#     fastmath=True,
# )
# def _moran_score_perms(
#     counts: np.ndarray,
#     indptr: np.ndarray,
#     indices: np.ndarray,
#     data: np.ndarray,
#     n_perms: ip,
# ) -> np.ndarray:

#     perms = np.empty(n_perms, dtype=ft)
#     res = np.empty(3, dtype=ft)

#     z = counts - counts.mean()
#     data_sum = data.sum()
#     z2ss = (z * z).sum()

#     res[0] = _compute_moran(counts, indptr, indices, data, z, data_sum, z2ss)

#     for p in range(len(perms)):
#         np.random.shuffle(z)
#         perms[p] = _compute_moran(counts, indptr, indices, data, z, data_sum, z2ss)
#     res[1] = (np.sum(perms > res[0]) + 1) / (n_perms + 1)
#     res[2] = np.var(perms)
#     return res

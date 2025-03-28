"""Moran's I global spatial autocorrelation."""

from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING

import numba
import numpy as np
from scipy import sparse

from .._compat import fullname, njit
from ..get import _get_obs_rep
from ._common import _check_vals, _resolve_vals

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import DaskArray


@singledispatch
def morans_i(
    adata: AnnData,
    *,
    vals: NDArray | sparse.spmatrix | DaskArray | None = None,
    use_graph: str | None = None,
    layer: str | None = None,
    obsm: str | None = None,
    obsp: str | None = None,
    use_raw: bool = False,
) -> np.ndarray | float:
    r"""Calculate Moran’s I Global Autocorrelation Statistic.

    Moran’s I is a global autocorrelation statistic for some measure on a graph. It is commonly used in
    spatial data analysis to assess autocorrelation on a 2D grid. It is closely related to Geary's C,
    but not identical. More info can be found `here <https://en.wikipedia.org/wiki/Moran%27s_I>`_.

    .. math::

        I =
            \frac{
                N \sum_{i, j} w_{i, j} z_{i} z_{j}
            }{
                S_{0} \sum_{i} z_{i}^{2}
            }

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
    Calculate Moran’s I for each components of a dimensionality reduction:

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
            msg = "Must run neighbors first."
            raise ValueError(msg)
    else:
        raise NotImplementedError()
    if vals is None:
        vals = _get_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return morans_i(g, vals)


###############################################################################
# Calculation
###############################################################################
# This is done in a very similar way to gearys_c. See notes there for details.


@njit
def _morans_i_vec(
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x: np.ndarray,
) -> float:
    W = g_data.sum()
    return _morans_i_vec_W(g_data, g_indices, g_indptr, x, W)


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _morans_i_vec_W(
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x: np.ndarray,
    W: np.float64,
) -> float:
    z = x - x.mean()
    z2ss = (z * z).sum()
    n = len(x)
    inum = 0.0

    for i in numba.prange(n):
        s = slice(g_indptr[i], g_indptr[i + 1])
        i_indices = g_indices[s]
        i_data = g_data[s]
        inum += (i_data * z[i_indices]).sum() * z[i]

    return len(x) / W * inum / z2ss


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _morans_i_vec_W_sparse(  # noqa: PLR0917
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x_data: np.ndarray,
    x_indices: np.ndarray,
    n: int,
    W: np.float64,
) -> float:
    x = np.zeros(n, dtype=x_data.dtype)
    x[x_indices] = x_data
    return _morans_i_vec_W(g_data, g_indices, g_indptr, x, W)


@njit
def _morans_i_mtx(
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    X: np.ndarray,
) -> np.ndarray:
    m, n = X.shape
    assert n == len(g_indptr) - 1
    W = g_data.sum()
    out = np.zeros(m, dtype=np.float64)
    for k in numba.prange(m):
        x = X[k, :]
        out[k] = _morans_i_vec_W(g_data, g_indices, g_indptr, x, W)
    return out


@njit
def _morans_i_mtx_csr(  # noqa: PLR0917
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x_data: np.ndarray,
    x_indices: np.ndarray,
    x_indptr: np.ndarray,
    x_shape: tuple,
) -> np.ndarray:
    m, n = x_shape
    W = g_data.sum()
    out = np.zeros(m, dtype=np.float64)
    x_data_list = np.split(x_data, x_indptr[1:-1])
    x_indices_list = np.split(x_indices, x_indptr[1:-1])
    for k in numba.prange(m):
        out[k] = _morans_i_vec_W_sparse(
            g_data,
            g_indices,
            g_indptr,
            x_data_list[k],
            x_indices_list[k],
            n,
            W,
        )
    return out


###############################################################################
# Interface (taken from gearys C)
###############################################################################


@morans_i.register(sparse.csr_matrix)
def _morans_i(
    g: sparse.csr_matrix, vals: NDArray | sparse.spmatrix | DaskArray
) -> np.ndarray:
    assert g.shape[0] == g.shape[1], "`g` should be a square adjacency matrix"
    vals = _resolve_vals(vals)
    g_data = g.data.astype(np.float64, copy=False)
    if isinstance(vals, sparse.csr_matrix):
        assert g.shape[0] == vals.shape[1]
        new_vals, idxer, full_result = _check_vals(vals)
        result = _morans_i_mtx_csr(
            g_data,
            g.indices,
            g.indptr,
            new_vals.data.astype(np.float64, copy=False),
            new_vals.indices,
            new_vals.indptr,
            new_vals.shape,
        )
        full_result[idxer] = result
        return full_result
    elif isinstance(vals, np.ndarray) and vals.ndim == 1:
        assert g.shape[0] == vals.shape[0]
        return _morans_i_vec(g_data, g.indices, g.indptr, vals)
    elif isinstance(vals, np.ndarray) and vals.ndim == 2:
        assert g.shape[0] == vals.shape[1]
        new_vals, idxer, full_result = _check_vals(vals)
        result = _morans_i_mtx(
            g_data, g.indices, g.indptr, new_vals.astype(np.float64, copy=False)
        )
        full_result[idxer] = result
        return full_result
    else:
        msg = (
            "Moran’s I metric not implemented for vals of type "
            f"{fullname(type(vals))} and ndim {vals.ndim}."
        )
        raise NotImplementedError(msg)

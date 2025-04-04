"""Geary's C autocorrelation."""

from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING, cast

import numba
import numpy as np

from .._compat import CSRBase, _register_union, njit
from ..get import _get_obs_rep
from ._common import _get_graph, _SparseMetric

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from ._common import _Vals


@singledispatch
def gearys_c(
    adata_or_graph: AnnData | CSRBase,
    /,
    vals: _Vals | None = None,
    *,
    use_graph: str | None = None,
    layer: str | None = None,
    obsm: str | None = None,
    obsp: str | None = None,
    use_raw: bool = False,
) -> np.ndarray | float:
    r"""Calculate `Geary's C <https://en.wikipedia.org/wiki/Geary's_C>`_.

    Specifically as used by `VISION <https://doi.org/10.1038/s41467-019-12235-0>`_.

    Geary's C is a measure of autocorrelation for some measure on a graph. This
    can be to whether measures are correlated between neighboring cells. Lower
    values indicate greater correlation.

    .. math::

        C =
        \frac{
            (N - 1)\sum_{i,j} w_{i,j} (x_i - x_j)^2
        }{
            2W \sum_i (x_i - \bar{x})^2
        }

    Params
    ------
    adata_or_graph
        AnnData object containing a graph (see ``use_graph``) or the graph itself.
        See the examples for more info.
    vals
        Values to calculate Geary's C for. If this is two dimensional, should
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

    Returns
    -------
    If vals is two dimensional, returns a 1 dimensional ndarray array. Returns
    a scalar if `vals` is 1d.

    Examples
    --------
    Calculate Geary’s C for each components of a dimensionality reduction:

    .. code:: python

        import scanpy as sc, numpy as np

        pbmc = sc.datasets.pbmc68k_processed()
        pc_c = sc.metrics.gearys_c(pbmc, obsm="X_pca")

    It's equivalent to call the function directly on the underlying arrays:

    .. code:: python

        alt = sc.metrics.gearys_c(pbmc.obsp["connectivities"], pbmc.obsm["X_pca"].T)
        np.testing.assert_array_equal(pc_c, alt)

    """
    adata = cast("AnnData", adata_or_graph)
    g = _get_graph(adata, use_graph=use_graph)
    if vals is None:
        vals = _get_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return gearys_c(g, vals)


@_register_union(gearys_c, CSRBase)
def _gearys_c(graph: CSRBase, /, vals: _Vals) -> NDArray:
    return _GearysC(graph, vals)()


class _GearysC(_SparseMetric):
    name = "Geary’s C"

    def mtx(self, vals_het: NDArray | CSRBase, /) -> NDArray:
        g_parts = (self.graph.data, self.graph.indices, self.graph.indptr)
        if isinstance(vals_het, np.ndarray):
            return _gearys_c_mtx(*g_parts, vals_het)
        v_parts = (vals_het.data, vals_het.indices, vals_het.indptr)
        return _gearys_c_mtx_csr(*g_parts, *v_parts, vals_het.shape)

    def vec(self) -> np.float64:
        W = self.graph.data.sum()
        g_parts = (self.graph.data, self.graph.indices, self.graph.indptr)
        return _gearys_c_vec_W(*g_parts, self._vals, W)


###############################################################################
# Calculation
###############################################################################
# Some notes on the implementation:
# * This could be phrased as tensor multiplication. However that does not get
#   parallelized, which boosts performance almost linearly with cores.
# * Due to the umap setting the default threading backend, a parallel numba
#   function that calls another parallel numba function can get stuck. This
#   ends up meaning code re-use will be limited until umap 0.4.
#   See: https://github.com/lmcinnes/umap/issues/306
# * There can be a fair amount of numerical instability here (big reductions),
#   so data is cast to float64. Removing these casts/ conversion will cause the
#   tests to fail.


@njit
def _gearys_c_vec_W(
    data: np.ndarray,
    indices: np.ndarray,
    indptr: np.ndarray,
    x: np.ndarray,
    W: np.float64,
) -> np.float64:
    n = len(indptr) - 1
    x = x.astype(np.float64)
    x_bar = x.mean()

    total = 0.0
    for i in numba.prange(n):
        s = slice(indptr[i], indptr[i + 1])
        i_indices = indices[s]
        i_data = data[s]
        total += np.sum(i_data * ((x[i] - x[i_indices]) ** 2))

    numer = (n - 1) * total
    denom = 2 * W * ((x - x_bar) ** 2).sum()
    return numer / denom


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Inner functions (per element C)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For calling gearys_c on collections.
# TODO: These are faster if we can compile them in parallel mode. However,
# `workqueue` does not allow nested functions to be parallelized.
# Additionally, there are currently problems with numba's compiler around
# parallelization of this code:
# https://github.com/numba/numba/issues/6774#issuecomment-788789663


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _gearys_c_inner_sparse_x_densevec(
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x: np.ndarray,
    W: np.float64,
) -> np.float64:
    x_bar = x.mean()
    total = 0.0
    n = len(x)
    for i in numba.prange(n):
        s = slice(g_indptr[i], g_indptr[i + 1])
        i_indices = g_indices[s]
        i_data = g_data[s]
        total += np.sum(i_data * ((x[i] - x[i_indices]) ** 2))
    numer = (n - 1) * total
    denom = 2 * W * ((x - x_bar) ** 2).sum()
    return numer / denom


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _gearys_c_inner_sparse_x_sparsevec(  # noqa: PLR0917
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x_data: np.ndarray,
    x_indices: np.ndarray,
    n: int,
    W: np.float64,
) -> np.float64:
    x = np.zeros(n, dtype=np.float64)
    x[x_indices] = x_data
    x_bar = np.sum(x_data) / n
    total = 0.0
    n = len(x)
    for i in numba.prange(n):
        s = slice(g_indptr[i], g_indptr[i + 1])
        i_indices = g_indices[s]
        i_data = g_data[s]
        total += np.sum(i_data * ((x[i] - x[i_indices]) ** 2))
    numer = (n - 1) * total
    # Expanded from 2 * W * ((x_k - x_k_bar) ** 2).sum(), but uses sparsity
    # to skip some calculations
    # fmt: off
    denom = (
        2 * W
        * (
            np.sum(x_data ** 2)
            - np.sum(x_data * x_bar * 2)
            + (x_bar ** 2) * n
        )
    )
    # fmt: on
    return numer / denom


@njit
def _gearys_c_mtx(
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
        x = X[k, :].astype(np.float64)
        out[k] = _gearys_c_inner_sparse_x_densevec(g_data, g_indices, g_indptr, x, W)
    return out


@njit
def _gearys_c_mtx_csr(  # noqa: PLR0917
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
        out[k] = _gearys_c_inner_sparse_x_sparsevec(
            g_data,
            g_indices,
            g_indptr,
            x_data_list[k],
            x_indices_list[k],
            n,
            W,
        )
    return out

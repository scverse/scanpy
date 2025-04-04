"""Moran's I global spatial autocorrelation."""

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
def morans_i(
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
    adata_or_graph
        AnnData object containing a graph (see ``use_graph``) or the graph itself.
        See the examples for more info.
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
    adata = cast("AnnData", adata_or_graph)
    g = _get_graph(adata, use_graph=use_graph)
    if vals is None:
        vals = _get_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return morans_i(g, vals)


@_register_union(morans_i, CSRBase)
def _morans_i(graph: CSRBase, /, vals: _Vals) -> NDArray:
    return _MoransI(graph, vals)()


class _MoransI(_SparseMetric):
    name = "Moran’s I"

    def mtx(self, vals_het: NDArray | CSRBase, /) -> NDArray:
        g_parts = (self.graph.data, self.graph.indices, self.graph.indptr)
        if isinstance(vals_het, np.ndarray):
            return _morans_i_mtx(*g_parts, vals_het)
        v_parts = (vals_het.data, vals_het.indices, vals_het.indptr)
        return _morans_i_mtx_csr(*g_parts, *v_parts, vals_het.shape)

    def vec(self) -> np.float64:
        W = self.graph.data.sum()
        g_parts = (self.graph.data, self.graph.indices, self.graph.indptr)
        return _morans_i_vec_W(*g_parts, self._vals, W)


###############################################################################
# Calculation
###############################################################################
# This is done in a very similar way to gearys_c. See notes there for details.


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _morans_i_vec_W(
    g_data: np.ndarray,
    g_indices: np.ndarray,
    g_indptr: np.ndarray,
    x: np.ndarray,
    W: np.float64,
) -> np.float64:
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
) -> np.float64:
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

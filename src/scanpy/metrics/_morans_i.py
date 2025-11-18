"""Moran's I global spatial autocorrelation."""

from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING, cast

import numba
import numpy as np

from .._compat import CSRBase, njit
from .._utils import _doc_params
from ..get import _get_obs_rep
from ..neighbors._doc import doc_neighbors_key
from ._common import _get_graph, _SparseMetric

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from ._common import _Vals


@singledispatch
@_doc_params(neighbors_key=doc_neighbors_key)
def morans_i(
    adata_or_graph: AnnData | CSRBase,
    /,
    vals: _Vals | None = None,
    *,
    use_graph: str | None = None,
    neighbors_key: str | None = None,
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
        Key to use for graph in anndata object.
        If not provided, default neighbors connectivities will be used instead.
        (See ``neighbors_key`` below.)
    {neighbors_key}
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
    g = _get_graph(adata, use_graph=use_graph, neighbors_key=neighbors_key)
    if vals is None:
        vals = _get_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return morans_i(g, vals)


@morans_i.register(CSRBase)
def _morans_i(graph: CSRBase, /, vals: _Vals) -> NDArray:
    return _MoransI(graph, vals)()


class _MoransI(_SparseMetric):
    name = "Moran’s I"

    def mtx(self, vals_het: NDArray | CSRBase, /) -> NDArray:
        if isinstance(vals_het, np.ndarray):
            return _morans_i_mtx(self.graph, vals_het)
        return _morans_i_mtx_csr(self.graph, vals_het)

    def vec(self) -> np.float64:
        w = self.graph.data.sum()
        return _morans_i_vec_w(self.graph, self._vals, w)


###############################################################################
# Calculation
###############################################################################
# This is done in a very similar way to gearys_c. See notes there for details.


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _morans_i_vec_w(g: CSRBase, x: np.ndarray, w: np.float64) -> np.float64:
    z = x - x.mean()
    z2ss = (z * z).sum()
    n = len(x)
    inum = 0.0

    for i in numba.prange(n):
        s = slice(g.indptr[i], g.indptr[i + 1])
        i_indices = g.indices[s]
        i_data = g.data[s]
        inum += (i_data * z[i_indices]).sum() * z[i]

    return len(x) / w * inum / z2ss


@numba.njit(cache=True, parallel=False)  # noqa: TID251
def _morans_i_vec_w_sparse(
    g: CSRBase, x_data: np.ndarray, x_indices: np.ndarray, n: int, w: np.float64
) -> np.float64:
    x_vec = np.zeros(n, dtype=x_data.dtype)
    x_vec[x_indices] = x_data
    return _morans_i_vec_w(g, x_vec, w)


@njit
def _morans_i_mtx(g: CSRBase, x: np.ndarray) -> np.ndarray:
    m, n = x.shape
    assert n == len(g.indptr) - 1
    w = g.data.sum()
    out = np.zeros(m, dtype=np.float64)
    for k in numba.prange(m):
        x_vec = x[k, :]
        out[k] = _morans_i_vec_w(g, x_vec, w)
    return out


@njit
def _morans_i_mtx_csr(g: CSRBase, x: CSRBase) -> np.ndarray:
    m, n = x.shape
    w = g.data.sum()
    out = np.zeros(m, dtype=np.float64)
    x_data_list = np.split(x.data, x.indptr[1:-1])
    x_indices_list = np.split(x.indices, x.indptr[1:-1])
    for k in numba.prange(m):
        out[k] = _morans_i_vec_w_sparse(g, x_data_list[k], x_indices_list[k], n, w)
    return out

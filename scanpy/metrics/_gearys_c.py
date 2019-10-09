# TODO: Calling this after UMAP has imported will hang, pending https://github.com/lmcinnes/umap/issues/306
from typing import Optional, Union

from anndata import AnnData
from multipledispatch import dispatch
import numba
import numpy as np
import pandas as pd
from scipy import sparse


def _choose_obs_rep(adata, *, use_raw=False, layer=None, obsm=None, obsp=None):
    """
    Choose array aligned with obs annotation.
    """
    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made <= 1
    if choices_made == 0:
        return adata.X
    elif is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    elif is_obsm:
        return adata.obsm[obsm]
    elif is_obsp:
        return adata.obsp[obsp]
    else:
        raise RuntimeError("You broke it. But how? Please report this.")


###############################################################################
# Calculation
###############################################################################


@numba.njit(cache=True, parallel=True)
def _gearys_c_vec(data, indices, indptr, x):
    W = data.sum()
    return _gearys_c_vec_W(data, indices, indptr, x, W)


@numba.njit(cache=True, parallel=True)
def _gearys_c_vec_W(data, indices, indptr, x, W):
    N = len(indptr) - 1
    x_bar = x.mean()

    total = 0.0
    for i in numba.prange(N):
        s = slice(indptr[i], indptr[i + 1])
        i_indices = indices[s]
        i_data = data[s]
        total += np.sum(i_data * ((x[i] - x[i_indices]) ** 2))

    numer = (N - 1) * total
    denom = 2 * W * ((x - x_bar) ** 2).sum()
    C = numer / denom

    return C


@numba.njit(cache=True, parallel=True)
def _gearys_c_mtx_csr(
    g_data, g_indices, g_indptr, x_data, x_indices, x_indptr, x_shape
):
    M, N = x_shape
    W = g_data.sum()
    out = np.zeros(M, dtype=np.float64)
    for k in numba.prange(M):
        x_arr = np.zeros(N, dtype=x_data.dtype)
        sk = slice(x_indptr[k], x_indptr[k + 1])
        x_arr[x_indices[sk]] = x_data[sk]
        outval = _gearys_c_vec_W(g_data, g_indices, g_indptr, x_arr, W)
        out[k] = outval
    return out


@numba.njit(cache=True, parallel=True)
def _gearys_c_mtx(g_data, g_indices, g_indptr, X):
    M, N = X.shape
    W = g_data.sum()
    out = np.zeros(M, dtype=np.float64)
    for k in numba.prange(M):
        outval = _gearys_c_vec_W(g_data, g_indices, g_indptr, X[k, :], W)
        out[k] = outval
    return out


###############################################################################
# Interface
###############################################################################


@dispatch(sparse.csr_matrix, sparse.csr_matrix)
def gearys_c(g, vals) -> np.ndarray:
    assert g.shape[0] == g.shape[1]
    assert g.shape[0] == vals.shape[1]
    return _gearys_c_mtx_csr(
        g.data,
        g.indices,
        g.indptr,
        vals.data,
        vals.indices,
        vals.indptr,
        vals.shape,
    )


@dispatch(sparse.spmatrix, np.ndarray)
def gearys_c(g, vals):
    assert g.shape[0] == g.shape[1], "`g` should be a square matrix."
    if not isinstance(g, sparse.csr_matrix):
        g = g.tocsr()
    if vals.ndim == 1:
        assert g.shape[0] == vals.shape[0]
        return _gearys_c_vec(g.data, g.indices, g.indptr, vals)
    elif vals.ndim == 2:
        assert g.shape[0] == vals.shape[1]
        return _gearys_c_mtx(g.data, g.indices, g.indptr, vals)
    else:
        raise ValueError()


@dispatch(sparse.spmatrix, (pd.DataFrame, pd.Series))
def gearys_c(g, vals):
    return gearys_c(g, vals.values)

@dispatch(sparse.spmatrix, sparse.spmatrix)
def gearys_c(g, vals) -> np.ndarray:
    if not isinstance(g, sparse.csr_matrix):
        g = g.tocsr()
    if not isinstance(vals, sparse.csc_matrix):
        vals = vals.tocsr()
    return gearys_c(g, vals)


# TODO: Document better
# TODO: Have scanpydoc work with multipledispatch
@dispatch(AnnData)
def gearys_c(
    adata: AnnData,
    *,
    vals: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    use_graph: Optional[str] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
    obsp: Optional[str] = None,
    use_raw: bool = False,
) -> Union[np.ndarray, float]:
    """
    Calculate `Geary's C` <https://en.wikipedia.org/wiki/Geary's_C>`_, as used
    by `VISION <https://doi.org/10.1038/s41467-019-12235-0>`_.

    Geary's C is a measure of autocorrelation for some measure on a graph. This
    can be to whether measures are correlated between neighboring cells. Lower 
    values indicate greater correlation.

    ..math

        C = 
        \frac{
            (N - 1)\sum_{i,j} w_{i,j} (x_i - x_j)^2
        }{
            2W \sum_i (x_i - \bar{x})^2
        }

    Params
    ------
    adata
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
    If vals is two dimensional, returns a 1 dimensional ndarray array. Returns a scalar if `vals` is 1d.
    """
    if use_graph is None:
        if "connectivities" in adata.obsp:
            g = adata.obsp["connectivities"]
        elif "neighbors" in adata.uns:
            g = adata.uns["neighbors"]["connectivities"]
        else:
            raise ValueError("Must run neighbors first.")
    else:
        raise NotImplementedError()
    if vals is None:
        vals = _choose_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp).T
    return gearys_c(g, vals)

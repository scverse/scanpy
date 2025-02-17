from __future__ import annotations

import warnings
from functools import singledispatch
from operator import truediv
from typing import TYPE_CHECKING

import numba
import numpy as np
from anndata import AnnData
from scipy.sparse import csc_matrix, csr_matrix, issparse

from .. import logging as logg
from .._compat import DaskArray, njit, old_positionals
from .._utils import (
    _check_array_function_arguments,
    axis_mul_or_truediv,
    raise_not_implemented_error_if_backed_type,
    renamed_arg,
    view_to_actual,
)
from ..get import _check_mask, _get_obs_rep, _set_obs_rep
from ._utils import _get_mean_var

# install dask if available
try:
    import dask.array as da
except ImportError:
    da = None

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

    from .._utils import _CSMatrix


@singledispatch
def clip(x: ArrayLike, *, max_value: float, zero_center: bool = True) -> None:
    return clip_array(x, max_value=max_value, zero_center=zero_center)


@clip.register(csr_matrix)
@clip.register(csc_matrix)
def _(x: _CSMatrix, *, max_value: float, zero_center: bool = True) -> _CSMatrix:
    x.data = clip(x.data, max_value=max_value, zero_center=zero_center)
    return x


@clip.register(DaskArray)
def _(x: DaskArray, *, max_value: float, zero_center: bool = True):
    return x.map_blocks(
        clip, max_value=max_value, zero_center=zero_center, dtype=x.dtype, meta=x._meta
    )


@njit
def clip_array(
    X: NDArray[np.floating], *, max_value: float, zero_center: bool
) -> NDArray[np.floating]:
    a_min, a_max = -max_value, max_value
    if X.ndim > 1:
        for r, c in numba.pndindex(X.shape):
            if X[r, c] > a_max:
                X[r, c] = a_max
            elif X[r, c] < a_min and zero_center:
                X[r, c] = a_min
    else:
        for i in numba.prange(X.size):
            if X[i] > a_max:
                X[i] = a_max
            elif X[i] < a_min and zero_center:
                X[i] = a_min
    return X


@renamed_arg("X", "data", pos_0=True)
@old_positionals("zero_center", "max_value", "copy", "layer", "obsm")
@singledispatch
def scale(
    data: AnnData | _CSMatrix | np.ndarray | DaskArray,
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    copy: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    mask_obs: NDArray[np.bool_] | str | None = None,
) -> AnnData | _CSMatrix | np.ndarray | DaskArray | None:
    """\
    Scale data to unit variance and zero mean.

    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and (for zero_center==True) set to 0
        during this operation. In the future, they might be set to NaNs.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        Whether this function should be performed inplace. If an AnnData object
        is passed, this also determines if a copy is returned.
    layer
        If provided, which element of layers to scale.
    obsm
        If provided, which element of obsm to scale.
    mask_obs
        Restrict both the derivation of scaling parameters and the scaling itself
        to a certain set of observations. The mask is specified as a boolean array
        or a string referring to an array in :attr:`~anndata.AnnData.obs`.
        This will transform data from csc to csr format if `issparse(data)`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an updated `AnnData` object. Sets the following fields:

    `adata.X` | `adata.layers[layer]` : :class:`numpy.ndarray` | :class:`scipy.sparse._csr.csr_matrix` (dtype `float`)
        Scaled count data matrix.
    `adata.var['mean']` : :class:`pandas.Series` (dtype `float`)
        Means per gene before scaling.
    `adata.var['std']` : :class:`pandas.Series` (dtype `float`)
        Standard deviations per gene before scaling.
    `adata.var['var']` : :class:`pandas.Series` (dtype `float`)
        Variances per gene before scaling.
    """
    _check_array_function_arguments(layer=layer, obsm=obsm)
    if layer is not None:
        msg = f"`layer` argument inappropriate for value of type {type(data)}"
        raise ValueError(msg)
    if obsm is not None:
        msg = f"`obsm` argument inappropriate for value of type {type(data)}"
        raise ValueError(msg)
    return scale_array(
        data, zero_center=zero_center, max_value=max_value, copy=copy, mask_obs=mask_obs
    )


@scale.register(np.ndarray)
@scale.register(DaskArray)
@scale.register(csc_matrix)
@scale.register(csr_matrix)
def scale_array(
    x: np.ndarray | DaskArray | _CSMatrix,
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    copy: bool = False,
    return_mean_std: bool = False,
    mask_obs: NDArray[np.bool_] | None = None,
) -> (
    np.ndarray
    | DaskArray
    | _CSMatrix
    | tuple[
        np.ndarray | DaskArray | _CSMatrix,
        NDArray[np.float64] | DaskArray,
        NDArray[np.float64],
    ]
):
    if copy:
        x = x.copy()

    if not zero_center and max_value is not None:
        logg.info(  # Be careful of what? This should be more specific
            "... be careful when using `max_value` without `zero_center`."
        )

    if np.issubdtype(x.dtype, np.integer):
        logg.info(
            "... as scaling leads to float results, integer "
            "input is cast to float, returning copy."
        )
        x = x.astype(float)

    mask_obs = _check_mask(x, mask_obs, "obs")
    if mask_obs is not None:
        return scale_array_masked(
            x,
            mask_obs,
            zero_center=zero_center,
            max_value=max_value,
            return_mean_std=return_mean_std,
        )

    mean, var = _get_mean_var(x)
    std = np.sqrt(var)
    std[std == 0] = 1
    if zero_center:
        if isinstance(x, csr_matrix | csc_matrix) or (
            isinstance(x, DaskArray) and issparse(x._meta)
        ):
            warnings.warn(
                "zero-center being used with `DaskArray` sparse chunks. "
                "This can be bad if you have large chunks or intend to eventually read the whole data into memory.",
                UserWarning,
            )
        x -= mean

    x = axis_mul_or_truediv(
        x,
        std,
        op=truediv,
        out=x if isinstance(x, np.ndarray | csr_matrix | csc_matrix) else None,
        axis=1,
    )

    # do the clipping
    if max_value is not None:
        x = clip(x, max_value=max_value, zero_center=zero_center)
    if return_mean_std:
        return x, mean, std
    else:
        return x


def scale_array_masked(
    x: np.ndarray | DaskArray | _CSMatrix,
    mask_obs: NDArray[np.bool_],
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    return_mean_std: bool = False,
) -> (
    np.ndarray
    | DaskArray
    | _CSMatrix
    | tuple[
        np.ndarray | DaskArray | _CSMatrix,
        NDArray[np.float64] | DaskArray,
        NDArray[np.float64],
    ]
):
    if isinstance(x, csr_matrix | csc_matrix) and not zero_center:
        if isinstance(x, csc_matrix):
            x = x.tocsr()
        mean, var = _get_mean_var(x[mask_obs, :])
        std = np.sqrt(var)
        std[std == 0] = 1

        scale_and_clip_csr(
            x.indptr,
            x.indices,
            x.data,
            std=std,
            mask_obs=mask_obs,
            max_value=max_value,
        )
    else:
        x[mask_obs, :], mean, std = scale_array(
            x[mask_obs, :],
            zero_center=zero_center,
            max_value=max_value,
            return_mean_std=True,
        )

    if return_mean_std:
        return x, mean, std
    else:
        return x


@njit
def scale_and_clip_csr(
    indptr: NDArray[np.integer],
    indices: NDArray[np.integer],
    data: NDArray[np.floating],
    *,
    std: NDArray[np.floating],
    mask_obs: NDArray[np.bool_],
    max_value: float | None,
) -> None:
    for i in numba.prange(len(indptr) - 1):
        if mask_obs[i]:
            for j in range(indptr[i], indptr[i + 1]):
                if max_value is not None:
                    data[j] = min(max_value, data[j] / std[indices[j]])
                else:
                    data[j] /= std[indices[j]]


@scale.register(AnnData)
def scale_anndata(
    adata: AnnData,
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    copy: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    mask_obs: NDArray[np.bool_] | str | None = None,
) -> AnnData | None:
    adata = adata.copy() if copy else adata
    str_mean_std = ("mean", "std")
    if mask_obs is not None:
        if isinstance(mask_obs, str):
            str_mean_std = (f"mean of {mask_obs}", f"std of {mask_obs}")
        else:
            str_mean_std = ("mean with mask", "std with mask")
        mask_obs = _check_mask(adata, mask_obs, "obs")
    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer, obsm=obsm)
    raise_not_implemented_error_if_backed_type(X, "scale")
    X, adata.var[str_mean_std[0]], adata.var[str_mean_std[1]] = scale(
        X,
        zero_center=zero_center,
        max_value=max_value,
        copy=False,  # because a copy has already been made, if it were to be made
        return_mean_std=True,
        mask_obs=mask_obs,
    )
    _set_obs_rep(adata, X, layer=layer, obsm=obsm)
    return adata if copy else None

from __future__ import annotations

from functools import singledispatch
from operator import truediv
from typing import TYPE_CHECKING

import numba
import numpy as np
from anndata import AnnData
from fast_array_utils.numba import njit
from fast_array_utils.stats import mean_var
from scverse_misc import Deprecation, deprecated_arg

from .. import logging as logg
from .._compat import CSBase, CSCBase, CSRBase, DaskArray, warn
from .._docs import DEPR_COPY, doc_mask, doc_out, doc_use
from .._settings import Default, settings
from .._utils import (
    _doc_params,
    axis_mul_or_truediv,
    check_array_function_arguments,
    dematrix,
    expose_dispatch,
    raise_not_implemented_error_if_backed_type,
    view_to_actual,
)
from ..get import _check_mask, _get_arr, _write_out
from ..get.get import AdRef, _mask_arg, _resolve_obs

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

    from ..get.get import Mask, RepAcc

type _Array = CSBase | np.ndarray | DaskArray


@singledispatch
def clip[A: _Array](
    x: ArrayLike | A, *, max_value: float, zero_center: bool = True
) -> A:
    return clip_array(x, max_value=max_value, zero_center=zero_center)


@clip.register(CSBase)
def _(x: CSBase, *, max_value: float, zero_center: bool = True) -> CSBase:
    x.data = clip(x.data, max_value=max_value, zero_center=zero_center)
    return x


@clip.register(DaskArray)
def _(x: DaskArray, *, max_value: float, zero_center: bool = True) -> DaskArray:
    return x.map_blocks(
        clip, max_value=max_value, zero_center=zero_center, dtype=x.dtype, meta=x._meta
    )


@njit
def clip_array(
    x: NDArray[np.floating], /, *, max_value: float, zero_center: bool
) -> NDArray[np.floating]:
    a_min, a_max = -max_value, max_value
    if x.ndim > 1:
        for r, c in numba.pndindex(x.shape):
            if x[r, c] > a_max:
                x[r, c] = a_max
            elif x[r, c] < a_min and zero_center:
                x[r, c] = a_min
    else:
        for i in numba.prange(x.size):
            if x[i] > a_max:
                x[i] = a_max
            elif x[i] < a_min and zero_center:
                x[i] = a_min
    return x


@singledispatch
def _scale(data, **kwargs):
    """Dispatch on array kind. `AnnData` goes to `scale_anndata`, see `scale`."""
    return scale_array(data, **kwargs)


@_doc_params(
    mask=doc_mask(
        "Restrict both the derivation of scaling parameters and the scaling itself\n"
        "    to a certain set of observations.",
        dim="obs",
        extra="This will transform data from csc to csr format if `issparse(data)`.",
    ),
    use=doc_use("Which matrix to scale."),
    out=doc_out(),
)
@expose_dispatch(_scale)
@deprecated_arg("mask_obs", Deprecation("1.13.0", "Use `mask` instead."))
@deprecated_arg("layer", Deprecation("1.13.0", "Use `use`/`out` instead."))
@deprecated_arg("obsm", Deprecation("1.13.0", "Use `use`/`out` instead."))
@deprecated_arg("copy", DEPR_COPY)
def scale[A: _Array](
    data: AnnData | A,
    *,
    zero_center: bool | Default = Default(preset=("scale", "zero_center")),
    max_value: float | None = None,
    use: RepAcc | str | None = None,
    out: RepAcc | str | bool = False,
    mask: Mask | None = None,
    # deprecated
    copy: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    mask_obs: Mask | None = None,
) -> AnnData | A | None:
    """Scale data to unit variance and zero mean.

    .. note::
        Variance and standard deviation are computed with Bessel's correction,
        i.e. dividing by ``n_obs - 1`` (``ddof=1``), matching :func:`numpy.std`
        with ``ddof=1`` rather than the numpy default of ``ddof=0`` (population
        variance). The difference is negligible for large `n_obs` but can matter for
        small datasets, where it also slightly shifts where `max_value`
        clipping takes effect.

    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and (for zero_center==True) set to 0
        during this operation. In the future, they might be set to NaNs.

    .. array-support:: pp.scale

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
        The default will be removed in scanpy 2.0.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    {use}
    {out}
    copy
        For array input, whether to leave the input unmodified.
    {mask}

    Returns
    -------
    Returns `None` if the result was written, else the scaled matrix.
    Sets the following fields:

    `adata.X` | `adata.layers[layer]` : :class:`numpy.ndarray` | :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Scaled count data matrix.
    `adata.var['mean']` : :class:`pandas.Series` (dtype `float`)
        Means per gene before scaling.
    `adata.var['std']` : :class:`pandas.Series` (dtype `float`)
        Standard deviations per gene before scaling.
    `adata.var['var']` : :class:`pandas.Series` (dtype `float`)
        Variances per gene before scaling.

    """
    # validated here, before dispatch, so array input gets these messages too
    if not isinstance(data, AnnData):
        check_array_function_arguments(layer=layer, obsm=obsm, use=use)
        if out is not False:
            msg = "Argument `out` is only valid if an AnnData object is passed."
            raise TypeError(msg)
        return _scale(
            data,
            zero_center=zero_center,
            max_value=max_value,
            copy=copy,
            mask=mask,
            mask_obs=mask_obs,
        )
    return scale_anndata(
        data,
        zero_center=zero_center,
        max_value=max_value,
        use=use,
        out=out,
        mask=mask,
        copy=copy,
        layer=layer,
        obsm=obsm,
        mask_obs=mask_obs,
    )


@_scale.register(np.ndarray)
@_scale.register(DaskArray)
@_scale.register(CSBase)
def scale_array[A: _Array](
    x: A,
    *,
    zero_center: bool | Default = Default(preset=("scale", "zero_center")),
    max_value: float | None = None,
    copy: bool = False,
    return_mean_std: bool = False,
    mask: NDArray[np.bool] | None = None,
    mask_obs: NDArray[np.bool] | None = None,
) -> (
    A
    | tuple[
        A,
        NDArray[np.float64] | DaskArray,
        NDArray[np.float64],
    ]
):
    if copy:
        x = x.copy()

    if isinstance(zero_center, Default):
        if settings.preset.scale.zero_center is None:
            msg = "scale() missing 1 required keyword argument: 'zero_center'"
            raise TypeError(msg)
        zero_center = settings.preset.scale.zero_center
    if not zero_center and max_value is not None:
        logg.info(  # Be careful of what? This should be more specific
            "... be careful when using `max_value` without `zero_center`."
        )

    if np.issubdtype(x.dtype, np.integer):
        logg.info(
            "... as scaling leads to float results, integer "
            "input is cast to float, returning copy."
        )
        x = x.astype(np.float64)

    mask = _mask_arg(mask, mask_obs, dim="obs")
    mask = (
        # For CSR matrices, default to a set mask to take the `scale_array_masked` path.
        # This is faster than the maskless `axis_mul_or_truediv` path.
        np.ones(x.shape[0], dtype=np.bool)
        if isinstance(x, CSRBase) and mask is None and not zero_center
        else _check_mask(x, mask, "obs")
    )
    if mask is not None:
        return scale_array_masked(
            x,
            mask,
            zero_center=zero_center,
            max_value=max_value,
            return_mean_std=return_mean_std,
        )

    mean, var = mean_var(x, axis=0, correction=1)
    std = np.sqrt(var)
    std[std == 0] = 1
    if zero_center:
        if isinstance(x, CSBase) or (
            isinstance(x, DaskArray) and isinstance(x._meta, CSBase)
        ):
            msg = "zero-centering a sparse array/matrix densifies it."
            warn(msg, UserWarning)
        x -= mean
        x = dematrix(x)

    x = axis_mul_or_truediv(
        x,
        std,
        op=truediv,
        out=x if isinstance(x, np.ndarray | CSBase) else None,
        axis=1,
    )

    # do the clipping
    if max_value is not None:
        x = clip(x, max_value=max_value, zero_center=zero_center)
    if return_mean_std:
        return x, mean, std
    else:
        return x


def scale_array_masked[A: _Array](
    x: A,
    mask_obs: NDArray[np.bool],
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    return_mean_std: bool = False,
) -> (
    A
    | tuple[
        A,
        NDArray[np.float64] | DaskArray,
        NDArray[np.float64],
    ]
):
    if isinstance(x, CSBase) and not zero_center:
        if isinstance(x, CSCBase):
            x = x.tocsr()
        mean, var = mean_var(x[mask_obs, :], axis=0, correction=1)
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
    mask_obs: NDArray[np.bool],
    max_value: float | None,
) -> None:
    for i in numba.prange(len(indptr) - 1):
        if mask_obs[i]:
            for j in range(indptr[i], indptr[i + 1]):
                if max_value is not None:
                    data[j] = min(max_value, data[j] / std[indices[j]])
                else:
                    data[j] /= std[indices[j]]


@_scale.register(AnnData)
def scale_anndata(
    adata: AnnData,
    *,
    zero_center: bool | Default = Default(preset=("scale", "zero_center")),
    max_value: float | None = None,
    use: RepAcc | str | None = None,
    out: RepAcc | str | bool = False,
    mask: Mask | None = None,
    copy: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    mask_obs: Mask | None = None,
) -> AnnData | _Array | None:
    if copy:
        if out is True:
            msg = "`copy=True` cannot be used with `out=True`."
            raise TypeError(msg)
        adata = adata.copy()
    use = _resolve_obs(use)
    mask = _mask_arg(mask, mask_obs, dim="obs")
    str_mean_std = ("mean", "std")
    if mask is not None:
        if isinstance(mask, str | AdRef):
            str_mean_std = (f"mean of {mask}", f"std of {mask}")
        else:
            str_mean_std = ("mean with mask", "std with mask")
        mask = _check_mask(adata, mask, "obs")
    view_to_actual(adata)
    x = _get_arr(adata, use, layer=layer, obsm=obsm)
    raise_not_implemented_error_if_backed_type(x, "scale")
    x, adata.var[str_mean_std[0]], adata.var[str_mean_std[1]] = _scale(
        x,
        zero_center=zero_center,
        max_value=max_value,
        # scaling is in-place, so leave the source alone unless we write back over it
        copy=out is not False,
        return_mean_std=True,
        mask=mask,
    )
    if (
        unwritten := _write_out(adata, x, out=out, use=use, layer=layer, obsm=obsm)
    ) is not None:
        return unwritten
    return adata if copy else None

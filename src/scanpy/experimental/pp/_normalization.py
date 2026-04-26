from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
from anndata import AnnData

from ... import logging as logg
from ... import settings
from ..._compat import CSBase, warn
from ..._settings import Default
from ..._utils import _doc_params, check_nonnegative_integers, view_to_actual
from ..._utils.random import _accepts_legacy_random_state
from ...experimental._docs import (
    doc_adata,
    doc_check_values,
    doc_copy,
    doc_dist_params,
    doc_inplace,
    doc_layer,
    doc_pca_chunk,
)
from ...get import _check_mask, _get_obs_rep, _set_obs_rep
from ...preprocessing._docs import doc_mask_var
from ...preprocessing._pca import pca

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from ..._utils.random import RNGLike, SeedLike


def _pearson_residuals(
    x: CSBase | np.ndarray, /, theta, clip, check_values, *, copy: bool = False
):
    x = x.copy() if copy else x

    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        msg = "Pearson residuals require theta > 0"
        raise ValueError(msg)
    # prepare clipping
    if clip is None:
        n = x.shape[0]
        clip = np.sqrt(n)
    if clip < 0:
        msg = "Pearson residuals require `clip>=0` or `clip=None`."
        raise ValueError(msg)

    if check_values and not check_nonnegative_integers(x):
        msg = "`normalize_pearson_residuals()` expects raw count data, but non-integers were found."
        warn(msg, UserWarning)

    if isinstance(x, CSBase):
        sums_genes = np.sum(x, axis=0)
        sums_cells = np.sum(x, axis=1)
        sum_total = np.sum(sums_genes).squeeze()
    else:
        sums_genes = np.sum(x, axis=0, keepdims=True)
        sums_cells = np.sum(x, axis=1, keepdims=True)
        sum_total = np.sum(sums_genes)

    mu = np.array(sums_cells @ sums_genes / sum_total)
    diff = np.array(x - mu)
    residuals = diff / np.sqrt(mu + mu**2 / theta)

    # clip
    residuals = np.clip(residuals, a_min=-clip, a_max=clip)

    return residuals


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    check_values=doc_check_values,
    layer=doc_layer,
    inplace=doc_inplace,
    copy=doc_copy,
)
def normalize_pearson_residuals(
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    check_values: bool = True,
    layer: str | None = None,
    obsm: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | dict[str, np.ndarray] | None:
    """Apply analytic Pearson residual normalization, based on :cite:t:`Lause2021`.

    The residuals are based on a negative binomial offset model with overdispersion
    `theta` shared across genes. By default, residuals are clipped to `sqrt(n_obs)`
    and overdispersion `theta=100` is used.

    Expects raw count input.

    Params
    ------
    {adata}
    {dist_params}
    {check_values}
    {layer}
    {inplace}
    {copy}

    Returns
    -------
    If `inplace=True`, `adata.X` or the selected layer in `adata.layers` is updated
    with the normalized values. `adata.uns` is updated with the following fields.
    If `inplace=False`, the same fields are returned as dictionary with the
    normalized values in `results_dict['X']`.

    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter.
    `.uns['pearson_residuals_normalization']['computed_on']`
         The name of the layer on which the residuals were computed.

    """
    if copy:
        if not inplace:
            msg = "`copy=True` cannot be used with `inplace=False`."
            raise ValueError(msg)
        adata = adata.copy()

    view_to_actual(adata)
    x = _get_obs_rep(adata, layer=layer, obsm=obsm)
    computed_on = layer or obsm or "adata.X"

    msg = f"computing analytic Pearson residuals on {computed_on}"
    start = logg.info(msg)

    residuals = _pearson_residuals(x, theta, clip, check_values, copy=not inplace)
    settings_dict = dict(theta=theta, clip=clip, computed_on=computed_on)

    if inplace:
        _set_obs_rep(adata, residuals, layer=layer, obsm=obsm)
        adata.uns["pearson_residuals_normalization"] = settings_dict
    else:
        results_dict = dict(X=residuals, **settings_dict)

    logg.info("    finished ({time_passed})", time=start)

    if copy:
        return adata
    elif not inplace:
        return results_dict


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    pca_chunk=doc_pca_chunk,
    mask_var=doc_mask_var,
    check_values=doc_check_values,
    inplace=doc_inplace,
)
@_accepts_legacy_random_state(0)
def normalize_pearson_residuals_pca(
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_comps: int | None = 50,
    rng: SeedLike | RNGLike | None = None,
    kwargs_pca: Mapping[str, Any] = MappingProxyType({}),
    mask_var: np.ndarray | str | None | Default = Default(
        "adata.var.get('highly_variable')"
    ),
    check_values: bool = True,
    inplace: bool = True,
) -> AnnData | None:
    """Apply analytic Pearson residual normalization and PCA, based on :cite:t:`Lause2021`.

    The residuals are based on a negative binomial offset model with overdispersion
    `theta` shared across genes. By default, residuals are clipped to `sqrt(n_obs)`,
    overdispersion `theta=100` is used, and PCA is run with 50 components.

    Operates on the subset of highly variable genes in `adata.var['highly_variable']`
    by default. Expects raw count input.

    Params
    ------
    {adata}
    {dist_params}
    {pca_chunk}
    {mask_var}
    {check_values}
    {inplace}

    Returns
    -------
    If `inplace=False`, returns the Pearson residual-based PCA results (as :class:`~anndata.AnnData`
    object). If `inplace=True`, updates `adata` with the following fields:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
        The subset of highly variable genes, normalized by Pearson residuals.
    `.uns['pearson_residuals_normalization']['theta']`
        The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
        The used value of the clipping parameter.

    `.obsm[kwargs_pca.get('key_added', 'X_pca')]`
        PCA representation of data after gene selection (if applicable) and Pearson
        residual normalization.
    `.varm[kwargs_pca.get('key_added', 'PCs')]`
        The principal components containing the loadings. When `inplace=True` and
        `mask_var is not None`, this will contain empty rows for the genes not
        selected.
    `.uns[kwargs_pca.get('key_added', 'pca')]['variance_ratio']`
        Ratio of explained variance.
    `.uns[kwargs_pca.get('key_added', 'pca')]['variance']`
        Explained variance, equivalent to the eigenvalues of the covariance matrix.

    """
    key_added = kwargs_pca.get("key_added", settings.preset.pca.key_added)
    key_obsm, key_varm, key_uns = (
        ("X_pca", "PCs", "pca") if key_added is None else [key_added] * 3
    )
    if isinstance(mask_var, Default):
        mask_var = "highly_variable" if "highly_variable" in adata.var else None
    mask_var = _check_mask(adata, mask_var, "var")

    if mask_var is not None:
        adata_sub = adata[:, mask_var].copy()
        adata_pca = AnnData(
            adata_sub.X.copy(), obs=adata_sub.obs[[]], var=adata_sub.var[[]]
        )
    else:
        adata_pca = AnnData(adata.X.copy(), obs=adata.obs[[]], var=adata.var[[]])

    normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, rng=rng, **kwargs_pca)
    n_comps = adata_pca.obsm[key_obsm].shape[1]  # might be None

    if inplace:
        norm_settings = adata_pca.uns["pearson_residuals_normalization"]
        norm_dict = dict(**norm_settings, pearson_residuals_df=adata_pca.to_df())
        if mask_var is not None:
            adata.varm[key_varm] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm[key_varm][mask_var] = adata_pca.varm[key_varm]
        else:
            adata.varm[key_varm] = adata_pca.varm[key_varm]
        adata.uns[key_uns] = adata_pca.uns[key_uns]
        adata.uns["pearson_residuals_normalization"] = norm_dict
        adata.obsm[key_obsm] = adata_pca.obsm[key_obsm]
        return None
    else:
        return adata_pca

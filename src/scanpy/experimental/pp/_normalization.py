from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 15):
    from types import MappingProxyType as frozendict  # noqa: N813

import numpy as np
from anndata import AnnData
from scverse_misc import Deprecation, deprecated_arg

from ... import logging as logg
from ..._compat import CSBase, warn
from ..._docs import doc_out, doc_use
from ..._keys import _embedding_keys
from ..._settings import Default
from ..._utils import _doc_params, check_nonnegative_integers, view_to_actual
from ..._utils.random import _accepts_legacy_random_state
from ...experimental._docs import (
    doc_adata,
    doc_check_values,
    doc_dist_params,
    doc_inplace,
    doc_input_matrix,
    doc_pca_chunk,
)
from ...get._out import (
    _legacy_out,
    _out_dests,
    _out_uns_key,
    _write_dests,
    _write_out,
    deprecated_copy,
)
from ...get._use import _check_mask, _get_arr, _mask_arg, _mask_hvg, _resolve_obs
from ...preprocessing._docs import doc_mask_var
from ...preprocessing._pca import PCA_KINDS, PcaOut, _pca_bases, _pca_dests, pca

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from ..._utils.random import RNGLike, SeedLike
    from ...get._out import Prefix
    from ...get._types import LayerAcc, Mask, RepAcc, _SetRep


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
    use=doc_input_matrix,
    out=doc_out("pearson_residuals"),
)
@deprecated_arg("layer", Deprecation("1.13.0", "Use `use` instead."))
@deprecated_arg("obsm", Deprecation("1.13.0", "Use `use` instead."))
@deprecated_arg(
    "inplace", Deprecation("1.13.0", "Use `out=True` instead of `inplace=False`.")
)
@deprecated_copy
def normalize_pearson_residuals(
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    check_values: bool = True,
    use: RepAcc | str | None = None,
    out: RepAcc | str | Prefix | bool = False,
    # deprecated
    layer: str | None = None,
    obsm: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | np.ndarray | dict[str, np.ndarray] | None:
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
    {use}
    {out}
    inplace
        `inplace=False` is `out=True`.
    copy
        Whether to return a modified copy instead of modifying `adata` in place.

    Returns
    -------
    If `out=True`, returns the normalized values.
    Otherwise, writes them to `out` (by default where `use` reads from),
    and updates `adata.uns` with the following fields.
    That `.uns` entry is named after the key `out` writes to
    (e.g. `'k'` for `out=A.layers['k']`, `'p_pearson_residuals'` for `out=sc.Prefix('p')`),
    and `'pearson_residuals_normalization'` by default.
    The deprecated `inplace=False` returns these fields as dictionary with the
    normalized values in `results_dict['X']`.

    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter.
    `.uns['pearson_residuals_normalization']['computed_on']`
         The name of the layer on which the residuals were computed.

    """
    out = _legacy_out(out, "inplace", None if inplace else True)

    view_to_actual(adata)
    use = _resolve_obs(use)
    x = _get_arr(adata, use, layer=layer, obsm=obsm)
    computed_on = layer or obsm or (repr(use) if use is not None else "adata.X")

    msg = f"computing analytic Pearson residuals on {computed_on}"
    start = logg.info(msg)

    residuals = _pearson_residuals(x, theta, clip, check_values, copy=out is not False)
    settings_dict = dict(theta=theta, clip=clip, computed_on=computed_on)
    logg.info("    finished ({time_passed})", time=start)

    if out is True:
        return residuals if inplace else dict(X=residuals, **settings_dict)
    _write_out(
        adata,
        residuals,
        out=out,
        base="pearson_residuals",
        use=use,
        layer=layer,
        obsm=obsm,
    )
    uns_key = _out_uns_key(
        out, base="pearson_residuals", default="pearson_residuals_normalization"
    )
    adata.uns[uns_key] = settings_dict
    return None


def _pca_kwargs_dests(
    kwargs_pca: Mapping[str, Any],
) -> tuple[str | Default | None, PcaOut[_SetRep | None, _SetRep | None]]:
    """Where `pca(**kwargs_pca)` writes, see `_write_pca_from`."""
    if kwargs_pca.get("out") is True:
        msg = "`kwargs_pca['out']=True` is not supported."
        raise ValueError(msg)
    return _out_dests(
        kwargs_pca.get("out", False),
        kwargs_pca.get("key_added"),
        dests=_pca_dests,
        kinds=PCA_KINDS,
        bases=_pca_bases(),
        uns="pca",
        default=Default(),
    )


def _write_pca_from(
    adata: AnnData,
    adata_pca: AnnData,
    dests: PcaOut[_SetRep | None, _SetRep | None],
    mask_var: np.ndarray | None,
) -> None:
    """Copy PCA outputs from `adata_pca` (`adata` subset to `mask_var`) to `adata`."""
    embedding, loadings = (
        None if dest is None else _get_arr(adata_pca, **dest) for dest in dests
    )
    if loadings is not None and mask_var is not None:
        # zero rows for the variables not selected
        padded = np.zeros(shape=(adata.n_vars, loadings.shape[1]))
        padded[mask_var] = loadings
        loadings = padded
    _write_dests(adata, PcaOut(embedding, loadings), dests)


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    pca_chunk=doc_pca_chunk,
    mask=doc_mask_var,
    check_values=doc_check_values,
    use=doc_use("Which count matrix to use as input.", legacy=()),
    inplace=doc_inplace,
)
@_accepts_legacy_random_state(0)
@deprecated_arg("mask_var", Deprecation("1.13.0", "Use `mask` instead."))
@deprecated_arg(
    "inplace",
    Deprecation("1.13.0", "Copy `adata` before calling this function instead."),
)
def normalize_pearson_residuals_pca(  # noqa: PLR0913
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_comps: int | None = 50,
    rng: SeedLike | RNGLike | None = None,
    kwargs_pca: Mapping[str, Any] = frozendict({}),
    mask: Mask | Default | None = Default("adata.var.get('highly_variable')"),
    check_values: bool = True,
    use: LayerAcc | str | None = None,
    # deprecated
    inplace: bool = True,
    mask_var: Mask | None = None,
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
    {mask}
    {check_values}
    {use}
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

    By default, the PCA results are stored in these fields
    (see `out` of :func:`~scanpy.pp.pca`, which can be passed via `kwargs_pca`):

    `.obsm['X_pca']`
        PCA representation of data after gene selection (if applicable) and Pearson
        residual normalization.
    `.varm['PCs']`
        The principal components containing the loadings. When `inplace=True` and
        `mask_var is not None`, this will contain empty rows for the genes not
        selected.
    `.uns['pca']['variance_ratio']`
        Ratio of explained variance.
    `.uns['pca']['variance']`
        Explained variance, equivalent to the eigenvalues of the covariance matrix.

    """
    key, dests = _pca_kwargs_dests(kwargs_pca)
    keys = _embedding_keys("pca", key)
    mask = _mask_hvg(adata, _mask_arg(mask, mask_var, dim="var"))
    mask_var = _check_mask(adata, mask, "var")

    x, var = _get_arr(adata, _resolve_obs(use)), adata.var[[]]
    if mask_var is not None:
        x, var = x[:, mask_var], var[mask_var]
    adata_pca = AnnData(x.copy(), obs=adata.obs[[]], var=var)

    normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, rng=rng, **kwargs_pca)

    if inplace:
        norm_settings = adata_pca.uns["pearson_residuals_normalization"]
        norm_dict = dict(**norm_settings, pearson_residuals_df=adata_pca.to_df())
        _write_pca_from(adata, adata_pca, dests, mask_var)
        adata.uns[keys.uns] = adata_pca.uns[keys.uns]
        adata.uns["pearson_residuals_normalization"] = norm_dict
        return None
    else:
        return adata_pca

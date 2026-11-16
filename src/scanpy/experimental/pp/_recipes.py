from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 15):
    from types import MappingProxyType as frozendict  # noqa: N813


from scverse_misc import Deprecation, deprecated_arg

from ... import experimental
from ..._docs import DEPR_BATCH_KEY, doc_out_key, doc_use
from ..._keys import _embedding_keys
from ..._utils import _doc_params, own_deprecations
from ..._utils.random import _accepts_legacy_random_state
from ...experimental._docs import (
    doc_adata,
    doc_batch_key,
    doc_check_values,
    doc_dist_params,
    doc_genes_batch_chunk,
    doc_pca_chunk,
)
from ...get._out import _legacy_out
from ...get._use import _get_arr, _resolve_layer
from ...preprocessing import pca
from ._normalization import _pca_kwargs_dests, _write_pca_from

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    import pandas as pd
    from anndata import AnnData

    from ..._utils.random import RNGLike, SeedLike
    from ...get._types import AdRef, LayerAcc


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    genes_batch_chunk=doc_genes_batch_chunk,
    pca_chunk=doc_pca_chunk,
    check_values=doc_check_values,
    use=doc_use("Which count matrix to use as input.", legacy=(), layers_only=True),
    out=doc_out_key(
        "the gene selection and Pearson residual-based PCA results",
        default="`adata` (see below)",
    ),
    batch_key=doc_batch_key,
)
@_accepts_legacy_random_state(0)
@deprecated_arg(
    "inplace", Deprecation("1.13.0", "Use `out=True` instead of `inplace=False`.")
)
@deprecated_arg("batch_key", DEPR_BATCH_KEY)
def recipe_pearson_residuals(  # noqa: PLR0913
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_top_genes: int = 1000,
    batch: AdRef | str | None = None,
    chunksize: int = 1000,
    n_comps: int | None = 50,
    rng: SeedLike | RNGLike | None = None,
    kwargs_pca: Mapping[str, Any] = frozendict({}),
    check_values: bool = True,
    use: LayerAcc | str | None = None,
    out: bool = False,
    # deprecated
    inplace: bool = True,
    batch_key: str | None = None,
) -> tuple[AnnData, pd.DataFrame] | None:
    """Full pipeline for HVG selection and normalization by analytic Pearson residuals :cite:p:`Lause2021`.

    Applies gene selection based on Pearson residuals. On the resulting subset,
    Pearson residual normalization and PCA are performed.

    Expects raw count input.

    Params
    ------
    {adata}
    {dist_params}
    {genes_batch_chunk}
    {pca_chunk}
    {check_values}
    {use}
    {out}
    inplace
        `inplace=False` is `out=True`.
    {batch_key}

    Returns
    -------
    If `out=True`, separately returns the Pearson residual-based PCA results (as
    :class:`~anndata.AnnData`) and the gene selection results (as
    :class:`~pandas.DataFrame`). Otherwise, updates `adata` with the
    following fields for gene selection results:

    `.var['highly_variable']` : bool
        boolean indicator of highly-variable genes.
    `.var['means']` : float
        means per gene.
    `.var['variances']` : float
        variances per gene.
    `.var['residual_variances']` : float
        Pearson residual variance per gene. Averaged in the case of multiple
        batches.
    `.var['highly_variable_rank']` : float
        Rank of the gene according to residual variance, median rank in the
        case of multiple batches.
    `.var['highly_variable_nbatches']` : int
        If `batch` is given, this denotes in how many batches genes are
        detected as HVG.
    `.var['highly_variable_intersection']` : bool
        If `batch` is given, this denotes the genes that are highly variable
        in all batches.

    The following fields contain Pearson residual-based PCA results and
    normalization settings:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
         The subset of highly variable genes, normalized by Pearson residuals.
    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter.

    By default, the PCA results are stored in these fields
    (see `out` of :func:`~scanpy.pp.pca`, which can be passed via `kwargs_pca`):

    `.obsm['X_pca']`
        PCA representation of data after gene selection and Pearson residual
        normalization.
    `.varm['PCs']`
         The principal components containing the loadings.
         This will contain empty rows for the genes not selected during HVG selection.
    `.uns['pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pca']['variance']`
         Explained variance, equivalent to the eigenvalues of the covariance matrix.

    """
    out = _legacy_out(out, "inplace", None if inplace else True)
    key, dests = _pca_kwargs_dests(kwargs_pca)
    keys = _embedding_keys("pca", key)
    hvg_args = dict(
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        # forward the deprecated `batch_key` as-is: its strings are column names
        **(dict(batch=batch) if batch_key is None else dict(batch_key=batch_key)),
        theta=theta,
        clip=clip,
        chunksize=chunksize,
        check_values=check_values,
        use=use,
    )

    if out is not True:
        with own_deprecations():
            experimental.pp.highly_variable_genes(adata, **hvg_args)
        # TODO: are these copies needed?
        adata_pca = adata[:, adata.var["highly_variable"]].copy()
    else:
        with own_deprecations():
            hvg = experimental.pp.highly_variable_genes(adata, **hvg_args, out=True)
        # TODO: are these copies needed?
        adata_pca = adata[:, hvg["highly_variable"]].copy()

    if (use := _resolve_layer(use)) is not None:
        adata_pca.X = _get_arr(adata_pca, use)
    experimental.pp.normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, rng=rng, **kwargs_pca)

    if out is True:
        return adata_pca, hvg
    normalization_param = adata_pca.uns["pearson_residuals_normalization"]
    normalization_dict = dict(
        **normalization_param, pearson_residuals_df=adata_pca.to_df()
    )
    adata.uns[keys.uns] = adata_pca.uns[keys.uns]
    _write_pca_from(adata, adata_pca, dests, adata.var["highly_variable"].to_numpy())
    adata.uns["pearson_residuals_normalization"] = normalization_dict
    return None

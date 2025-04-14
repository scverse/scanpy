from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np

from scanpy import experimental
from scanpy._utils import _doc_params
from scanpy.experimental._docs import (
    doc_adata,
    doc_check_values,
    doc_dist_params,
    doc_genes_batch_chunk,
    doc_inplace,
    doc_pca_chunk,
)
from scanpy.preprocessing import pca

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    import pandas as pd
    from anndata import AnnData


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    genes_batch_chunk=doc_genes_batch_chunk,
    pca_chunk=doc_pca_chunk,
    check_values=doc_check_values,
    inplace=doc_inplace,
)
def recipe_pearson_residuals(  # noqa: PLR0913
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_top_genes: int = 1000,
    batch_key: str | None = None,
    chunksize: int = 1000,
    n_comps: int | None = 50,
    random_state: float | None = 0,
    kwargs_pca: Mapping[str, Any] = MappingProxyType({}),
    check_values: bool = True,
    inplace: bool = True,
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
    {inplace}

    Returns
    -------
    If `inplace=False`, separately returns the gene selection results (as
    :class:`~pandas.DataFrame`) and Pearson residual-based PCA results (as
    :class:`~anndata.AnnData`). If `inplace=True`, updates `adata` with the
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
        If batch_key is given, this denotes in how many batches genes are
        detected as HVG.
    `.var['highly_variable_intersection']` : bool
        If batch_key is given, this denotes the genes that are highly variable
        in all batches.

    The following fields contain Pearson residual-based PCA results and
    normalization settings:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
         The subset of highly variable genes, normalized by Pearson residuals.
    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter.

    `.obsm['X_pca']`
        PCA representation of data after gene selection and Pearson residual
        normalization.
    `.varm['PCs']`
         The principal components containing the loadings. When `inplace=True` this
         will contain empty rows for the genes not selected during HVG selection.
    `.uns['pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pca']['variance']`
         Explained variance, equivalent to the eigenvalues of the covariance matrix.

    """
    hvg_args = dict(
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        theta=theta,
        clip=clip,
        chunksize=chunksize,
        check_values=check_values,
    )

    if inplace:
        experimental.pp.highly_variable_genes(adata, **hvg_args, inplace=True)
        # TODO: are these copies needed?
        adata_pca = adata[:, adata.var["highly_variable"]].copy()
    else:
        hvg = experimental.pp.highly_variable_genes(adata, **hvg_args, inplace=False)
        # TODO: are these copies needed?
        adata_pca = adata[:, hvg["highly_variable"]].copy()

    experimental.pp.normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, random_state=random_state, **kwargs_pca)

    if inplace:
        normalization_param = adata_pca.uns["pearson_residuals_normalization"]
        normalization_dict = dict(
            **normalization_param, pearson_residuals_df=adata_pca.to_df()
        )

        adata.uns["pca"] = adata_pca.uns["pca"]
        adata.varm["PCs"] = np.zeros(shape=(adata.n_vars, n_comps))
        adata.varm["PCs"][adata.var["highly_variable"]] = adata_pca.varm["PCs"]
        adata.uns["pearson_residuals_normalization"] = normalization_dict
        adata.obsm["X_pca"] = adata_pca.obsm["X_pca"]
        return None
    else:
        return adata_pca, hvg

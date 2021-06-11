"""Preprocessing recipes from the literature"""
from typing import Optional, Union, Literal, Tuple

from anndata import AnnData

import pandas as pd

from .. import preprocessing as pp
from ._deprecated.highly_variable_genes import (
    filter_genes_dispersion,
    filter_genes_cv_deprecated,
)
from ._normalization import normalize_total
from .. import logging as logg
from .._utils import AnyRandom


def recipe_weinreb17(
    adata: AnnData,
    log: bool = True,
    mean_threshold: float = 0.01,
    cv_threshold: int = 2,
    n_pcs: int = 50,
    svd_solver='randomized',
    random_state: AnyRandom = 0,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Normalization and filtering as of [Weinreb17]_.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    log
        Logarithmize data?
    copy
        Return a copy if true.
    """
    from ._deprecated import normalize_per_cell_weinreb16_deprecated, zscore_deprecated
    from scipy.sparse import issparse

    if issparse(adata.X):
        raise ValueError('`recipe_weinreb16 does not support sparse matrices.')
    if copy:
        adata = adata.copy()
    if log:
        pp.log1p(adata)
    adata.X = normalize_per_cell_weinreb16_deprecated(
        adata.X, max_fraction=0.05, mult_with_mean=True
    )
    gene_subset = filter_genes_cv_deprecated(adata.X, mean_threshold, cv_threshold)
    adata._inplace_subset_var(gene_subset)  # this modifies the object itself
    X_pca = pp.pca(
        zscore_deprecated(adata.X),
        n_comps=n_pcs,
        svd_solver=svd_solver,
        random_state=random_state,
    )
    # update adata
    adata.obsm['X_pca'] = X_pca
    return adata if copy else None


def recipe_seurat(
    adata: AnnData, log: bool = True, plot: bool = False, copy: bool = False
) -> Optional[AnnData]:
    """\
    Normalization and filtering as of Seurat [Satija15]_.

    This uses a particular preprocessing.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.
    """
    if copy:
        adata = adata.copy()
    pp.filter_cells(adata, min_genes=200)
    pp.filter_genes(adata, min_cells=3)
    normalize_total(adata, target_sum=1e4)
    filter_result = filter_genes_dispersion(
        adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5, log=not log
    )
    if plot:
        from ..plotting import (
            _preprocessing as ppp,
        )  # should not import at the top of the file

        ppp.filter_genes_dispersion(filter_result, log=not log)
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    if log:
        pp.log1p(adata)
    pp.scale(adata, max_value=10)
    return adata if copy else None


def recipe_zheng17(
    adata: AnnData,
    n_top_genes: int = 1000,
    log: bool = True,
    plot: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Normalization and filtering as of [Zheng17]_.

    Reproduces the preprocessing of [Zheng17]_ – the Cell Ranger R Kit of 10x
    Genomics.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.

    The recipe runs the following steps

    .. code:: python

        sc.pp.filter_genes(adata, min_counts=1)         # only consider genes with more than 1 count
        sc.pp.normalize_per_cell(                       # normalize with total UMI count per cell
             adata, key_n_counts='n_counts_all'
        )
        filter_result = sc.pp.filter_genes_dispersion(  # select highly-variable genes
            adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False
        )
        adata = adata[:, filter_result.gene_subset]     # subset the genes
        sc.pp.normalize_per_cell(adata)                 # renormalize after filtering
        if log: sc.pp.log1p(adata)                      # log transform: adata.X = log(adata.X + 1)
        sc.pp.scale(adata)                              # scale to unit variance and shift to zero mean


    Parameters
    ----------
    adata
        Annotated data matrix.
    n_top_genes
        Number of genes to keep.
    log
        Take logarithm.
    plot
        Show a plot of the gene dispersion vs. mean relation.
    copy
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    start = logg.info('running recipe zheng17')
    if copy:
        adata = adata.copy()
    # only consider genes with more than 1 count
    pp.filter_genes(adata, min_counts=1)
    # normalize with total UMI count per cell
    normalize_total(adata, key_added='n_counts_all')
    filter_result = filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False
    )
    if plot:  # should not import at the top of the file
        from ..plotting import _preprocessing as ppp

        ppp.filter_genes_dispersion(filter_result, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter_result.gene_subset]
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    normalize_total(adata)  # renormalize after filtering
    if log:
        pp.log1p(adata)  # log transform: X = log(X + 1)
    pp.scale(adata)
    logg.info('    finished', time=start)
    return adata if copy else None


def recipe_pearson_residuals(
    adata: AnnData,
    n_top_genes: int = 1000,
    theta: float = 100,
    clip: Optional[float] = None,
    chunksize: int = 1000,
    batch_key: Optional[str] = None,
    n_comps_pca: Optional[int] = 50,
    random_state_pca: Optional[float] = 0,
    kwargs_pca: dict = {},
    check_values: bool = True,
    inplace: bool = True,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """\
    Applies gene selection based on Pearson residuals. On the resulting subset,
    Pearson residual normalization and PCA are performed.

    This recipe is based on "Analytic Pearson residuals for normalization of
    single-cell RNA-seq UMI data", bioRxiv, [Lause20]_.


    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    n_top_genes
        Number of highly-variable genes to keep. Mandatory if
        `flavor='seurat_v3'` or `flavor='pearson_residuals'`.
    chunksize
        This dertermines how many genes are processed at once while computing
        the Pearson residual variance. Choosing a smaller value will reduce
        the required memory.
    theta
        This is the NB overdispersion parameter theta for Pearson residual
        computations. Higher values correspond to less overdispersion
        (var = mean + mean^2/theta), and `theta=np.Inf` corresponds to a
        Poisson model.
    clip
        This determines if and how Pearson residuals are clipped:

        * If `None`, residuals are clipped to the interval
        [-sqrt(n), sqrt(n)], where n is the number of cells in the dataset
        (default behavior).
        * If any scalar c, residuals are clipped to the interval [-c, c]. Set
        `clip=np.Inf` for no clipping.
    batch_key
        If specified, highly-variable genes are selected within each batch
        separately and merged. This simple process avoids the selection of
        batch-specific genes and acts as a lightweight batch correction
        method. For all flavors, genes are first sorted by how many batches
        they are a HVG. Ties are broken by the median rank (across batches)
        based on within-batch residual variance.

    n_comps_pca
        Number of principal components to compute.
    random_state_pca
        Change to use different initial states for the optimization.
    kwargs_pca
        Dictionary of further keyword arguments passed on to `sc.pp.pca()`.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to True.
    inplace
        Whether to place results in `adata` or return them.

    Returns
    ------
    If `inplace=False`, separately returns the gene selection results (`hvg`)
    and Pearson residual-based PCA results (`adata_pca`). If `inplace=True`,
    updates `adata` with the following fields for gene selection results…:

    `.var['highly_variable']`
        boolean indicator of highly-variable genes.
    `.var['means']`
        means per gene.
    `.var['variances']`
        variances per gene.
    `.var['residual_variances']`
        Pearson residual variance per gene. Averaged in the case of multiple
        batches.
    `.var['highly_variable_rank']`
        Rank of the gene according to residual variance, median rank in the
        case of multiple batches.
    `.var['highly_variable_nbatches']`
        If batch_key is given, this denotes in how many batches genes are
        detected as HVG.
    `.var['highly_variable_intersection']`
        If batch_key is given, this denotes the genes that are highly variable
        in all batches.

    …and the following fields for Pearson residual-based PCA results and
    normalization settings:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
         The hvg-subset, normalized by Pearson residuals.
    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta.
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter.

    `.obsm['pearson_residuals_X_pca']`
        PCA representation of data after gene selection and Pearson residual
        normalization.
    `.uns['pearson_residuals_pca']['PCs']`
         The principal components containing the loadings.
    `.uns['pearson_residuals_pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pearson_residuals_pca']['variance']`
         Explained variance, equivalent to the eigenvalues of the
         covariance matrix.

    """

    hvg_args = dict(
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        theta=theta,
        clip=clip,
        chunksize=chunksize,
        check_values=check_values,
    )

    if inplace:
        pp.highly_variable_genes(adata, **hvg_args, inplace=True)
        # TODO: are these copies needed?
        adata_pca = adata[:, adata.var['highly_variable']].copy()
    else:
        hvg = pp.highly_variable_genes(adata, **hvg_args, inplace=False)
        # TODO: are these copies needed?
        adata_pca = adata[:, hvg['highly_variable']].copy()

    pp.normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pp.pca(adata_pca, n_comps=n_comps_pca, random_state=random_state_pca, **kwargs_pca)

    if inplace:
        normalization_param = adata_pca.uns['pearson_residuals_normalization']
        normalization_dict = dict(
            **normalization_param, pearson_residuals_df=adata_pca.to_df()
        )
        pca_param = adata_pca.uns['pca']
        pca_dict = dict(**pca_param, PCs=adata_pca.varm['PCs'])
        adata.uns['pearson_residuals_pca'] = pca_dict
        adata.uns['pearson_residuals_normalization'] = normalization_dict
        adata.obsm['X_pearson_residuals_pca'] = adata_pca.obsm['X_pca']
        return None
    else:
        adata_pca.obsm['X_pearson_residuals_pca'] = adata_pca.obsm['X_pca'].copy()
        adata_pca.uns['pearson_residuals_pca'] = adata_pca.uns['pca'].copy()
        del adata_pca.obsm['X_pca']
        del adata_pca.uns['pca']
        return adata_pca, hvg

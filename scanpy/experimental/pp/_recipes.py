from typing import Optional, Tuple
from anndata import AnnData
import pandas as pd
import numpy as np
from scanpy import experimental
from scanpy.preprocessing import pca


def recipe_pearson_residuals(
    adata: AnnData,
    *,
    theta: float = 100,
    clip: Optional[float] = None,
    n_top_genes: int = 1000,
    batch_key: Optional[str] = None,
    chunksize: int = 1000,
    n_comps: Optional[int] = 50,
    random_state: Optional[float] = 0,
    kwargs_pca: dict = {},
    check_values: bool = True,
    inplace: bool = True,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """\
    Gene selection and normalization based on [Lause21]_.

    Applies gene selection based on Pearson residuals. On the resulting subset,
    Pearson residual normalization and PCA are performed.

    Expects raw count input.


    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    theta
        The negative binomial overdispersion parameter theta for Pearson residuals.
        Higher values correspond to less overdispersion (var = mean + mean^2/theta),
        and `theta=np.Inf` corresponds to a Poisson model.
    clip
        Determines if and how residuals are clipped:

            * If `None`, residuals are clipped to the interval [-sqrt(n), sqrt(n)], \
            where n is the number of cells in the dataset (default behavior).
            * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
            `clip=np.Inf` for no clipping.

    n_top_genes
        Number of highly-variable genes to keep.
    batch_key
        If specified, highly-variable genes are selected within each batch separately
        and merged. This simple process avoids the selection of batch-specific genes
        and acts as a lightweight batch correction method. Genes are first sorted by
        how many batches they are a HVG. Ties are broken by the median rank (across
        batches) based on within-batch residual variance.
    chunksize
        This dertermines how many genes are processed at once while computing
        the Pearson residual variance. Choosing a smaller value will reduce
        the required memory.
    n_comps
        Number of principal components to compute in the PCA step.
    random_state
        Change to use different initial states for the optimization in the PCA step.
    kwargs_pca
        Dictionary of further keyword arguments passed on to `scanpy.pp.pca()`.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to
        True.
    inplace
        Whether to place results in `adata` or return them.


    Returns
    -------
    If `inplace=False`, separately returns the gene selection results (`hvg`)
    and Pearson residual-based PCA results (`adata_pca`). If `inplace=True`,
    updates `adata` with the following fields for gene selection results:

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

    …and the following fields for Pearson residual-based PCA results and
    normalization settings:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
         The hvg-subset, normalized by Pearson residuals.
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
        flavor='pearson_residuals',
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
        adata_pca = adata[:, adata.var['highly_variable']].copy()
    else:
        hvg = experimental.pp.highly_variable_genes(adata, **hvg_args, inplace=False)
        # TODO: are these copies needed?
        adata_pca = adata[:, hvg['highly_variable']].copy()

    experimental.pp.normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, random_state=random_state, **kwargs_pca)

    if inplace:
        normalization_param = adata_pca.uns['pearson_residuals_normalization']
        normalization_dict = dict(
            **normalization_param, pearson_residuals_df=adata_pca.to_df()
        )

        adata.uns['pca'] = adata_pca.uns['pca']
        adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
        adata.varm['PCs'][adata.var['highly_variable']] = adata_pca.varm['PCs']
        adata.uns['pearson_residuals_normalization'] = normalization_dict
        adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
        return None
    else:
        return adata_pca, hvg

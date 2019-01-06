"""Preprocessing recipes from the literature
"""

from . import _simple as pp
from ._deprecated.highly_variable_genes import filter_genes_dispersion, filter_genes_cv_deprecated
from .. import logging as logg


def recipe_weinreb17(adata, log=True, mean_threshold=0.01, cv_threshold=2,
                     n_pcs=50, svd_solver='randomized', random_state=0,
                     copy=False):
    """Normalization and filtering as of [Weinreb17]_.

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    copy : bool (default: False)
        Return a copy if true.
    """
    from scipy.sparse import issparse
    if issparse(adata.X):
        raise ValueError('`recipe_weinreb16 does not support sparse matrices.')
    if copy: adata = adata.copy()
    if log: pp.log1p(adata)
    adata.X = pp.normalize_per_cell_weinreb16_deprecated(adata.X,
                                                         max_fraction=0.05,
                                                         mult_with_mean=True)
    gene_subset = filter_genes_cv_deprecated(adata.X, mean_threshold, cv_threshold)
    adata._inplace_subset_var(gene_subset)  # this modifies the object itself
    X_pca = pp.pca(pp.zscore_deprecated(adata.X),
                   n_comps=n_pcs, svd_solver=svd_solver, random_state=random_state)
    # update adata
    adata.obsm['X_pca'] = X_pca
    return adata if copy else None


def recipe_seurat(adata, log=True, plot=False, copy=False):
    """Normalization and filtering as of Seurat [Satija15]_.

    This uses a particular preprocessing.

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.
    """
    if copy: adata = adata.copy()
    pp.filter_cells(adata, min_genes=200)
    pp.filter_genes(adata, min_cells=3)
    pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    filter_result = filter_genes_dispersion(
        adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5, log=not log)
    if plot:
        from ..plotting import _preprocessing as ppp  # should not import at the top of the file
        ppp.filter_genes_dispersion(filter_result, log=not log)
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    if log: pp.log1p(adata)
    pp.scale(adata, max_value=10)
    return adata if copy else None


def recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False):
    """Normalization and filtering as of [Zheng17]_.

    Reproduces the preprocessing of [Zheng17]_ - the Cell Ranger R Kit of 10x
    Genomics.

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    The recipe runs the following steps

    .. code:: python

        sc.pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
        sc.pp.normalize_per_cell(                # normalize with total UMI count per cell
             adata, key_n_counts='n_counts_all')
        filter_result = sc.pp.filter_genes_dispersion(  # select highly-variable genes
            adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
        adata = adata[:, filter_result.gene_subset]     # subset the genes
        sc.pp.normalize_per_cell(adata)          # renormalize after filtering
        if log: sc.pp.log1p(adata)               # log transform: adata.X = log(adata.X + 1)
        sc.pp.scale(adata)                       # scale to unit variance and shift to zero mean


    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_top_genes : `int`, optional (default: 1000)
        Number of genes to keep.
    log : `bool`, optional (default: `True`)
        Take logarithm.
    plot : `bool`, optional (default: `True`)
        Show a plot of the gene dispersion vs. mean relation.
    copy : `bool`, optional (default: `False`)
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    logg.info('running recipe zheng17', reset=True)
    if copy: adata = adata.copy()
    pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
    pp.normalize_per_cell(adata,  # normalize with total UMI count per cell
                          key_n_counts='n_counts_all')
    filter_result = filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
    if plot:
        from ..plotting import _preprocessing as ppp  # should not import at the top of the file
        ppp.filter_genes_dispersion(filter_result, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter_result.gene_subset]
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    pp.normalize_per_cell(adata)  # renormalize after filtering
    if log: pp.log1p(adata)  # log transform: X = log(X + 1)
    pp.scale(adata)
    logg.info('    finished', time=True)
    return adata if copy else None

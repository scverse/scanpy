"""Preprocessing recipes from the literature."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ... import logging as logg
from ... import preprocessing as pp
from ..._compat import CSBase
from ..._utils.random import _accepts_legacy_random_state
from .._normalization import normalize_total

if TYPE_CHECKING:
    from anndata import AnnData

    from ..._utils.random import RNGLike, SeedLike


@_accepts_legacy_random_state(0)
def recipe_weinreb17(
    adata: AnnData,
    *,
    log: bool = True,
    mean_threshold: float = 0.01,
    cv_threshold: int = 2,
    n_pcs: int = 50,
    svd_solver="randomized",
    rng: SeedLike | RNGLike | None = None,
    copy: bool = False,
) -> AnnData | None:
    """Normalize and filter as of :cite:p:`Weinreb2017`.

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
    from .weinreb import filter_genes_cv, normalize_per_cell, zscore

    if isinstance(adata.X, CSBase):
        msg = "`recipe_weinreb16 does not support sparse matrices."
        raise ValueError(msg)
    if copy:
        adata = adata.copy()
    if log:
        pp.log1p(adata)
    adata.X = normalize_per_cell(adata.X, max_fraction=0.05, mult_with_mean=True)
    gene_subset = filter_genes_cv(adata.X, mean_threshold, cv_threshold)
    adata._inplace_subset_var(gene_subset)  # this modifies the object itself
    x_pca = pp.pca(
        zscore(adata.X),
        n_comps=n_pcs,
        svd_solver=svd_solver,
        rng=rng,
    )
    # update adata
    adata.obsm["X_pca"] = x_pca
    return adata if copy else None


def recipe_seurat(
    adata: AnnData, *, log: bool = True, plot: bool = False, copy: bool = False
) -> AnnData | None:
    """Normalize and filter as of Seurat :cite:p:`Satija2015`.

    This uses a particular preprocessing.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    log
        Logarithmize data?
    plot
        Show a plot of the gene dispersion vs. mean relation.
    copy
        Return a copy if true.

    """
    if copy:
        adata = adata.copy()
    pp.filter_cells(adata, min_genes=200)
    pp.filter_genes(adata, min_cells=3)
    normalize_total(adata, target_sum=1e4)
    adata.layers[layer_log := "log1p"] = adata.X
    pp.log1p(adata, layer=layer_log)
    filter_result = pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer=layer_log, inplace=False
    )
    assert filter_result is not None
    if plot:
        from ...plotting import _preprocessing as ppp

        ppp.highly_variable_genes(filter_result, log=not log)
    adata._inplace_subset_var(filter_result["highly_variable"])  # filter genes
    if log:
        adata.X = adata.layers[layer_log]
    del adata.layers[layer_log]
    pp.scale(adata, max_value=10)
    return adata if copy else None


def recipe_zheng17(
    adata: AnnData,
    *,
    n_top_genes: int = 1000,
    log: bool = True,
    plot: bool = False,
    copy: bool = False,
) -> AnnData | None:
    """Normalize and filter as of :cite:t:`Zheng2017`.

    Reproduces the preprocessing of :cite:t:`Zheng2017` – the Cell Ranger R Kit of 10x
    Genomics.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.


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
    start = logg.info("running recipe zheng17")
    if copy:
        adata = adata.copy()
    # only consider genes with more than 1 count
    pp.filter_genes(adata, min_counts=1)
    # normalize with total UMI count per cell
    normalize_total(adata, key_added="n_counts_all")
    adata.layers[layer_log := "log1p"] = adata.X
    pp.log1p(adata, layer=layer_log)
    filter_result = pp.highly_variable_genes(
        adata,
        flavor="cell_ranger",
        n_top_genes=n_top_genes,
        layer=layer_log,
        inplace=False,
    )
    assert filter_result is not None
    if plot:  # should not import at the top of the file
        from ...plotting import _preprocessing as ppp

        ppp.highly_variable_genes(filter_result, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter_result["gene_subset"]]
    adata._inplace_subset_var(filter_result["highly_variable"])  # filter genes
    normalize_total(adata)  # renormalize after filtering
    if log:
        adata.X = adata.layers[layer_log]
    del adata.layers[layer_log]
    pp.scale(adata)
    logg.info("    finished", time=start)
    return adata if copy else None

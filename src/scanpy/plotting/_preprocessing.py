from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib import rcParams

from .._settings import settings
from ._utils import savefig_or_show

# --------------------------------------------------------------------------------
# Plot result of preprocessing functions
# --------------------------------------------------------------------------------


def highly_variable_genes(  # noqa: PLR0912
    adata_or_result: AnnData | pd.DataFrame | np.recarray,
    *,
    log: bool = False,
    show: bool | None = None,
    highly_variable_genes: bool = True,
    # deprecated
    save: bool | str | None = None,
) -> None:
    """Plot dispersions or normalized variance versus means for genes.

    Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() and
    VariableFeaturePlot() of Seurat.

    Parameters
    ----------
    adata
        Result of :func:`~scanpy.pp.highly_variable_genes`.
    log
        Plot on logarithmic axes.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.

    Examples
    --------
    Compute and plot highly variable genes from raw PBMC data.

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata)

    Plot on logarithmic axes.

    .. plot::
        :context: close-figs

        sc.pl.highly_variable_genes(adata, log=True)

    """
    if isinstance(adata_or_result, AnnData):
        result = adata_or_result.var
        seurat_v3_flavor = adata_or_result.uns["hvg"]["flavor"] == "seurat_v3"
    else:
        result = adata_or_result
        if isinstance(result, pd.DataFrame):
            seurat_v3_flavor = "variances_norm" in result.columns
        else:
            seurat_v3_flavor = False
    if highly_variable_genes:
        gene_subset = result["highly_variable"]
    else:
        gene_subset = result["gene_subset"]
    means = result["means"]

    if seurat_v3_flavor:
        var_or_disp = result["variances"]
        var_or_disp_norm = result["variances_norm"]
    else:
        var_or_disp = result["dispersions"]
        var_or_disp_norm = result["dispersions_norm"]
    size = rcParams["figure.figsize"]
    plt.figure(figsize=(2 * size[0], size[1]))
    plt.subplots_adjust(wspace=0.3)
    for idx, d in enumerate([var_or_disp_norm, var_or_disp]):
        plt.subplot(1, 2, idx + 1)
        for label, color, mask in zip(
            ["highly variable genes", "other genes"],
            ["black", "grey"],
            [gene_subset, ~gene_subset],
            strict=True,
        ):
            if False:
                means_, var_or_disps_ = np.log10(means[mask]), np.log10(d[mask])
            else:
                means_, var_or_disps_ = means[mask], d[mask]
            plt.scatter(means_, var_or_disps_, label=label, c=color, s=1)
        if log:  # there's a bug in autoscale
            plt.xscale("log")
            plt.yscale("log")
            y_min = np.min(var_or_disp)
            y_min = 0.95 * y_min if y_min > 0 else 1e-1
            plt.xlim(0.95 * np.min(means), 1.05 * np.max(means))
            plt.ylim(y_min, 1.05 * np.max(var_or_disp))
        if idx == 0:
            plt.legend()
        plt.xlabel(f"{'$log_{10}$ ' if False else ''}mean expressions of genes")
        data_type = "dispersions" if not seurat_v3_flavor else "variances"
        plt.ylabel(
            f"{'$log_{10}$ ' if False else ''}{data_type} "
            f"of genes ({'' if idx == 0 else 'not '}normalized)"
        )

    show = settings.autoshow if show is None else show
    savefig_or_show("highly_variable_genes", show=show, save=save)
    if show:
        return None
    return plt.gca()

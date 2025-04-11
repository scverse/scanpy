from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib import rcParams

from .._compat import deprecated, old_positionals
from .._settings import settings
from ._utils import savefig_or_show

# --------------------------------------------------------------------------------
# Plot result of preprocessing functions
# --------------------------------------------------------------------------------


@old_positionals("log", "show", "save", "highly_variable_genes")
def highly_variable_genes(  # noqa: PLR0912
    adata_or_result: AnnData | pd.DataFrame | np.recarray,
    *,
    log: bool = False,
    show: bool | None = None,
    save: bool | str | None = None,
    highly_variable_genes: bool = True,
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
        gene_subset = result.highly_variable
    else:
        gene_subset = result.gene_subset
    means = result.means

    if seurat_v3_flavor:
        var_or_disp = result.variances
        var_or_disp_norm = result.variances_norm
    else:
        var_or_disp = result.dispersions
        var_or_disp_norm = result.dispersions_norm
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
        plt.xlabel(("$log_{10}$ " if False else "") + "mean expressions of genes")
        data_type = "dispersions" if not seurat_v3_flavor else "variances"
        plt.ylabel(
            ("$log_{10}$ " if False else "")
            + f"{data_type} of genes"
            + (" (normalized)" if idx == 0 else " (not normalized)")
        )

    show = settings.autoshow if show is None else show
    savefig_or_show("filter_genes_dispersion", show=show, save=save)
    if show:
        return None
    return plt.gca()


# backwards compat
@deprecated("Use sc.pl.highly_variable_genes instead")
@old_positionals("log", "show", "save")
def filter_genes_dispersion(
    result: np.recarray,
    *,
    log: bool = False,
    show: bool | None = None,
    save: bool | str | None = None,
) -> None:
    """Plot dispersions versus means for genes.

    Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() of Seurat.

    Parameters
    ----------
    result
        Result of :func:`~scanpy.pp.filter_genes_dispersion`.
    log
        Plot on logarithmic axes.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.

    """
    highly_variable_genes(
        result, log=log, show=show, save=save, highly_variable_genes=False
    )

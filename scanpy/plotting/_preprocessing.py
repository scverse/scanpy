from typing import Optional, Union

import numpy as np
from matplotlib import pyplot as pl
from matplotlib import rcParams
from anndata import AnnData
from . import _utils

# --------------------------------------------------------------------------------
# Plot result of preprocessing functions
# --------------------------------------------------------------------------------


def highly_variable_genes(
    adata_or_result: Union[AnnData, np.recarray],
    log: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    highly_variable_genes: bool = True,
):
    """Plot dispersions versus means for genes.

    Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() of Seurat.

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
    else:
        result = adata_or_result
    if highly_variable_genes:
        gene_subset = result.highly_variable
    else:
        gene_subset = result.gene_subset
    means = result.means
    dispersions = result.dispersions
    dispersions_norm = result.dispersions_norm
    size = rcParams['figure.figsize']
    pl.figure(figsize=(2*size[0], size[1]))
    pl.subplots_adjust(wspace=0.3)
    for idx, d in enumerate([dispersions_norm, dispersions]):
        pl.subplot(1, 2, idx + 1)
        for label, color, mask in zip(['highly variable genes', 'other genes'],
                                      ['black', 'grey'],
                                      [gene_subset, ~gene_subset]):
            if False: means_, disps_ = np.log10(means[mask]), np.log10(d[mask])
            else: means_, disps_ = means[mask], d[mask]
            pl.scatter(means_, disps_, label=label, c=color, s=1)
        if log:  # there's a bug in autoscale
            pl.xscale('log')
            pl.yscale('log')
            min_dispersion = np.min(dispersions)
            y_min = 0.95*min_dispersion if min_dispersion > 0 else 1e-1
            pl.xlim(0.95*np.min(means), 1.05*np.max(means))
            pl.ylim(y_min, 1.05*np.max(dispersions))
        if idx == 0: pl.legend()
        pl.xlabel(('$log_{10}$ ' if False else '') + 'mean expressions of genes')
        pl.ylabel(('$log_{10}$ ' if False else '') + 'dispersions of genes'
                  + (' (normalized)' if idx == 0 else ' (not normalized)'))
    _utils.savefig_or_show('filter_genes_dispersion', show=show, save=save)
    if show is False:
        return pl.gca()


# backwards compat
def filter_genes_dispersion(
    result: np.recarray,
    log: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
):
    """\
    Plot dispersions versus means for genes.

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

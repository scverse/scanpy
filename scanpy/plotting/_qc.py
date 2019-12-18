<<<<<<< HEAD
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from . import _utils as utils
from ..preprocessing._simple import normalize_per_cell
from ..utils import doc_params
from ._docs import doc_show_save_ax


@doc_params(show_save_ax=doc_show_save_ax)
def highest_expr_genes(
        adata, n_top=30, show=None, save=None,
        ax=None, gene_symbols=None, **kwds
    ):
=======
from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.axes import Axes

from . import _utils
from ._docs import doc_show_save_ax
from ..preprocessing._normalization import normalize_total
from .._utils import _doc_params


@_doc_params(show_save_ax=doc_show_save_ax)
def highest_expr_genes(
    adata: AnnData,
    n_top: int = 30,
    show: Optional[bool] = None,
    save: Optional[Union[str, bool]] = None,
    ax: Optional[Axes] = None,
    gene_symbols: Optional[str] = None,
    log: bool = False,
    **kwds,
):
>>>>>>> upstream/master
    """\
    Fraction of counts assigned to each gene over all cells.

    Computes, for each gene, the fraction of counts assigned to that gene within
    a cell. The `n_top` genes with the highest mean fraction over all cells are
    plotted as boxplots.

    This plot is similar to the `scater` package function `plotHighestExprs(type
    = "highest-expression")`, see `here
    <https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-qc.html>`__. Quoting
    from there:

        *We expect to see the “usual suspects”, i.e., mitochondrial genes, actin,
        ribosomal protein, MALAT1. A few spike-in transcripts may also be
        present here, though if all of the spike-ins are in the top 50, it
        suggests that too much spike-in RNA was added. A large number of
        pseudo-genes or predicted genes may indicate problems with alignment.*
        -- Davis McCarthy and Aaron Lun

    Parameters
    ----------
<<<<<<< HEAD
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_top : `int`, optional (default:30)
        Number of top
    {show_save_ax}
    gene_symbols : `str`, optional (default:None)
        Key for field in .var that stores gene symbols if you do not want to use .var_names.
    **kwds : keyword arguments
        Are passed to `seaborn.boxplot`.

    Returns
    -------
    If `show==False` a `matplotlib.Axis`.
    """
    from scipy.sparse import issparse

    # compute the percentage of each gene per cell
    dat = normalize_per_cell(adata, counts_per_cell_after=100, copy=True)

    # identify the genes with the highest mean
    if issparse(dat.X):
        dat.var['mean_percent'] = dat.X.mean(axis=0).A1
    else:
        dat.var['mean_percent'] = dat.X.mean(axis=0)

    top = dat.var.sort_values('mean_percent', ascending=False).index[:n_top]
    dat = dat[:, top]
    columns = dat.var_names if gene_symbols is None else dat.var[gene_symbols]
    dat = pd.DataFrame(dat.X.toarray(), index=dat.obs_names, columns=columns)

    if not ax:
        # figsize is hardcoded to produce a tall image. To change the fig size,
        # a matplotlib.Axes object needs to be passed.
        height = (n_top * 0.2) + 1.5
        fig, ax = plt.subplots(figsize=(5, height))
    sns.boxplot(data=dat, orient='h', ax=ax, fliersize=1, **kwds)
    ax.set_xlabel('% of total counts')
    utils.savefig_or_show('highest_expr_genes', show=show, save=save)
    return ax if show == False else None
=======
    adata
        Annotated data matrix.
    n_top
        Number of top
    {show_save_ax}
    gene_symbols
        Key for field in .var that stores gene symbols if you do not want to use .var_names.
    log
        Plot x-axis in log scale
    **kwds
        Are passed to :func:`~seaborn.boxplot`.

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes`.
    """
    import seaborn as sns  # Slow import, only import if called
    from scipy.sparse import issparse

    # compute the percentage of each gene per cell
    norm_dict = normalize_total(adata, target_sum=100, inplace=False)

    # identify the genes with the highest mean
    if issparse(norm_dict['X']):
        mean_percent = norm_dict['X'].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx].A
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx]
    columns = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols][top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )

    if not ax:
        # figsize is hardcoded to produce a tall image. To change the fig size,
        # a matplotlib.axes.Axes object needs to be passed.
        height = (n_top * 0.2) + 1.5
        fig, ax = plt.subplots(figsize=(5, height))
    sns.boxplot(data=counts_top_genes, orient='h', ax=ax, fliersize=1, **kwds)
    ax.set_xlabel('% of total counts')
    if log:
        ax.set_xscale('log')
    _utils.savefig_or_show('highest_expr_genes', show=show, save=save)
    return ax if not show else None
>>>>>>> upstream/master

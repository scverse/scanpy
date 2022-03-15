from typing import Optional
from anndata import AnnData
from scanpy.plotting import scatter
import numpy as np


def pearson_residuals_hvg_scatter(
    adata: AnnData,
    marker_names=None,
    gene_name_key: Optional[str] = None,
    x: str = 'means',
    y: str = 'residual_variances',
    hvg_key: str = 'highly_variable',
    kwargs_sc_pl_scatter: dict = dict(),
    return_ax: bool = False,
):
    '''\
    A plot inspect the gene selection by Pearson residuals and sanity-check that the
    underlying model is appropriate.

    Expects that `sc.experimental.pp.highly_variable_genes(flavor='pearson_residuals')` has
    been run before with `inplace=True` and `subset=False`. Gene expression mean is plotted
    against gene residual variance, and the HVG selection is highlighted.

    If the underlying model is appropriate for the data at hand, the following patterns can
    be observed:

    * Genes that behave as predicted by the null model (i.e. that are homogeneous across
    cells in all cells) will have residual variance around ~1
    * Genes that are heterogeneous / differentially expressed across cells will have residual
    variance >1
    * little to no relationship between mean and residual variance, indicating
    successful variance stabilization.


    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    marker_names
        An optional list of gene names. Those genes will be highlighted additionally.
    gene_name_key
        If specified, the function will look in `adata.var[gene_name_key]` for gene names.
        Otherwise, `adata.var_names` is used to find the marker genes (default).
    x
        A key to `adata.var` to be plotted on the x-axis. Default is `'means'`, as returned by
        `sc.experimental.pp.highly_variable_genes(flavor='pearson_residuals')`.
    y
        A key to `adata.var` to be plotted on the y-axis. Default is `'residual_variances'`,
        as returned by `sc.experimental.pp.highly_variable_genes(flavor='pearson_residuals')`.
    hvg_key
        A key to `adata.var` that holds the indicator variable for the HVG selection
    kwargs_sc_pl_scatter
        Further keyword arguments passed on to `sc.pl.scatter()`.
    return_ax
        If `True`, an axis object is returned after plotting.

    Returns
    -------
    Depending on `return_ax`, the axis object is returned.
    '''

    def clean_helper_fields(ad):
        del ad.var['hvg_marker_status']
        if 'hvg_marker_status_colors' in ad.uns.keys():
            del ad.uns['hvg_marker_status_colors']

    if gene_name_key is None:
        gene_names = adata.var_names
    else:
        gene_names = adata.var[gene_name_key]

    try:
        adata.var['hvg_marker_status'] = 'others'
        adata.var['hvg_marker_status'].loc[adata.var[hvg_key]] = 'HVGs'
        if marker_names:
            is_marker = np.isin(gene_names, marker_names)
            adata.var['hvg_marker_status'].loc[is_marker] = 'marker genes'
        ax = scatter(
            adata,
            x=x,
            y=y,
            color='hvg_marker_status',
            show=False,
            **kwargs_sc_pl_scatter,
        )

    except Exception as e:
        # clean up helper field also in case of error
        clean_helper_fields(adata)
        raise e
    # clean up helper field
    clean_helper_fields(adata)

    # label axis flexibly depending on chosen x/y field to plot
    if x == 'means':
        ax.set_xlabel('Mean expression')
    else:
        ax.set_xlabel(x)
    if y == 'residual_variances':
        ax.set_ylabel('Residual variance')
    else:
        ax.set_ylabel(y)

    ax.set_title('Pearson residuals:\nGene selection')
    ax.set_yscale('log')
    ax.set_xscale('log')

    if return_ax:
        return ax

"""
Computes a dendrogram based on a given categorical observation.
"""

from typing import Optional, List
import pandas as pd
from anndata import AnnData
from pandas.api.types import is_categorical_dtype

from .. utils import doc_params
from ..tools._utils import choose_representation, doc_use_rep, doc_n_pcs


@doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def dendrogram(adata: AnnData, groupby: str,
               n_pcs: Optional[int]=None,
               use_rep: Optional[str]=None,
               var_names: Optional[List[str]]=None,
               use_raw: Optional[bool]=None,
               cor_method: Optional[str]='pearson',
               linkage_method: Optional[str]='complete',
               key_added: Optional[str]='dendrogram') -> None:

    """
    Computes a hierarchical clustering for the given `groupby` categories. Be default the PCA
    components are used unless .X has less than 50 variables.

    Alternatively, a list of var_names (e.g genes) can be given.

    Average values of either var_names or components are then used to compute a correlation matrix.

    The hierarchical clustering can be visualized using `sc.pl.dendrogram` or multiple other
    visualizations that can include a dendrogram: `matrixplot`, `heatmap`, `dotplot` and `stacked_violin`

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix
    {n_pcs}
    {use_rep}
    var_names : `list of str` (default: None)
        List of var_names to use for computing the hierarchical clustering. If `var_names` is given,
        then `use_rep` and `n_pcs` is ignored.
    use_raw : `bool`, optional (default: None)

    cor_method : `str`, optional (default: `"pearson"`)
        correlation method to use.
    linkage_method : `str`, optional (default: `"complete"`)
        linkage method to use. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
        for more information.
    key_added : : `str`, optional (default: `"dendrogram"`)
        By default, the dendrogram information is added to the `.uns['dendrogram']`.

    Returns
    -------
    adata.uns['dendrogram'] (or instead of 'dendrogram' the value selected for `key_added`) is updated
    with the dendrogram information

    Examples
    --------

    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, groupby='bulk_labels')
    >>> sc.pl.dendrogram(adata)
    >>> sc.pl.dotplot(adata, ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ'],
    ...               groupby='bulk_labels', dendrogram=True)
    """
    if groupby not in adata.obs_keys():
        raise ValueError('groupby has to be a valid observation. Given value: {}, '
                         'valid observations: {}'.format(groupby, adata.obs_keys()))
    if not is_categorical_dtype(adata.obs[groupby]):
        # if the groupby column is not categorical, turn it into one
        # by subdividing into  `num_categories` categories
        raise ValueError('groupby has to be a categorical observation. Given value: {}, '
                         'Column type: {}'.format(groupby, adata.obs[groupby].dtype))

    if var_names is None:
        rep_df = pd.DataFrame(choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs))
        rep_df.set_index(adata.obs[groupby], inplace=True)
        categories = rep_df.index.categories
    else:
        if use_raw is None and adata.raw is not None: use_raw = True
        gene_names = adata.raw.var_names if use_raw else adata.var_names
        from ..plotting._anndata import _prepare_dataframe
        categories, rep_df = _prepare_dataframe(adata, gene_names, groupby, use_raw)

    # aggregate values within categories using 'mean'
    mean_df = rep_df.groupby(level=0).mean()

    import scipy.cluster.hierarchy as sch

    corr_matrix = mean_df.T.corr(method=cor_method)
    z_var = sch.linkage(corr_matrix, method=linkage_method)
    dendro_info = sch.dendrogram(z_var, labels=categories, no_plot=True)

    # order of groupby categories
    categories_idx_ordered = dendro_info['leaves']

    adata.uns[key_added] = {'linkage': z_var,
                            'groupby': groupby,
                            'use_rep': use_rep,
                            'cor_method': cor_method,
                            'linkage_method': linkage_method,
                            'categories_idx_ordered': categories_idx_ordered,
                            'dendrogram_info': dendro_info}

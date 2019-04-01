"""Calculate overlaps of rank_genes_groups marker genes with marker gene dictionaries
"""
import numpy as np
import pandas as pd

from typing import Union, Optional, Dict
from anndata import AnnData

from .. import logging as logg

def _calc_overlap_count(
    markers1: dict,
    markers2: dict,
):
    """Calculate overlap count between the values of two dictionaries

    Note: dict values must be sets
    """
    overlaps=np.zeros((len(markers1), len(markers2)))

    j=0
    for marker_group in markers1:
        tmp = [len(markers2[i].intersection(markers1[marker_group])) for i in markers2.keys()]
        overlaps[j,:] = tmp
        j += 1

    return overlaps


def _calc_overlap_coef(
    markers1: dict,
    markers2: dict,
):
    """Calculate overlap coefficient between the values of two dictionaries

    Note: dict values must be sets
    """
    overlap_coef=np.zeros((len(markers1), len(markers2)))

    j=0
    for marker_group in markers1:
        tmp = [len(markers2[i].intersection(markers1[marker_group]))/
               max(min(len(markers2[i]), len(markers1[marker_group])),1) for i in markers2.keys()]
        overlap_coef[j,:] = tmp
        j += 1

    return overlap_coef


def _calc_jaccard(
    markers1: dict,
    markers2: dict,
):
    """Calculate jaccard index between the values of two dictionaries

    Note: dict values must be sets
    """
    jacc_results=np.zeros((len(markers1), len(markers2)))

    j=0
    for marker_group in markers1:
        tmp = [len(markers2[i].intersection(markers1[marker_group]))/
               len(markers2[i].union(markers1[marker_group])) for i in markers2.keys()]
        jacc_results[j,:] = tmp
        j += 1

    return jacc_results


def marker_gene_overlap(
    adata: AnnData,
    reference_markers: Union[Dict[str, set], Dict[str,list]],
    *,
    key: str = 'rank_genes_groups',
    method: Optional[str] = 'overlap_count',
    normalize: Union[str, None] = None,
    top_n_markers: Optional[int] = None,
    adj_pval_threshold: Optional[float] = None,
    key_added: Optional[str] = 'marker_gene_overlap',
    inplace: Optional[bool] = False
):
    """Calculate an overlap score between data-deriven marker genes and 
    provided markers

    Marker gene overlap scores can be quoted as overlap counts, overlap 
    coefficients, or jaccard indices. The method returns a pandas dataframe
    which can be used to annotate clusters based on marker gene overlaps.

    This function was written by Malte Luecken.

    Parameters
    ----------
    adata
        The annotated data matrix.
    reference_markers
        A marker gene dictionary object. Keys should be strings with the 
        cell identity name and values are sets or lists of strings which match 
        format of `adata.var_name`.
    key
        The key in `adata.uns` where the rank_genes_groups output is stored.
        By default this is `'rank_genes_groups'`.
    method : `{'overlap_count', 'overlap_coef', 'jaccard'}`, optional 
        (default: `overlap_count`)
        Method to calculate marker gene overlap. `'overlap_count'` uses the
        intersection of the gene set, `'overlap_coef'` uses the overlap 
        coefficient, and `'jaccard'` uses the Jaccard index.
    normalize : `{'reference', 'data', 'None'}`, optional (default: `None`)
        Normalization option for the marker gene overlap output. This parameter
        can only be set when `method` is set to `'overlap_count'`. `'reference'`
        normalizes the data by the total number of marker genes given in the 
        reference annotation per group. `'data'` normalizes the data by the
        total number of marker genes used for each cluster.
    top_n_markers
        The number of top data-derived marker genes to use. By default all 
        calculated marker genes are used. If `adj_pval_threshold` is set along
        with `top_n_markers`, then `adj_pval_threshold` is ignored.
    adj_pval_threshold
        A significance threshold on the adjusted p-values to select marker 
        genes. This can only be used when adjusted p-values are calculated by 
        `sc.tl.rank_genes_groups()`. If `adj_pval_threshold` is set along with 
        `top_n_markers`, then `adj_pval_threshold` is ignored.
    key_added
        Name of the `.uns` field that will contain the marker overlap scores.
    inplace
        Return a marker gene dataframe or store it inplace in `adata.uns`.


    Returns
    -------
    A pandas dataframe with the marker gene overlap scores if `inplace=False`.
    For `inplace=True` `adata.uns` is updated with an additional field 
    specified by the `key_added` parameter (default = 'marker_gene_overlap'). 


    Examples
    --------
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.pp.pca(adata, svd_solver='arpack')
    >>> sc.pp.neighbors(adata)
    >>> sc.tl.louvain(adata)
    >>> sc.tl.rank_genes_groups(adata, groupby='louvain')
    >>> marker_genes = {'CD4 T cells':{'IL7R'},'CD14+ Monocytes':{'CD14', 
    ...                 'LYZ'}, 'B cells':{'MS4A1'}, 'CD8 T cells':{'CD8A'}, 
    ...                 'NK cells':{'GNLY', 'NKG7'}, 'FCGR3A+ Monocytes':
    ...                 {'FCGR3A', 'MS4A7'}, 'Dendritic Cells':{'FCER1A', 
    ...                 'CST3'}, 'Megakaryocytes':{'PPBP'}}
    >>> marker_matches = sc.tl.marker_gene_overlap(adata, marker_genes)
    """
    # Test user inputs
    if inplace:
        raise NotImplementedError('Writing Pandas dataframes to h5ad is '
                                  'currently under development.\n'
                                  'Please use `inplace=False`.')

    if key not in adata.uns:
        raise ValueError('Could not find marker gene data. '
                         'Please run `sc.tl.rank_genes_groups()` first.')

    avail_methods = {'overlap_count', 'overlap_coef', 'jaccard', 'enrich'}
    if method not in avail_methods:
        raise ValueError('Method must be one of {}.'.format(avail_methods))
    
    if normalize == 'None':
        normalize = None

    avail_norm = {'reference', 'data', None}
    if normalize not in avail_norm:
        raise ValueError('Normalize must be one of {}.'.format(avail_norm))
    
    if normalize is not None and method != 'overlap_count':
        raise ValueError('Can only normalize with method=`overlap_count`.')

    if not np.all([isinstance(val, set) for val in reference_markers.values()]):
        try:
            reference_markers = {key:set(val) for key,val
                                 in reference_markers.items()}
        except:
            raise ValueError('Please ensure that `reference_markers` contains '
                             'sets or lists of markers as values.')

    if adj_pval_threshold is not None:
        if 'pvals_adj' not in adata.uns[key]:
            raise ValueError('Could not find adjusted p-value data. '
                             'Please run `sc.tl.rank_genes_groups()` with a '
                             'method that outputs adjusted p-values.')

        if adj_pval_threshold < 0:
            logg.warn('`adj_pval_threshold` was set below 0. '
                      'Threshold will be set to 0.')
            adj_pval_threshold = 0

        if adj_pval_threshold > 1:
            logg.warn('`adj_pval_threshold` was set above 1. '
                      'Threshold will be set to 1.')
            adj_pval_threshold = 1

        if top_n_markers is not None:
            logg.warn('Both `adj_pval_threshold` and `top_n_markers` is set. '
                      '`adj_pval_threshold` will be ignored.')
            
    if top_n_markers is not None:
        if top_n_markers < 1:
            logg.warn('`top_n_markers` was set below 1. '
                      '`top_n_markers` will be set to 1.')
            top_n_markers = 1
            
    # Get data-derived marker genes in a dictionary of sets
    data_markers = dict()
    cluster_ids = adata.uns[key]['names'].dtype.names

    for group in cluster_ids:

        if top_n_markers is not None:
            n_genes = min(top_n_markers, adata.uns[key]['names'].shape[0])
            data_markers[group] = set(adata.uns[key]['names'][group][:n_genes])

        elif adj_pval_threshold is not None:
            n_genes = (adata.uns[key]['pvals_adj'][group] < adj_pval_threshold).sum()
            data_markers[group] = set(adata.uns[key]['names'][group][:n_genes])

            if n_genes == 0:
                logg.warn('No marker genes passed the significance threshold of'
                          ' {} for cluster {!r}.'.format(adj_pval_threshold, 
                                                        group))
        else:
            data_markers[group] = set(adata.uns[key]['names'][group])

    # Find overlaps
    if method == 'overlap_count':
        marker_match = _calc_overlap_count(reference_markers, data_markers)

        if normalize == 'reference':
            # Ensure rows sum to 1
            ref_lengths = np.array([len(reference_markers[m_group]) 
                                    for m_group in reference_markers])
            marker_match = marker_match/ref_lengths[:,np.newaxis]
            marker_match = np.nan_to_num(marker_match)

        elif normalize == 'data':
            #Ensure columns sum to 1
            data_lengths = np.array([len(data_markers[dat_group]) 
                                     for dat_group in data_markers])
            marker_match = marker_match/data_lengths
            marker_match = np.nan_to_num(marker_match)
            
    elif method == 'overlap_coef':
        marker_match = _calc_overlap_coef(reference_markers, data_markers)

    elif method == 'jaccard':
        marker_match = _calc_jaccard(reference_markers, data_markers)
        
    #Note:
    # Could add an 'enrich' option here (fisher's exact test or hypergeometric 
    # test), but that would require knowledge of the size of the space from which
    # the reference marker gene set was taken. This is at best approximately 
    # known. 
    
    # Create a pandas dataframe with the results
    marker_groups = list(reference_markers.keys())
    clusters = list(cluster_ids)
    marker_matching_df = pd.DataFrame(marker_match, index=marker_groups, 
                                      columns=clusters)

    # Store the results
    if inplace: 
        adata.uns[key_added] = marker_matching_df

        logg.hint('added\n'
                  '    \'{}\', marker overlap scores (adata.uns)'
                  .format(key_added))

    else:
        return marker_matching_df

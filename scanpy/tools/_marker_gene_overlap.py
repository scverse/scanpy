"""Calculate overlaps of rank_genes_groups marker genes with marker gene dictionaries
"""
import numpy as np

from typing import Union, Optional
from anndata import AnnData


def _calc_overlap_count(
        markers1: dict,
        markers2: dict):
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
        markers2: dict):
    """Calculate overlap coefficient between the values of two dictionaries

    Note: dict values must be sets
    """
    overlap_coef=np.zeros((len(markers1), len(markers2)))

    j=0
    for marker_group in markers1:
        tmp = [len(markers2[i].intersection(markers1[marker_group]))/
               min(markers2[i], markers1[marker_group]) for i in markers2.keys()]
        overlap_coef[j,:] = tmp
        j += 1

    return overlap_coef


def _calc_jaccard(
        markers1: dict,
        markers2: dict):
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
        key: str,
        reference_markers: dict,
        method: Optional[str] = 'overlap',
        normalize: Union[str, None] = None,
        key_added: Optional[str] = 'marker_gene_overlap'):
    """Calculate an overlap score between data-deriven marker genes and provided markers

    Marker gene overlap scores can be quoted as overlap counts, overlap coefficients, or
    jaccard indices. The method returns a pandas dataframe which can be used to annotate
    clusters based on marker gene overlaps.

    This function was written by Malte Luecken.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    key : `str`, optional (default: `rank_genes_groups`)
        The key in `adata.uns` where the rank_genes_groups output is stored. This field
        should contain a dictionary with a `numpy.recarray()` under the key 'names'.
    reference_markers : `dict`, optional (default: `None`)

    key_added : `str`, optional (default: `None`)

    method : `{'overlap_count', 'overlap_coef', 'jaccard'}`, optional (default: `overlap_count`)

    normalize : `{'reference', 'data', 'None'}`, optional (default: `None`)

    Returns
    -------
    Updates `adata.uns` with an additional field specified by the `key_added`
    parameter (default = 'marker_gene_overlap'). 

    Examples
    --------
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> 
    >>> 
    >>> 
    >>> 
    """
    # Test user inputs
    if key not in adata.uns:
        raise ValueError()

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
        raise ValueError('Please ensure that `reference_markers` contains sets '
                         'of markers as values.')
        

    # Get data-derived marker genes in a dictionary of sets
    data_markers = dict()
    cluster_ids = adata.uns[key]['names'].dtype.names

    for group in cluster_ids:
        data_markers[group] = set(adata.uns[key]['names'][group])

    # To do:
    # - allow to only use the top X marker genes calculated
    # - allow using a p-value cutoff for the genes

    # Find overlaps
    if method == 'overlap_count':
        marker_match = _calc_overlap_count(reference_markers, data_markers)

        if normalize == 'reference':
            # Ensure rows sum to 1
            marker_match = marker_match/marker_match.sum(1)[:,np.newaxis]

        elif noramlize == 'data':
            #Ensure columns sum to 1
            marker_match = marker_match/marker_match.sum(0)
            
    elif method == 'overlap_coef':
        marker_match = _calc_overlap_coef(reference_markers, data_markers)

    elif method == 'jaccard':
        marker_match = _calc_jaccard(reference_markers, data_markers)
        
    #Note:
    # Could add an 'enrich' option here (fisher's exact test or hypergeometric test),
    # but that would require knowledge of the size of the space from which the reference
    # marker gene set was taken. This is at best approximately known.
        
    # Create a pandas dataframe with the results
    marker_groups = list(reference_markers.keys())
    clusters = list(cluster_ids)
    marker_matching_df = pd.DataFrame(marker_match, index=marker_groups, columns=clusters)

    # Store the results
    adata.uns[key_added] = marker_matching_df

    logg.hint('added\n'
              '    \'{}\', marker overlap scores (adata.uns)'.format(key_added))

    return None

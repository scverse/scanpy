# Author: F. Alex Wolf (http://falexwolf.de)
"""Differential Gene Expression Analysis

This is a Beta Version of a tool for differential gene expression testing
between sets detected in previous tools. Tools such as dpt, cluster,...
"""

import numpy as np
from scipy.sparse import issparse
from .. import utils
from .. import logging as logg
from ..preprocessing import simple

def rank_genes_groups(
        adata,
        groupings,
        groups='all',
        n_genes=50,
        compute_distribution=False,
        only_positive=True,
        copy=False):
    """Compare groups by ranking genes according to differential expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    grouping : str
        The key of the sample grouping to consider.
    groups : str, list, np.ndarray, optional (default: 'all')
        Subset of categories - e.g. 'C1,C2,C3' or ['C1', 'C2', 'C3'] - to which
        comparison shall be restricted. If not provided all categories will be
        compared to all other categories.
    n_genes : int (default: 50)
        How many genes to rank by default.
    compute_distribution : bool
        If True, also computes the distribution for top-ranked genes,
        which can be visualized using sc.pl.rank_genes_groups_violin(adata).

    Writes to adata
    ---------------
    rank_genes_groups_zscores : np.ndarray
        Array of shape (number of comparisons) x (number of genes) storing the
        zscore of the each gene for each test.
    rank_genes_groups_rankings_names : np.ndarray of dtype str
        Array of shape (number of comparisons). Stores the labels for each comparison,
        for example "C1 vs. C2" when comparing category 'C1' with 'C2'.
    rank_genes_groups_rankings_geneidcs : np.ndarray
        Array of shape (number of comparisons) x (number of genes) storing gene
        indices that sort them according to decreasing absolute value of the
        zscore.
    """
    logg.m('find differentially expressed genes', r=True)
    adata = adata.copy() if copy else adata
    n_genes_user = n_genes
    utils.check_adata(adata)
    # for clarity, rename variable
    group_key = groupings
    groups_names = groups
    if isinstance(groups_names, list) and isinstance(groups_names[0], int):
        groups_names = [str(n) for n in groups_names]
    groups_names, groups_masks = utils.select_groups(adata, groups_names, group_key)
    adata.add['rank_genes_groups'] = group_key
    adata.add['rank_genes_groups_names'] = groups_names
    X = adata.X

    # loop over all masks and compute means, variances and sample numbers
    n_groups = groups_masks.shape[0]
    n_genes = X.shape[1]
    means = np.zeros((n_groups, n_genes))
    vars = np.zeros((n_groups, n_genes))
    ns = np.zeros(n_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        means[imask], vars[imask] = simple._get_mean_var(X[mask])
        ns[imask] = np.where(mask)[0].size
    logg.info('... consider "{}":'.format(group_key), groups_names,
              'with sample numbers', ns)

    # test each group against the rest of the data
    rankings_gene_zscores = []
    rankings_gene_names = []
    reference_indices = np.arange(adata.n_vars, dtype=int)
    for igroup in range(n_groups):
        mask_rest = ~groups_masks[igroup]
        mean_rest, var_rest = simple._get_mean_var(X[mask_rest])
        # Make a more conservative assumption on the variance reduction
        # in the reference. Instead of this
        # ns_rest = np.where(mask)[0].size
        # use this
        ns_rest = ns[igroup]
        denominator = np.sqrt(vars[igroup]/ns[igroup] + var_rest/ns_rest)
        denominator[np.flatnonzero(denominator == 0)] = np.nan
        zscores = (means[igroup] - mean_rest) / denominator
        zscores[np.isnan(zscores)] = 0
        zscores = zscores if only_positive else np.abs(zscores)
        partition = np.argpartition(zscores, -n_genes_user)[-n_genes_user:]
        partial_indices = np.argsort(zscores[partition])[::-1]
        global_indices = reference_indices[partition][partial_indices]
        rankings_gene_zscores.append(zscores[global_indices])
        rankings_gene_names.append(adata.var_names[global_indices])
        if compute_distribution:
            mask = groups_masks[igroup]
            for gene_counter in range(n_genes_user):
                gene_idx = global_indices[gene_counter]
                X_col = X[mask, gene_idx]
                if issparse(X): X_col = X_col.toarray()[:, 0]
                identifier = _build_identifier(group_key, groups_names[igroup],
                                               gene_counter, adata.var_names[gene_idx])
                full_col = np.empty(adata.n_smps)
                full_col[:] = np.nan
                full_col[mask] = (X_col - mean_rest[gene_idx])/denominator[gene_idx]
                adata.smp[identifier] = full_col
                
    adata.add['rank_genes_groups_gene_scores'] = np.rec.fromarrays(
        [n for n in rankings_gene_zscores],
        dtype=[(rn, 'float32') for rn in groups_names])
    adata.add['rank_genes_groups_gene_names'] = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, 'U50') for rn in groups_names])
    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "rank_genes_groups_gene_names", np.recarray to be indexed by the `groups` (adata.add)\n'
           '    "rank_genes_groups_gene_zscores", the scores (adata.add)\n'
           '    "rank_genes_...", distributions of top-ranked genes (adata.smp)')
    return adata if copy else None


def _build_identifier(group_key, name, gene_counter, gene_name):
    return 'rank_genes_{}_{}_{}_{}'.format(
        group_key, name, gene_counter, gene_name)

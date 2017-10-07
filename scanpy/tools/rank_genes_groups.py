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
        groupby,
        groups='all',
        group_reference=None,
        n_genes=100,
        compute_distribution=False,
        only_positive=True,
        copy=False):
    """Rank genes according to differential expression [Wolf17]_.

    Rank genes by differential expression. By default, a t-test-like ranking is
    used, in which means are normalized with variances. Soon, a Wilcoxon-rank
    test and other alternatives will be provided.

    Parameters
    ----------
    adata : `AnnData`
        Annotated data matrix.
    groupby : `str`
        The key of the sample grouping to consider.
    groups : `str`, `list`, optional (default: `'all'`)
        Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which comparison shall
        be restricted. If not passed, a ranking will be generated for all
        groups.
    group_reference : `str` or `None`, optional (default: `None`)
        If `None`, compare each group to the union of the rest of the group.  If
        a group identifier, the comparison will be with respect to this group.
    n_genes : `int` (default: 100)
        How many genes to rank by default.
    compute_distribution : `bool`
        If `True`, also computes the distribution for top-ranked genes, which
        can be visualized using `sc.pl.rank_genes_groups_violin(adata)`.

    Returns
    -------
    rank_genes_groups_gene_zscores : np.ndarray of dtype float (adata.add)
        Array of shape (number of comparisons) × (number of genes) storing the
        zscore of the each gene for each test.
    rank_genes_groups_gene_names : np.ndarray of dtype str (adata.add)
        Array of shape (number of comparisons). Stores the labels for each comparison,
        for example "C1 vs. C2" when comparing category 'C1' with 'C2'.
    """
    logg.info('find differentially expressed genes', r=True)
    adata = adata.copy() if copy else adata
    n_genes_user = n_genes
    utils.check_adata(adata)
    # for clarity, rename variable
    groups_order = groups
    if isinstance(groups_order, list) and isinstance(groups_order[0], int):
        groups_order = [str(n) for n in groups_order]
    if group_reference is not None and group_reference not in set(groups_order):
        groups_order += [group_reference]
    if (group_reference is not None
        and group_reference not in set(adata.add[groupby + '_order'])):
        raise ValueError('group_reference = {} needs to be one of groupby = {}.'
                         .format(group_reference, groupby))
    groups_order, groups_masks = utils.select_groups(
        adata, groups_order, groupby)
    adata.add['rank_genes_groups'] = groupby
    adata.add['rank_genes_groups_order'] = groups_order
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
    logg.info('... consider "{}":'.format(groupby), groups_order,
              'with sample numbers', ns)

    if group_reference is not None:
        ireference = np.where(groups_order == group_reference)[0][0]
    
    # test each either against the union of all other groups
    # or against a specific group
    rankings_gene_zscores = []
    rankings_gene_names = []
    reference_indices = np.arange(adata.n_vars, dtype=int)
    for igroup in range(n_groups):
        if group_reference is None:
            mask_rest = ~groups_masks[igroup]
        else:
            if igroup == ireference: continue
            else: mask_rest = groups_masks[ireference]
        mean_rest, var_rest = simple._get_mean_var(X[mask_rest])
        # Make a more conservative assumption on the variance reduction
        # in the reference. Instead of this
        ns_rest = np.where(mask_rest)[0].size
        # use this
        # ns_rest = ns[igroup]
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
                identifier = _build_identifier(groupby, groups_order[igroup],
                                               gene_counter, adata.var_names[gene_idx])
                full_col = np.empty(adata.n_smps)
                full_col[:] = np.nan
                full_col[mask] = (X_col - mean_rest[gene_idx])/denominator[gene_idx]
                adata.smp[identifier] = full_col

    groups_order_save = groups_order
    if group_reference is not None:
        groups_order_save = [g for g in groups_order if g != group_reference]
                
    adata.add['rank_genes_groups_gene_scores'] = np.rec.fromarrays(
        [n for n in rankings_gene_zscores],
        dtype=[(rn, 'float32') for rn in groups_order_save])
    adata.add['rank_genes_groups_gene_names'] = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, 'U50') for rn in groups_order_save])
    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "rank_genes_groups_gene_names", np.recarray to be indexed by the `groups` (adata.add)\n'
           '    "rank_genes_groups_gene_zscores", the scores (adata.add)\n'
           '    "rank_genes_...", distributions of top-ranked genes (adata.smp)')
    return adata if copy else None


def _build_identifier(groupby, name, gene_counter, gene_name):
    return 'rank_genes_{}_{}_{}_{}'.format(
        groupby, name, gene_counter, gene_name)

# Author: F. Alex Wolf (http://falexwolf.de)
"""Differential Gene Expression Analysis

This is a Beta Version of a tool for differential gene expression testing
between sets detected in previous tools. Tools such as dpt, cluster,...
"""

import numpy as np
from itertools import combinations
from scipy.stats.distributions import norm
from scipy.sparse import issparse
from .. import utils
from .. import logging as logg
from .. import settings as sett
from ..preprocessing import simple

def rank_genes_groups(
        adata,
        groupings,
        compute_distribution=False,
        groups='all',
        n_genes=30,
        pairings=False,
        only_positive=None,
        sig_level=0.05,
        correction='Bonferroni',
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
    n_genes : int (default: 30)
        How many genes to rank by default.
    compute_distribution : bool
        If True, also computes the distribution for top-ranked genes,
        which can be visualized using sc.pl.rank_genes_groups_violin(adata).
    pairings : bool, optional (default: False)
        Test pairings of groups instead of group vs. the rest of the data.

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
        # means[imask] = X[mask].mean(axis=0)
        # vars[imask] = X[mask].var(axis=0)
        ns[imask] = np.where(mask)[0].size
    logg.info('... consider "{}":'.format(group_key), groups_names,
              'with sample numbers', ns)

    # each test provides a ranking of genes
    # we store the name of the ranking, i.e. the name of the test,
    # in the following list
    adata.add['rank_genes_groups_rankings_names'] = []
    distribution_per_gene = np.recarray((adata.n_smps, n_groups*n_genes_user),
                                        dtype=[('{}_{}'.format(group, gene), np.float32)
                                               for group in range(n_groups)
                                               for gene in range(n_genes_user)])

    # test each group against the rest of the data
    if not pairings:

        only_positive = True if only_positive is None else only_positive

        zscores_all = np.zeros((n_groups, n_genes))
        rankings_geneidcs = np.zeros((n_groups, n_genes), dtype=int)

        for igroup in range(n_groups):
            mask_rest = ~groups_masks[igroup]
            mean_rest, var_rest = simple._get_mean_var(X[mask_rest])
            ns_rest = np.where(mask)[0].size
            # z-scores
            denom = np.sqrt(vars[igroup]/ns[igroup] + var_rest/ns_rest)
            zeros = np.flatnonzero(denom == 0)
            denom[zeros] = np.nan
            zscores = (means[igroup] - mean_rest) / denom
            # the following is equivalent with
            # zscores = np.ma.masked_invalid(zscores)
            zscores = np.ma.masked_array(zscores, mask=np.isnan(zscores))

            zscores_all[igroup] = zscores
            abs_zscores = zscores if only_positive else np.abs(zscores)

            # sort genes according to score
            ranking_geneidcs = np.argsort(abs_zscores)[::-1]
            # move masked values to the end of the index array
            masked = abs_zscores[ranking_geneidcs].mask
            len_not_masked = len(ranking_geneidcs[masked == False])
            save_masked_idcs = np.copy(ranking_geneidcs[masked])
            ranking_geneidcs[:len_not_masked] = ranking_geneidcs[masked == False]
            ranking_geneidcs[len_not_masked:] = save_masked_idcs
            # write to global rankings_genedics
            rankings_geneidcs[igroup] = ranking_geneidcs
            if compute_distribution:
                mask = groups_masks[igroup]
                for gene_counter in range(n_genes_user):
                    gene_idx = ranking_geneidcs[gene_counter]
                    X_col = X[mask, gene_idx]
                    if issparse(X): X_col = X_col.toarray()[:, 0]
                    # distribution_per_gene[:, igroup*n_genes_user+gene_counter] = X_col - mean_rest[gene_idx]
                    identifier = _build_identifier(group_key, groups_names[igroup],
                                                   gene_counter, adata.var_names[gene_idx])
                    full_col = np.empty(adata.n_smps)
                    full_col[:] = np.NAN
                    full_col[mask] = X_col - mean_rest[gene_idx]
                    adata.smp[identifier] = full_col
            # names
            ranking_name = groups_names[igroup]
            adata.add['rank_genes_groups_rankings_names'].append(ranking_name)

    # test all group pairings # TODO: needs to be updated
    else:
        logg.warn('`pairings` option to be used with care!')
        igroups_masks = np.arange(len(groups_masks), dtype=int)
        pairs = list(combinations(igroups_masks, 2))
        pvalues_all = np.zeros((len(pairs), n_genes))
        zscores_all = np.zeros((len(pairs), n_genes))
        rankings_geneidcs = np.zeros((len(pairs), n_genes), dtype=int)

        # test all combinations of groups against each other
        for ipair, (i, j) in enumerate(pairs):
            # z-scores
            denom = np.sqrt(vars[i]/ns[i] + vars[j]/ns[j])
            zeros = np.flatnonzero(denom == 0)
            denom[zeros] = np.nan
            zscores = (means[i] - means[j]) / denom
            # the following is equivalent with
            # zscores = np.ma.masked_invalid(zscores)
            zscores = np.ma.masked_array(zscores, mask=np.isnan(zscores))

            zscores_all[ipair] = zscores
            abs_zscores = np.abs(zscores)

            # sort genes according to score
            ranking_geneidcs = np.argsort(abs_zscores)[::-1]
            # move masked values to the end of the index array
            masked = abs_zscores[ranking_geneidcs].mask
            len_not_masked = len(ranking_geneidcs[masked == False])
            save_masked_idcs = np.copy(ranking_geneidcs[masked])
            ranking_geneidcs[:len_not_masked] = ranking_geneidcs[masked == False]
            ranking_geneidcs[len_not_masked:] = save_masked_idcs
            # write to global rankings_genedics
            rankings_geneidcs[ipair] = ranking_geneidcs
            # names
            ranking_name = groups_names[i] + ' vs ' + groups_names[j]
            adata.add['rank_genes_groups_rankings_names'].append(ranking_name)

    adata.add['rank_genes_groups_zscores'] = zscores_all
    adata.add['rank_genes_groups_rankings_geneidcs'] = rankings_geneidcs
    adata.add['rank_genes_groups_scoreskey'] = 'zscores'

    names_of_top_ranked_genes = []
    for irank in range(len(adata.add['rank_genes_groups_rankings_names'])):
        names_of_top_ranked_genes.append([])
        for ig, g in enumerate(rankings_geneidcs[irank, :n_genes_user]):
            names_of_top_ranked_genes[-1].append(adata.var_names[g])
    adata.add['rank_genes_groups_names_of_top_ranked_genes'] = np.rec.fromarrays(
        [n for n in names_of_top_ranked_genes],
        dtype=[(rn, 'U50') for rn in adata.add['rank_genes_groups_rankings_names']])

    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "rank_genes_groups_names_of_top_ranked_genes", a np.recarray of top-ranked genes to be indexed by the str keys in `names` (adata.add)\n'
           '    "rank_genes_...", distributions of top-ranked genes (adata.smp)')
    return adata if copy else None


def _build_identifier(group_key, name, gene_counter, gene_name):
    return 'rank_genes_{}_{}_{}_{}'.format(
        group_key, name, gene_counter, gene_name)

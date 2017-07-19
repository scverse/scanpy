# Author: F. Alex Wolf (http://falexwolf.de)
"""Differential Gene Expression Analysis

This is a Beta Version of a tool for differential gene expression testing
between sets detected in previous tools. Tools such as dpt, cluster,...
"""

import numpy as np
from itertools import combinations
from scipy.stats.distributions import norm
from .. import utils
from .. import logging as logg
from .. import settings as sett
from ..preprocessing import simple

def diffrank(adata,
             key,
             names='all',
             n_genes=30,
             pairings=False,
             only_positive=None,
             sig_level=0.05,
             correction='Bonferroni',
             log=False,
             copy=False):
    """Compare groups by ranking genes according to differential expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    key : str
        The key of the sample grouping to consider.
    names : str, list, np.ndarray, optional (default: 'all')
        Subset of categories - e.g. 'C1,C2,C3' or ['C1', 'C2', 'C3'] - to which
        comparison shall be restricted. If not provided all categories will be
        compared to all other categories.
    pairings : bool, optional (default: False)
        Test pairings of groups instead of group vs. the rest of the data.

    Writes to adata
    ---------------
    diffrank_zscores : np.ndarray
        Array of shape (number of comparisons) x (number of genes) storing the
        zscore of the each gene for each test.
    diffrank_rankings_names : np.ndarray of dtype str
        Array of shape (number of comparisons). Stores the labels for each comparison,
        for example "C1 vs. C2" when comparing category 'C1' with 'C2'.
    diffrank_rankings_geneidcs : np.ndarray
        Array of shape (number of comparisons) x (number of genes) storing gene
        indices that sort them according to decreasing absolute value of the
        zscore.
    """
    logg.m('find differentially expressed genes', r=True)
    adata = adata.copy() if copy else adata
    n_genes_user = n_genes
    utils.check_adata(adata)
    # for clarity, rename variable
    groups_names = names
    groups_names, groups_masks = utils.select_groups(adata, groups_names, key)
    adata.add['diffrank_groups'] = key
    adata.add['diffrank_groups_names'] = groups_names
    X = adata.X
    if log:
        # TODO: treat negativity explicitly
        X = np.abs(X)
        X = np.log(X) / np.log(2)

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
    logg.info('consider "{}":'.format(key), groups_names,
              'with sample numbers', ns)
    sett.m(2, 'means', means)
    sett.m(2, 'variances', vars)

    # each test provides a ranking of genes
    # we store the name of the ranking, i.e. the name of the test,
    # in the following list
    adata.add['diffrank_rankings_names'] = []

    # test each group against the rest of the data
    if not pairings:

        only_positive = True if only_positive is None else only_positive

        zscores_all = np.zeros((n_groups, n_genes))
        rankings_geneidcs = np.zeros((n_groups, n_genes), dtype=int)

        for igroup in range(n_groups):
            mask = ~groups_masks[igroup]
            mean_rest, var_rest = simple._get_mean_var(X[mask])
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

            # p-values
            if False:
                pvalues = 2 * norm.sf(abs_zscores)  # two-sided test
                pvalues = np.ma.masked_invalid(pvalues)
                sig_genes = np.flatnonzero(pvalues < 0.05/zscores.shape[0])
                pvalues_all[igroup] = pvalues

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
            # names
            ranking_name = groups_names[igroup]
            adata.add['diffrank_rankings_names'].append(ranking_name)

    # test all group pairings
    else:
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

            # p-values
            if False:
                pvalues = 2 * norm.sf(abs_zscores)  # two-sided test
                pvalues = np.ma.masked_invalid(pvalues)
                sig_genes = np.flatnonzero(pvalues < 0.05/zscores.shape[0])
                pvalues_all[ipair] = pvalues

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
            adata.add['diffrank_rankings_names'].append(ranking_name)

    if False:
        adata.add['diffrank_pvalues'] = -np.log10(pvalues_all)

    adata.add['diffrank_zscores'] = zscores_all
    adata.add['diffrank_rankings_geneidcs'] = rankings_geneidcs
    adata.add['diffrank_scoreskey'] = 'zscores'

    names_of_top_ranked_genes = []
    for irank in range(len(adata.add['diffrank_rankings_names'])):
        names_of_top_ranked_genes.append([])
        for ig, g in enumerate(rankings_geneidcs[irank, :n_genes_user]):
            names_of_top_ranked_genes[-1].append(adata.var_names[g])

    adata.add['diffrank_names_of_top_ranked_genes'] = np.array(names_of_top_ranked_genes)

    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "diffrank_names_of_top_ranked_genes", the ordered names of topranked genes (adata.add)\n'
           '    "diffrank_zscores", the rankings for each genes (adata.add)')
    return adata if copy else None

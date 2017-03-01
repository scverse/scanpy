# Author: F. Alex Wolf (http://falexwolf.de)
"""
Differential Gene Expression Analysis

This is a Beta Version of a tool for differential gene expression testing
between sets detected in previous tools. Tools such as dpt, cluster,...
"""   

from itertools import combinations
import numpy as np
from scipy.stats.distributions import norm
from .. import utils
from .. import plotting as plott
from .. import settings as sett

def diffrank(adata,
             smp='groups',
             names='all',
             sig_level=0.05,
             correction='Bonferroni',
             log=False):
    """
    Compare groups by ranking genes according to differential expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: 'exp_groups')
        Specify the name of the grouping to consider.
    names : str, list, np.ndarray, optional (default: 'all')
        Subset of categories - e.g. 'C1,C2,C3' or ['C1', 'C2', 'C3'] - to which
        comparison shall be restricted. If not provided all categories will be
        compared to all other categories.

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
    # for clarity, rename variable
    groups_names = names
    groups_names, groups_masks = utils.select_groups(adata, groups_names, smp)
    adata['diffrank_groups'] = smp
    adata['diffrank_groups_names'] = groups_names
    X = adata.X
    if log:
        # TODO: treat negativity explicitly
        X = np.abs(X)
        X = np.log(X) / np.log(2)

    # loop over all masks and compute means, variances and sample numbers
    nr_groups = groups_masks.shape[0]
    nr_genes = X.shape[1]
    means = np.zeros((nr_groups, nr_genes))
    vars = np.zeros((nr_groups, nr_genes))
    ns = np.zeros(nr_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        means[imask] = X[mask].mean(axis=0)
        vars[imask] = X[mask].var(axis=0)
        ns[imask] = np.where(mask)[0].size
    sett.m(0, 'testing', smp, groups_names, 'with sample numbers', ns)
    sett.m(2, 'means', means) 
    sett.m(2, 'variances', vars)

    igroups_masks = np.arange(len(groups_masks), dtype=int)
    pairs = list(combinations(igroups_masks, 2))
    pvalues_all = np.zeros((len(pairs), nr_genes))
    zscores_all = np.zeros((len(pairs), nr_genes))
    rankings_geneidcs = np.zeros((len(pairs), nr_genes),dtype=int)
    # each test provides a ranking of genes
    # we store the name of the ranking, i.e. the name of the test, 
    # in the following list
    adata['diffrank_rankings_names'] = []
    
    # test all combinations of groups against each other
    for ipair, (i,j) in enumerate(pairs):
        # z-scores
        denom = np.sqrt(vars[i]/ns[i] + vars[j]/ns[j])
        zeros = np.flatnonzero(denom==0)
        denom[zeros] = np.nan
        zscores = (means[i] - means[j]) / denom
        # the following is equivalent with 
        # zscores = np.ma.masked_invalid(zscores)
        zscores = np.ma.masked_array(zscores, mask=np.isnan(zscores))
        
        zscores_all[ipair] = zscores
        abs_zscores = np.abs(zscores)

        # p-values
        if False:
            pvalues = 2 * norm.sf(abs_zscores) # two-sided test
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
        ranking_name = groups_names[i] + ' vs '+ groups_names[j]
        adata['diffrank_rankings_names'].append(ranking_name)

    if False:
        adata['diffrank_pvalues'] = -np.log10(pvalues_all)

    adata['diffrank_zscores'] = zscores_all
    adata['diffrank_rankings_geneidcs'] = rankings_geneidcs
    adata['diffrank_scoreskey'] = 'zscores'

    return adata

def plot_diffrank(adata, n_genes=20):
    """
    Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_genes : int
        Number of genes to show.
    """
    plott.ranking(adata, toolkey='diffrank', n_genes=n_genes)
    writekey = sett.basekey + '_diffrank_' + adata['diffrank_groups'] + sett.plotsuffix
    plott.savefig(writekey)
    if not sett.savefigs and sett.autoshow:
        from ..compat.matplotlib import pyplot as pl
        pl.show()


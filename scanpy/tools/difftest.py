# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
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

def difftest(dgroups, ddata=None,
             groups_names='all',
             sig_level=0.05,
             correction='Bonferroni',
             log=True):
    """
    Perform differential gene expression test for groups defined in dgroups.

    Parameters
    ----------
    dgroups (or ddata) : dict containing
        groups_names : list, np.ndarray of dtype str
            Array of shape (number of groups) that names the groups.
        groups : list, np.ndarray of dtype str
            Array of shape (number of samples) that names the groups.
    ddata : dict, optional
        Data dictionary containing gene names.
    groups_names : list, np.ndarray of dtype str
        Subset of names in dgroups['groups_names'] to which comparison shall be
        restricted.

    Returns
    -------
    ddifftest : dict containing
        zscores : np.ndarray
            Array of shape (number of tests) x (number of genes) storing the
            zscore of the each gene for each test.
        testlabels : np.ndarray of dtype str
            Array of shape (number of tests). Stores the labels for each test.
        genes_sorted : np.ndarray
            Array of shape (number of tests) x (number of genes) storing genes
            sorted according the decreasing absolute value of the zscore.
    """    
    params = locals(); del params['ddata']; del params['dgroups']
    # if ddata is empty, assume that ddata dgroups also contains
    # the data file elements
    if not ddata:
        sett.m(0, 'testing experimental groups')
        ddata = dgroups
    # TODO: treat negativity explicitly
    X = np.abs(ddata['X'])
    # Convert X to log scale
    if params['log']:
        XL = np.log(X) / np.log(2)
    else:
        XL = X

    # select groups
    groups_names, groups_masks = utils.select_groups(dgroups, groups_names)
    sett.m(0, 'testing groups', groups_names)

    # loop over all groups_masks and compute means, variances and sample numbers
    # in groups_masks
    means = np.zeros((groups_masks.shape[0],X.shape[1]))
    vars = np.zeros((groups_masks.shape[0],X.shape[1]))
    ns = np.zeros(groups_masks.shape[0],dtype=int)
    for igroup, group in enumerate(groups_masks):
        means[igroup] = XL[group].mean(axis=0)
        vars[igroup] = XL[group].var(axis=0)
        ns[igroup] = np.where(group)[0].size

    ddifftest = {'type' : 'difftest'}
    igroups_masks = np.arange(len(groups_masks),dtype=int)
    pairs = list(combinations(igroups_masks,2))
    pvalues_all = np.zeros((len(pairs), X.shape[1]))
    zscores_all = np.zeros((len(pairs), X.shape[1]))
    genes_sorted = np.zeros((len(pairs), X.shape[1]),dtype=int)
    ddifftest['testnames'] = []
    
    # test all combinations of groups against each other
    for ipair, (i,j) in enumerate(pairs):
        # z-scores
        denom = np.sqrt(vars[i]/ns[i] + vars[j]/ns[j])
        zeros = np.flatnonzero(denom==0)
        denom[zeros] = np.nan
        zscores = (means[i] - means[j]) / denom
        zscores = np.ma.masked_invalid(zscores)
        zscores_all[ipair] = zscores

        abszscores = np.abs(zscores)
        # p-values
        if False:
            pvalues = 2*norm.sf(abszscores) # two-sided test
            pvalues = np.ma.masked_invalid(pvalues)
            sig_genes = np.flatnonzero(pvalues < 0.05/zscores.shape[0])
            pvalues_all[ipair] = pvalues

        # sort genes according to score
        genes_sorted[ipair] = np.argsort(abszscores)[::-1]

        # names
        testlabel = groups_names[i] + ' vs '+ groups_names[j]
        ddifftest['testnames'].append(testlabel)

    if False:
        ddifftest['pvalues'] = -np.log10(pvalues_all)

    ddifftest['zscores'] = zscores_all
    ddifftest['genes_sorted'] = genes_sorted
    ddifftest['scoreskey'] = 'zscores'
    return ddifftest

def plot(ddifftest, ddata, params=None):
    """
    Plot ranking of genes for all tested comparisons.
    """
    plott.ranking(ddifftest, ddata)
    plott.savefig(ddifftest['writekey'])
    if not sett.savefigs and sett.autoshow:
        from ..compat.matplotlib import pyplot as pl
        pl.show()


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

def difftest(dprev, ddata=None,
             groups_ids='all',
             groups_names='all',
             sig_level=0.05,
             correction='Bonferroni',
             log=True):
    """
    Perform differential gene expression test for groups defined in dprev.

    Parameters
    ----------
    dprev (or ddata) : dict containing
        groups_names : np.ndarray of dtype str
            Array of shape (number of groups) that names the groups.
        groups : np.ndarray of dtype str
            Array of shape (number of samples) that names the groups.
    ddata : dict
        Data dictionary containing gene names.
    params: dict, optional. possible keys are
        groups_names : list of str or int
            Subset of names in dprev['groups_names'] to which comparison shall be
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
    params = locals(); del params['ddata']; del params['dprev']
    # if ddata is empty, assume that ddata dprev also contains
    # the data file elements
    if not ddata:
        ddata = dprev
    # TODO: treat negativity explicitly
    X = np.abs(ddata['X'])
    # Convert X to log scale
    if params['log']:
        XL = np.log(X) / np.log(2)
    else:
        XL = X

    # select groups
    groups_names = dprev['groups_names']
    groups_ids = list(range(len(groups_names)))
    groups_masks = dprev['groups_masks']
    if params['groups_names'] != 'all':
        groups_names = np.array(params['groups_names'])
        groups_ids = np.where(np.in1d(dprev['groups_names'], groups_names))[0]
        if not np.any(groups_ids):
            sett.m(0, 'specify valid groups_names for testing, one of',
                   dprev['groups_names'])
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
    sett.m(0, 'testing groups', groups_names, 'with ids', groups_ids)

    # loop over all groups_masks and compute means, variances and sample numbers in groups_masks
    means = np.zeros((groups_masks.shape[0],X.shape[1]))
    vars = np.zeros((groups_masks.shape[0],X.shape[1]))
    ns = np.zeros(groups_masks.shape[0],dtype=int)
    for igroup,group in enumerate(groups_masks):
        means[igroup] = XL[group].mean(axis=0)
        vars[igroup] = XL[group].var(axis=0)
        ns[igroup] = np.where(group)[0].size

    ddifftest = {'type' : 'difftest'}
    igroups_masks = np.arange(len(groups_masks),dtype=int)
    pairs = list(combinations(igroups_masks,2))
    pvalues_all = np.zeros((len(pairs),X.shape[1]))
    zscores_all = np.zeros((len(pairs),X.shape[1]))
    genes_sorted = np.zeros((len(pairs),X.shape[1]),dtype=int)
    ddifftest['testnames'] = []
    
    # test all combinations of groups against each other
    for ipair,(i,j) in enumerate(pairs):
        # z-scores
        denom = np.sqrt(vars[i]/ns[i] + vars[j]/ns[j])
        zeros = np.flatnonzero(denom==0)
        denom[zeros] = np.nan
        zscores = (means[i] - means[j]) / denom
        zscores = np.ma.masked_invalid(zscores)
        zscores_all[ipair] = zscores

        abszscores = np.abs(zscores)
        # pvalues
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
    from ..compat.matplotlib import pyplot as pl
    plott.ranking(ddifftest, ddata)
    if not sett.savefigs:
        pl.show()
    else:
        pl.savefig(sett.figdir+ddifftest['writekey']+'.'+sett.extf)

    
def difftest_shedden(ddata=None, params=None):
    """
    Perform differential gene expression test for groups defined in ddata.

    Parameters
    ----------
    dgroup : dict from tool containing at least a key
        groups : array of index arrays that denotes subgroups

    Returns
    -------
    ddifftest : dict containing
        indices : ....
    
    Note
    ----
    The function is based on a script by Kerby Shedden.
    http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html
    """

    if len(params) == 0:
        params = default_params

    X = ddata['X']
    GID = ddata['colnames']
    STP = ddata['STP']
    UC = ddata['UC']
    CD = ddata['CD']

    # Convert X to log scale
    XL = np.log(X) / np.log(2)

    # Z-test statistics
    MUC = XL[UC].mean(axis=0) # means in ulcerative colitis samples
    MCD = XL[CD].mean(axis=0) # means in Crohn's disease samples

    VUC = XL[UC].var(axis=0)  # variances in ulcerative colitis samples
    VCD = XL[CD].var(axis=0)  # variances in Crohn's disease samples

    nUC = len(UC)          # Number of ulcerative colitis samples
    nCD = len(CD)          # Number of Crohn's disease samples

    Z = (MUC - MCD) / np.sqrt(VUC/nUC + VCD/nCD)

    # Split the sample in half based on marginal standard deviation of
    # gene expression.
    SD = X.std(axis=0)

    isdl = np.flatnonzero(SD < np.median(SD))
    isdh = np.flatnonzero(SD > np.median(SD))

    # Split the sample in half based on the marginal mean gene expression
    SD = X.mean(axis=0)
    imnl = np.flatnonzero(SD < np.median(SD))
    imnh = np.flatnonzero(SD > np.median(SD))

    # Z-score threshold under Bonferroni correction
    zst = -norm.ppf(params['sig_level']/2/Z.shape[0])

    # Find the genes that meet this condition
    ii = np.flatnonzero(np.abs(Z) > zst)

    with open("out/bonferroni_genes.csv", "w") as OUT:
        for i in ii:
            OUT.write("%.4f,%s\n" % (Z[i], GID[i]))

    ddifftest = {}
    ddifftest['indices'] = ii

    print(Z.mean())
    print(len(ii))

    return ddifftest
    
# def add_args(p):
#     """
#     Update parser.
#     """
#     # dictionary for adding arguments
#     dadd_args = {
#         '--prev' : {
#             'type' : str,
#             'default' : 'dpt',
#             'help' : 'Specify the "previous" tool ' 
#                      '- the one you used to generate subgroups of the '
#                      'data. For example, "scct" or "dpt" (default: dpt).'
#             }
#         }
#     p = utils.add_args(p, dadd_args)
#     return p

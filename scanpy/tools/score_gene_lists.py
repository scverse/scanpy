"""Calculate scores based on the expression of gene lists.
"""

import numpy as np
import pandas as pd
import scipy.sparse
from .. import logging as logg


def score_gene_list(
        adata,
        gene_list,
        gene_pool=None,
        n_bins=25,
        ctrl_size=50,
        score_name='Score',
        random_state=0,
        copy=False):  # we use the scikit-learn convention of calling the seed "random_state"
    """Score a set of genes.

    This function calculates a score using a reference of control genes
    sampled matching bins of expression

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    gene_list : iterable
        The list of gene names used for score calculation.
    gene_pool : `list` or `None`, optional (default: `None`)
        A list of genes to be sampled as control set.
    n_bins : `int`, optional (default: 25)
        Number of expression level cuts for sampling.
    ctrl_size : `int`, optional (default: 100)
        Number of genes to be sampled for each bin of expression.
    score_name : `str`, optional (default: `Score`)
        Name of the slot to be added in obs.
    random_state : `int`, optional (default: 0)
        Change random seed.
    copy : `bool`, optional (default: `False`)
        Copy `adata` or modify it inplace.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with an additional field
    `score_name`.
    """
    logg.info('Adding score', r=True)

    if random_state:
        np.random.seed(random_state)

    adata = adata.copy() if copy else adata

    gene_list = set([x for x in gene_list if x in adata.var_names])

    if not gene_list:
        logg.error('A gene list must be passed', r=True)

    if not gene_pool:
        gene_pool = list(adata.var_names)
    else:
        gene_pool = [x for x in gene_pool if x in adata.var_names]

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes in our cells against random genes in a
    # matched interval of expression...

    if scipy.sparse.issparse(adata.X):
        obs_avg = pd.Series(
            np.nanmean(
                adata[:, gene_pool].X.toarray(), axis=0), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(adata[:, gene_pool].X, axis=0), index=gene_pool)  # average expression of genes

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick 100 genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))  # if ctrl_size > len(r_genes) is not a problem for numpy...

    # To index, we need a list - indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)

    score = np.mean(adata[:, gene_list].X, axis=1) - np.mean(adata[:, control_genes].X, axis=1)
    adata.obs[score_name] = pd.Series(np.array(score).ravel(), index=adata.obs_names)

    return adata if copy else None


def score_cell_cycle_genes(
        adata,
        s_genes,
        g2m_genes,
        copy=False):
    """Score cell cycle genes.

    Given two lists of genes associated to S phase and G2M phase, calculates scores
    and assigns a cell cycle phase (G1, S or G2M).

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    s_genes : `list`
        List of genes associated to S phase
    g2m_genes : `list`
        List of genes associated to G2M phase
    copy : `bool`, optional (default: `False`)
        Copy adata or modify it inplace.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    S_score : `pd.Series` (``adata.obs``, dtype `object`)
        Array of dim (number of samples) that stores the score for S phase for each cell.
    G2M_score : `pd.Series` (``adata.obs``, dtype `object`)
        Array of dim (number of samples) that stores the score for G2M phase for each cell.
    Phase : `pd.Series` (``adata.obs``, dtype `object`)
        Array of dim (number of samples) that stores the phase (`S`, `G2M` or `G1`)
        for each cell.
    """
    logg.info('Calculating cell cycle scores', r=True)

    adata = adata.copy() if copy else adata

    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    add_score(adata, gene_list=s_genes, score_name='S_score', ctrl_size=ctrl_size)

    # add g2m-score
    add_score(adata, gene_list=g2m_genes, score_name='G2M_score', ctrl_size=ctrl_size)

    scores = adata.obs[['S_score', 'G2M_score']]

    # default phase is S
    phase = pd.Series('S', index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores.G2M_score > scores.S_score] = 'G2M'

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = 'G1'

    adata.obs['Phase'] = phase

    return adata if copy else None

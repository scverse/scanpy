"""Calculate scores based on the expression of gene lists.
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from .. import logging as logg
from ..logging import _settings_verbosity_greater_or_equal_than


def score_genes(
        adata,
        gene_list,
        ctrl_size=50,
        gene_pool=None,
        n_bins=25,
        score_name='score',
        random_state=0,
        copy=False,
        use_raw=False):  # we use the scikit-learn convention of calling the seed "random_state"
    """Score a set of genes [Satija15]_.

    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    gene_list : iterable
        The list of gene names used for score calculation.
    ctrl_size : `int`, optional (default: 50)
        Number of reference genes to be sampled. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool : `list` or `None`, optional (default: `None`)
        Genes for sampling the reference set. Default is all genes.
    n_bins : `int`, optional (default: 25)
        Number of expression level bins for sampling.
    score_name : `str`, optional (default: `'score'`)
        Name of the field to be added in `.obs`.
    random_state : `int`, optional (default: 0)
        The random seed for sampling.
    copy : `bool`, optional (default: `False`)
        Copy `adata` or modify it inplace.
    use_raw : `bool`, optional (default: `False`)
        Use `raw` attribute of `adata` if present.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with an additional field
    `score_name`.

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    logg.info('computing score \'{}\''.format(score_name), r=True)
    adata = adata.copy() if copy else adata

    if random_state:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            logg.warn('gene: {} is not in adata.var_names and will be ignored'.format(gene))
    gene_list = set(gene_list_in_var[:])

    if not gene_pool:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    _adata = adata.raw if use_raw else adata
    # TODO: this densifies the whole data matrix for `gene_pool`
    if issparse(_adata.X):
        obs_avg = pd.Series(
            np.nanmean(
                _adata[:, gene_pool].X.toarray(), axis=0), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata[:, gene_pool].X, axis=0), index=gene_pool)  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)] # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))  # uses full r_genes if ctrl_size > len(r_genes)

    # To index, we need a list - indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)


    X_list = _adata[:, gene_list].X
    if issparse(X_list): X_list = X_list.toarray()
    X_control = _adata[:, control_genes].X
    if issparse(X_control): X_control = X_control.toarray()
    X_control = np.nanmean(X_control, axis=1)

    if len(gene_list) == 0:
        # We shouldn't even get here, but just in case
        logg.hint(
            'could not add \n'
            '    \'{}\', score of gene set (adata.obs)'.format(score_name))
        return adata if copy else None
    elif len(gene_list) == 1:
        score = _adata[:, gene_list].X - X_control
    else:
        score = np.nanmean(X_list, axis=1) - X_control

    adata.obs[score_name] = pd.Series(np.array(score).ravel(), index=adata.obs_names)

    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    logg.hint('added\n'
              '    \'{}\', score of gene set (adata.obs)'.format(score_name))
    return adata if copy else None


def score_genes_cell_cycle(
        adata,
        s_genes,
        g2m_genes,
        copy=False,
        **kwargs):
    """Score cell cycle genes [Satija15]_.

    Given two lists of genes associated to S phase and G2M phase, calculates
    scores and assigns a cell cycle phase (G1, S or G2M). See
    :func:`~scanpy.api.score_genes` for more explanation.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    s_genes : `list`
        List of genes associated with S phase.
    g2m_genes : `list`
        List of genes associated with G2M phase.
    copy : `bool`, optional (default: `False`)
        Copy `adata` or modify it inplace.
    **kwargs : optional keyword arguments
        Are passed to :func:`~scanpy.api.score_genes`. `ctrl_size` is not
        possible, as it's set as `min(len(s_genes), len(g2m_genes))`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **S_score** : `adata.obs`, dtype `object`
        The score for S phase for each cell.
    **G2M_score** : `adata.obs`, dtype `object`
        The score for G2M phase for each cell.
    **phase** : `adata.obs`, dtype `object`
        The cell cycle phase (`S`, `G2M` or `G1`) for each cell.

    See also
    --------
    score_genes

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    logg.info('calculating cell cycle phase')

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    score_genes(adata, gene_list=s_genes, score_name='S_score', ctrl_size=ctrl_size, **kwargs)
    # add g2m-score
    score_genes(adata, gene_list=g2m_genes, score_name='G2M_score', ctrl_size=ctrl_size, **kwargs)
    scores = adata.obs[['S_score', 'G2M_score']]

    # default phase is S
    phase = pd.Series('S', index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores.G2M_score > scores.S_score] = 'G2M'

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = 'G1'

    adata.obs['phase'] = phase
    logg.hint('    \'phase\', cell cycle phase (adata.obs)')
    return adata if copy else None

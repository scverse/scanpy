"""Rank genes according to differential expression.
"""

import numpy as np
import pandas as pd
from math import sqrt, floor
from scipy.sparse import issparse

from .. import utils
from .. import settings
from .. import logging as logg
from ..preprocessing import simple


def rank_genes_groups(
        adata,
        groupby,
        use_raw=True,
        groups='all',
        reference='rest',
        n_genes=100,
        only_positive=True,
        key_added=None,
        copy=False,
        method='t-test_overestim_var',
        **kwds):
    """Rank genes for characterizing groups.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    groupby : `str`
        The key of the observations grouping to consider.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    groups : `str`, `list`, optional (default: `'all'`)
        Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which comparison shall
        be restricted. If not passed, a ranking will be generated for all
        groups.
    reference : `str`, optional (default: `'rest'`)
        If `'rest'`, compare each group to the union of the rest of the group.  If
        a group identifier, compare with respect to this group.
    n_genes : `int`, optional (default: 100)
        The number of genes that appear in the returned tables.
    method : {'logreg', 't-test', 'wilcoxon', 't-test_overestim_var'}, optional (default: 't-test_overestim_var')
        If 't-test', uses t-test, if 'wilcoxon', uses Wilcoxon-Rank-Sum. If
        't-test_overestim_var', overestimates variance of each group. If
        'logreg' uses logistic regression, see [Ntranos18]_, `here
        <https://github.com/theislab/scanpy/issues/95>`_ and `here
        <http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`_, for
        why this is meaningful.
    only_positive : bool, optional (default: `True`)
        Only consider positive differences.
    **kwds : keyword parameters
        Are passed to test methods. Currently this affects only parameters that
        are passed to `sklearn.linear_model.LogisticRegression
        <http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html>`_.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).

    Returns
    -------
    Updates `adata` with the following fields.
    names : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    scores : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the score for each
        gene for each group. Ordered according to scores.
    """
    logg.info('ranking genes', r=True)
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    # for clarity, rename variable
    groups_order = groups
    if isinstance(groups_order, list) and isinstance(groups_order[0], int):
        groups_order = [str(n) for n in groups_order]
    if reference != 'rest' and reference not in set(groups_order):
        groups_order += [reference]
    if (reference != 'rest'
        and reference not in set(adata.obs[groupby].cat.categories)):
        raise ValueError('reference = {} needs to be one of groupby = {}.'
                         .format(reference,
                                 adata.obs[groupby].cat.categories.tolist()))
    groups_order, groups_masks = utils.select_groups(
        adata, groups_order, groupby)
    
    if key_added is None:
        key_added = 'rank_genes_groups'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = {
        'groupby': groupby,
        'reference': reference,
        'method': method,
        'use_raw': use_raw,
    }

    # adata_comp mocks an AnnData object if use_raw is True
    # otherwise it's just the AnnData object
    adata_comp = adata
    if adata.raw is not None and use_raw:
        adata_comp = adata.raw
    X = adata_comp.X

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    if n_genes_user > X.shape[1]:
        n_genes_user = X.shape[1]
    # in the following, n_genes is simply another name for the total number of genes
    n_genes = X.shape[1]

    rankings_gene_scores = []
    rankings_gene_names = []
    n_groups = groups_masks.shape[0]
    ns = np.zeros(n_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        ns[imask] = np.where(mask)[0].size
    logg.msg('consider \'{}\' groups:'.format(groupby), groups_order, v=4)
    logg.msg('with sizes:', ns, v=4)
    if reference != 'rest':
        ireference = np.where(groups_order == reference)[0][0]
    reference_indices = np.arange(adata_comp.n_vars, dtype=int)

    avail_methods = {'t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'}
    if method not in avail_methods:
        raise ValueError('Method must be one of {}.'.format(avail_methods))

    # TODO: all of this is probably going to be removed
    if method is 't-test_correction_factors':
        if correction_factors is None:
            raise ValueError('For this test type, you need to enter correction factors manually.')
        if len(correction_factors) != 2:
            raise ValueError('We need exactly 2 correction factors, accessible via correction_factors[i], i=0,1')
        if correction_factors[0] < 0 or correction_factors[1] < 0:
            raise ValueError('Correction factors need to be positive numbers!')

    if method in {'t-test', 't-test_overestim_var',
                  't-test_double_overestim_var', 't-test_correction_factors'}:
        # loop over all masks and compute means, variances and sample numbers
        means = np.zeros((n_groups, n_genes))
        vars = np.zeros((n_groups, n_genes))
        for imask, mask in enumerate(groups_masks):
            means[imask], vars[imask] = simple._get_mean_var(X[mask])
        # test each either against the union of all other groups or against a
        # specific group
        for igroup in range(n_groups):
            if reference == 'rest':
                mask_rest = ~groups_masks[igroup]
            else:
                if igroup == ireference: continue
                else: mask_rest = groups_masks[ireference]
            mean_rest, var_rest = simple._get_mean_var(X[mask_rest])
            if method == 't-test':
                ns_rest = np.where(mask_rest)[0].size
            elif method in {'t-test_overestim_var', 't-test_double_overestim_var'}:
                ns_rest = ns[igroup]  # hack for overestimating the variance
            elif method == 't-test_correction_factors':
                # The tendency is as follows: For the comparison group (rest), overesimate variance --> smaller ns_rest
                ns_rest = np.where(mask_rest)[0].size/correction_factors[1]

            if method in {'t-test', 't-test_overestim_var'}:
                ns_group = ns[igroup]
            elif method == 't-test_correction_factors':
                # We underestimate group variance by increasing denominator, i.e. ns_group
                ns_group = ns[igroup] * correction_factors[0]
            else:
                # We do the opposite of t-test_overestim_var
                ns_group = np.where(mask_rest)[0].size

            denominator = np.sqrt(vars[igroup]/ns_group + var_rest/ns_rest)
            denominator[np.flatnonzero(denominator == 0)] = np.nan
            scores = (means[igroup] - mean_rest) / denominator
            scores[np.isnan(scores)] = 0
            scores = scores if only_positive else np.abs(scores)
            partition = np.argpartition(scores, -n_genes_user)[-n_genes_user:]
            partial_indices = np.argsort(scores[partition])[::-1]
            global_indices = reference_indices[partition][partial_indices]
            rankings_gene_scores.append(scores[global_indices])
            rankings_gene_names.append(adata_comp.var_names[global_indices])
    elif method == 'logreg':
        from sklearn.linear_model import LogisticRegression
        if reference != 'rest':
            raise ValueError('\'logreg\' is only implemented for `reference==\'rest\'`.')
        clf = LogisticRegression(**kwds)
        clf.fit(X, adata.obs[groupby])
        scores_all = clf.coef_
        for igroup, group in enumerate(groups_order):
            scores = scores_all[igroup]
            partition = np.argpartition(scores, -n_genes_user)[-n_genes_user:]
            partial_indices = np.argsort(scores[partition])[::-1]
            global_indices = reference_indices[partition][partial_indices]
            rankings_gene_scores.append(scores[global_indices])
            rankings_gene_names.append(adata_comp.var_names[global_indices])
    elif method == 'wilcoxon':
        # Wilcoxon-rank-sum test is usually more powerful in detecting marker genes
        # Limit maximal RAM that is required by the calculation. Currently set fixed to roughly 100 MByte
        CONST_MAX_SIZE = 10000000
        ns_rest = np.zeros(n_groups, dtype=int)
        # initialize space for z-scores
        scores = np.zeros(n_genes)
        # First loop: Loop over all genes
        if reference != 'rest':
            for imask, mask in enumerate(groups_masks):
                if imask == ireference: continue
                else: mask_rest = groups_masks[ireference]
                ns_rest[imask] = np.where(mask_rest)[0].size
                if ns_rest[imask] <= 25 or ns[imask] <= 25:
                    logg.hint('Few observations in a group for '
                              'normal approximation (<=25). Lower test accuracy.')
                n_active = ns[imask]
                m_active = ns_rest[imask]
                # Now calculate gene expression ranking in chunkes:
                chunk = []
                # Calculate chunk frames
                n_genes_max_chunk = floor(CONST_MAX_SIZE / (n_active + m_active))
                if n_genes_max_chunk < n_genes - 1:
                    chunk_index = n_genes_max_chunk
                    while chunk_index < n_genes - 1:
                        chunk.append(chunk_index)
                        chunk_index = chunk_index + n_genes_max_chunk
                    chunk.append(n_genes - 1)
                else:
                    chunk.append(n_genes - 1)
                left = 0
                # Calculate rank sums for each chunk for the current mask
                for chunk_index, right in enumerate(chunk):
                    # Check if issparse is true: AnnData objects are currently sparse.csr or ndarray.
                    if issparse(X):
                        df1 = pd.DataFrame(data=X[mask, left:right].todense())
                        df2 = pd.DataFrame(data=X[mask_rest, left:right].todense(),
                                           index=np.arange(start=n_active, stop=n_active + m_active))
                    else:
                        df1 = pd.DataFrame(data=X[mask, left:right])
                        df2 = pd.DataFrame(data=X[mask_rest, left:right],
                                           index=np.arange(start=n_active, stop=n_active + m_active))
                    df1 = df1.append(df2)
                    ranks = df1.rank()
                    # sum up adjusted_ranks to calculate W_m,n
                    scores[left:right] = np.sum(ranks.loc[0:n_active, :])
                    left = right + 1
                scores = (scores - (n_active * (n_active + m_active + 1) / 2)) / sqrt(
                    (n_active * m_active * (n_active + m_active + 1) / 12))
                scores = scores if only_positive else np.abs(scores)
                scores[np.isnan(scores)] = 0
                partition = np.argpartition(scores, -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(scores[partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_scores.append(scores[global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])
        # If no reference group exists, ranking needs only to be done once (full mask)
        else:
            scores = np.zeros((n_groups, n_genes))
            chunk = []
            n_cells = X.shape[0]
            n_genes_max_chunk = floor(CONST_MAX_SIZE / n_cells)
            if n_genes_max_chunk < n_genes - 1:
                chunk_index = n_genes_max_chunk
                while chunk_index < n_genes - 1:
                    chunk.append(chunk_index)
                    chunk_index = chunk_index + n_genes_max_chunk
                chunk.append(n_genes - 1)
            else:
                chunk.append(n_genes - 1)
            left = 0
            for chunk_index, right in enumerate(chunk):
                # Check if issparse is true
                if issparse(X):
                    df1 = pd.DataFrame(data=X[:, left:right].todense())
                else:
                    df1 = pd.DataFrame(data=X[:, left:right])
                ranks = df1.rank()
                # sum up adjusted_ranks to calculate W_m,n
                for imask, mask in enumerate(groups_masks):
                    scores[imask, left:right] = np.sum(ranks.loc[mask, :])
                left = right + 1

            for imask, mask in enumerate(groups_masks):
                scores[imask, :] = (scores[imask, :] - (ns[imask] * (n_cells + 1) / 2)) / sqrt(
                    (ns[imask] * (n_cells - ns[imask]) * (n_cells + 1) / 12))
                scores = scores if only_positive else np.abs(scores)
                scores[np.isnan(scores)] = 0
                partition = np.argpartition(scores[imask, :], -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(scores[imask, partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_scores.append(scores[imask, global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])

    groups_order_save = [str(g) for g in groups_order]
    if reference != 'rest':
        groups_order_save = [g for g in groups_order if g != reference]
    adata.uns[key_added]['scores'] = np.rec.fromarrays(
        [n for n in rankings_gene_scores],
        dtype=[(rn, 'float32') for rn in groups_order_save])
    adata.uns[key_added]['names'] = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, 'U50') for rn in groups_order_save])
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns[\'{}\']`\n'
        '    \'names\', sorted np.recarray to be indexed by group ids\n'
        '    \'scores\', sorted np.recarray to be indexed by group ids'
        .format(key_added))
    return adata if copy else None

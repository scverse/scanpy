# Author: Alex Wolf (http://falexwolf.de)
#         T. Callies
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
        group_by,
        use_raw=True,
        groups='all',
        reference='rest',
        n_genes=100,
        compute_distribution=False,
        only_positive=True,
        copy=False,
        test_type='t-test_overestim_var',
        correction_factors=None):
    """Rank genes according to differential expression [Wolf17]_.

    Rank genes by differential expression. By default, a t-test-like ranking is
    used, in which means are normalized with variances.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    group_by : `str`
        The key of the sample grouping to consider.
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
    test_type : {'t-test_overestim_var', 't-test', 'wilcoxon', , 't-test_double_overestim_var',
                   't-test_correction_factors'}, optional (default: 't-test_overestim_var')
        If 't-test', use t-test to calculate test statistics. If 'wilcoxon', use
        Wilcoxon-Rank-Sum to calculate test statistic. If
        't-test_overestim_var', overestimate variance.
        't-test_double_overestim_var', additionally, underestimate variance of the rest
        't-test_correction_factors', define correction factors manually
    only_positive : bool, optional (default: `True`)
        Only consider positive differences.
    correction_factors: [a,b], optional (default: None)
        Only for the test-type 't-test_correction_factors'. Then, a determines correction factor for group variance,
        b determines correction factor for variance of the comparison group
    Returns
    -------
    rank_genes_groups_gene_scores : structured `np.ndarray` (adata.uns)
        Structured array to be indexed by group id of shape storing the zscore
        for each gene for each group.
    rank_genes_groups_gene_names : structured `np.ndarray` (adata.uns)
        Structured array to be indexed by group id for storing the gene names.
    """
    logg.info('rank differentially expressed genes', r=True)
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    if compute_distribution:
        logg.warn('`compute_distribution` is deprecated, as it requires storing'
                  'a shifted and rescaled disribution for each gene'
                  'You can now run `sc.pl.rank_genes_groups_violin` without it, '
                  'which will show the original distribution of the gene.')
    # for clarity, rename variable
    groups_order = groups
    if isinstance(groups_order, list) and isinstance(groups_order[0], int):
        groups_order = [str(n) for n in groups_order]
    if reference != 'rest' and reference not in set(groups_order):
        groups_order += [reference]
    if (reference != 'rest'
        and reference not in set(adata.obs[group_by].cat.categories)):
        raise ValueError('reference = {} needs to be one of group_by = {}.'
                         .format(reference,
                                 adata.obs[group_by].cat.categories.tolist()))
    groups_order, groups_masks = utils.select_groups(
        adata, groups_order, group_by)
    adata.uns['rank_genes_groups_params'] = np.array(
        (group_by, reference, test_type, use_raw),
        dtype=[('group_by', 'U50'), ('reference', 'U50'), ('test_type', 'U50'), ('use_raw', np.bool_)])

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

    rankings_gene_zscores = []
    rankings_gene_names = []
    n_groups = groups_masks.shape[0]
    ns = np.zeros(n_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        ns[imask] = np.where(mask)[0].size
    logg.info('    consider \'{}\':'.format(group_by), groups_order,
              'with sample numbers', ns)
    if reference != 'rest':
        ireference = np.where(groups_order == reference)[0][0]
    reference_indices = np.arange(adata_comp.n_vars, dtype=int)

    avail_tests = {'t-test', 't-test_overestim_var', 'wilcoxon', 't-test_double_overestim_var',
                   't-test_correction_factors'}
    if test_type not in avail_tests:
        raise ValueError('test_type should be one of {}.'
                         '"t-test_overestim_var" is being used as default.'
                         .format(avail_tests))

    if test_type is 't-test_correction_factors':
        if correction_factors is None:
            raise ValueError('For this test type, you need to enter correction factors manually.')
        if len(correction_factors) != 2:
            raise ValueError('We need exactly 2 correction factors, accessible via correction_factors[i], i=0,1')
        if correction_factors[0]<0 or correction_factors[1]<0:
            raise ValueError('Correction factors need to be positive numbers!')

    if test_type in {'t-test', 't-test_overestim_var', 't-test_double_overestim_var',
                   't-test_correction_factors'}:
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
            if test_type == 't-test':
                ns_rest = np.where(mask_rest)[0].size
            elif test_type == 't-test_correction_factors':
                # The tendency is as follows: For the comparison group (rest), overesimate variance --> smaller ns_rest
                ns_rest = np.where(mask_rest)[0].size/correction_factors[1]
            else:  # hack for overestimating the variance
                ns_rest = ns[igroup]

            if test_type in {'t-test', 't-test_overestim_var'}:
                ns_group=ns[igroup]
            elif test_type == 't-test_correction_factors':
                # We underestimate group variance by increasing denominator, i.e. ns_group
                ns_group=ns[igroup]*correction_factors[0]
            else :
                # We do the opposite of t-test_overestim_var
                ns_group=np.where(mask_rest)[0].size

            denominator = np.sqrt(vars[igroup]/ns_group + var_rest/ns_rest)
            denominator[np.flatnonzero(denominator == 0)] = np.nan
            zscores = (means[igroup] - mean_rest) / denominator
            zscores[np.isnan(zscores)] = 0
            zscores = zscores if only_positive else np.abs(zscores)
            partition = np.argpartition(zscores, -n_genes_user)[-n_genes_user:]
            partial_indices = np.argsort(zscores[partition])[::-1]
            global_indices = reference_indices[partition][partial_indices]
            rankings_gene_zscores.append(zscores[global_indices])
            rankings_gene_names.append(adata_comp.var_names[global_indices])
            if compute_distribution:
                mask = groups_masks[igroup]
                for gene_counter in range(n_genes_user):
                    gene_idx = global_indices[gene_counter]
                    X_col = X[mask, gene_idx]
                    if issparse(X): X_col = X_col.toarray()[:, 0]
                    identifier = _build_identifier(group_by, groups_order[igroup],
                                                   gene_counter, adata_comp.var_names[gene_idx])
                    full_col = np.empty(adata.n_obs)
                    full_col[:] = np.nan
                    full_col[mask] = (X_col - mean_rest[gene_idx]) / denominator[gene_idx]
                    adata.obs[identifier] = full_col
    elif test_type == 'wilcoxon':
        # Wilcoxon-rank-sum test is usually more powerful in detecting marker genes
        # Limit maximal RAM that is required by the calculation. Currently set fixed to roughly 100 MByte
        CONST_MAX_SIZE = 10000000
        ns_rest = np.zeros(n_groups, dtype=int)
        # initialize space for z-scores
        zscores = np.zeros(n_genes)
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
                    zscores[left:right] = np.sum(ranks.loc[0:n_active, :])
                    left = right + 1
                zscores = (zscores - (n_active * (n_active + m_active + 1) / 2)) / sqrt(
                    (n_active * m_active * (n_active + m_active + 1) / 12))
                zscores = zscores if only_positive else np.abs(zscores)
                zscores[np.isnan(zscores)] = 0
                partition = np.argpartition(zscores, -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(zscores[partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_zscores.append(zscores[global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])
                if compute_distribution:
                    # Add calculation of means, var: (Unnecessary for wilcoxon if compute distribution=False)
                    mean, vars = simple._get_mean_var(X[mask])
                    mean_rest, var_rest = simple._get_mean_var(X[mask_rest])
                    denominator = np.sqrt(vars / ns[imask] + var_rest / ns_rest[imask])
                    denominator[np.flatnonzero(denominator == 0)] = np.nan
                    for gene_counter in range(n_genes_user):
                        gene_idx = global_indices[gene_counter]
                        X_col = X[mask, gene_idx]
                        if issparse(X): X_col = X_col.toarray()[:, 0]
                        identifier = _build_identifier(group_by, groups_order[imask],
                                                       gene_counter, adata_comp.var_names[gene_idx])
                        full_col = np.empty(adata.n_obs)
                        full_col[:] = np.nan
                        full_col[mask] = (X_col - mean_rest[gene_idx]) / denominator[gene_idx]
                        adata.obs[identifier] = full_col

        # If no reference group exists, ranking needs only to be done once (full mask)
        else:
            zscores = np.zeros((n_groups, n_genes))
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
                    zscores[imask, left:right] = np.sum(ranks.loc[mask, :])
                left = right + 1

            for imask, mask in enumerate(groups_masks):
                zscores[imask, :] = (zscores[imask, :] - (ns[imask] * (n_cells + 1) / 2)) / sqrt(
                    (ns[imask] * (n_cells - ns[imask]) * (n_cells + 1) / 12))
                zscores = zscores if only_positive else np.abs(zscores)
                zscores[np.isnan(zscores)] = 0
                partition = np.argpartition(zscores[imask, :], -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(zscores[imask, partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_zscores.append(zscores[imask, global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])
                if compute_distribution:
                    mean, vars = simple._get_mean_var(X[mask])
                    mean_rest, var_rest = simple._get_mean_var(X[~mask])
                    denominator = np.sqrt(vars / ns[imask] + var_rest / (n_cells-ns[imask]))
                    denominator[np.flatnonzero(denominator == 0)] = np.nan
                    for gene_counter in range(n_genes_user):
                        gene_idx = global_indices[gene_counter]
                        X_col = X[mask, gene_idx]
                        if issparse(X): X_col = X_col.toarray()[:, 0]
                        identifier = _build_identifier(group_by, groups_order[imask],
                                                       gene_counter, adata_comp.var_names[gene_idx])
                        full_col = np.empty(adata.n_obs)
                        full_col[:] = np.nan
                        full_col[mask] = (X_col - mean_rest[gene_idx]) / denominator[gene_idx]
                        adata.obs[identifier] = full_col

    groups_order_save = [str(g) for g in groups_order]
    if reference != 'rest':
        groups_order_save = [g for g in groups_order if g != reference]
    adata.uns['rank_genes_groups_gene_scores'] = np.rec.fromarrays(
        [n for n in rankings_gene_zscores],
        dtype=[(rn, 'float32') for rn in groups_order_save])
    adata.uns['rank_genes_groups_gene_names'] = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, 'U50') for rn in groups_order_save])
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
           '    \'rank_genes_groups_gene_names\', np.recarray to be indexed by group ids (adata.uns)\n'
           '    \'rank_genes_groups_gene_scores\', np.recarray to be indexed by group ids (adata.uns)')
    return adata if copy else None


def _build_identifier(group_by, name, gene_counter, gene_name):
    return 'rank_genes_{}_{}_{}_{}'.format(
        group_by, name, gene_counter, gene_name)

"""Rank genes according to differential expression.
"""
from math import floor
from typing import Iterable, Union, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse, vstack

from .. import _utils
from .. import logging as logg
from ..preprocessing._simple import _get_mean_var
from .._compat import Literal
from ..get import _get_obs_rep
from .._utils import check_nonnegative_integers


_Method = Optional[Literal['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']]
_CorrMethod = Literal['benjamini-hochberg', 'bonferroni']


def _select_top_n(scores, n_top):
    n_from = scores.shape[0]
    reference_indices = np.arange(n_from, dtype=int)
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]

    return global_indices


def _ranks(X, mask=None, mask_rest=None):
    CONST_MAX_SIZE = 10000000

    n_genes = X.shape[1]

    if issparse(X):
        merge = lambda tpl: vstack(tpl).toarray()
        adapt = lambda X: X.toarray()
    else:
        merge = np.vstack
        adapt = lambda X: X

    masked = mask is not None and mask_rest is not None

    if masked:
        n_cells = np.count_nonzero(mask) + np.count_nonzero(mask_rest)
        get_chunk = lambda X, left, right: merge(
            (X[mask, left:right], X[mask_rest, left:right])
        )
    else:
        n_cells = X.shape[0]
        get_chunk = lambda X, left, right: adapt(X[:, left:right])

    # Calculate chunk frames
    max_chunk = floor(CONST_MAX_SIZE / n_cells)

    for left in range(0, n_genes, max_chunk):
        right = min(left + max_chunk, n_genes)

        df = pd.DataFrame(data=get_chunk(X, left, right))
        ranks = df.rank()
        yield ranks, left, right


def _tiecorrect(ranks):
    size = np.float64(ranks.shape[0])
    if size < 2:
        return np.repeat(ranks.shape[1], 1.0)

    arr = np.sort(ranks, axis=0)
    tf = np.insert(arr[1:] != arr[:-1], (0, arr.shape[0] - 1), True, axis=0)
    idx = np.where(tf, np.arange(tf.shape[0])[:, None], 0)
    idx = np.sort(idx, axis=0)
    cnt = np.diff(idx, axis=0).astype(np.float64)

    return 1.0 - (cnt**3 - cnt).sum(axis=0) / (size**3 - size)


class _RankGenes:
    def __init__(
        self,
        adata,
        groups,
        groupby,
        reference='rest',
        use_raw=True,
        layer=None,
        comp_pts=False,
    ):

        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            self.expm1_func = lambda x: np.expm1(x * np.log(adata.uns['log1p']['base']))
        else:
            self.expm1_func = np.expm1

        self.groups_order, self.groups_masks = _utils.select_groups(
            adata, groups, groupby
        )

        # Singlet groups cause division by zero errors
        invalid_groups_selected = set(self.groups_order) & set(
            adata.obs[groupby].value_counts().loc[lambda x: x < 2].index
        )

        if len(invalid_groups_selected) > 0:
            raise ValueError(
                "Could not calculate statistics for groups {} since they only "
                "contain one sample.".format(', '.join(invalid_groups_selected))
            )

        adata_comp = adata
        if layer is not None:
            if use_raw:
                raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
            X = adata_comp.layers[layer]
        else:
            if use_raw and adata.raw is not None:
                adata_comp = adata.raw
            X = adata_comp.X

        # for correct getnnz calculation
        if issparse(X):
            X.eliminate_zeros()

        self.X = X
        self.var_names = adata_comp.var_names

        self.ireference = None
        if reference != 'rest':
            self.ireference = np.where(self.groups_order == reference)[0][0]

        self.means = None
        self.vars = None

        self.means_rest = None
        self.vars_rest = None

        self.comp_pts = comp_pts
        self.pts = None
        self.pts_rest = None

        self.stats = None

        # for logreg only
        self.grouping_mask = adata.obs[groupby].isin(self.groups_order)
        self.grouping = adata.obs.loc[self.grouping_mask, groupby]

    def _basic_stats(self):
        n_genes = self.X.shape[1]
        n_groups = self.groups_masks.shape[0]

        self.means = np.zeros((n_groups, n_genes))
        self.vars = np.zeros((n_groups, n_genes))
        self.pts = np.zeros((n_groups, n_genes)) if self.comp_pts else None

        if self.ireference is None:
            self.means_rest = np.zeros((n_groups, n_genes))
            self.vars_rest = np.zeros((n_groups, n_genes))
            self.pts_rest = np.zeros((n_groups, n_genes)) if self.comp_pts else None
        else:
            mask_rest = self.groups_masks[self.ireference]
            X_rest = self.X[mask_rest]
            self.means[self.ireference], self.vars[self.ireference] = _get_mean_var(
                X_rest
            )
            # deleting the next line causes a memory leak for some reason
            del X_rest

        if issparse(self.X):
            get_nonzeros = lambda X: X.getnnz(axis=0)
        else:
            get_nonzeros = lambda X: np.count_nonzero(X, axis=0)

        for imask, mask in enumerate(self.groups_masks):
            X_mask = self.X[mask]

            if self.comp_pts:
                self.pts[imask] = get_nonzeros(X_mask) / X_mask.shape[0]

            if self.ireference is not None and imask == self.ireference:
                continue

            self.means[imask], self.vars[imask] = _get_mean_var(X_mask)

            if self.ireference is None:
                mask_rest = ~mask
                X_rest = self.X[mask_rest]
                self.means_rest[imask], self.vars_rest[imask] = _get_mean_var(X_rest)
                # this can be costly for sparse data
                if self.comp_pts:
                    self.pts_rest[imask] = get_nonzeros(X_rest) / X_rest.shape[0]
                # deleting the next line causes a memory leak for some reason
                del X_rest

    def t_test(self, method):
        from scipy import stats

        self._basic_stats()

        for group_index, mask in enumerate(self.groups_masks):
            if self.ireference is not None and group_index == self.ireference:
                continue

            mean_group = self.means[group_index]
            var_group = self.vars[group_index]
            ns_group = np.count_nonzero(mask)

            if self.ireference is not None:
                mean_rest = self.means[self.ireference]
                var_rest = self.vars[self.ireference]
                ns_other = np.count_nonzero(self.groups_masks[self.ireference])
            else:
                mean_rest = self.means_rest[group_index]
                var_rest = self.vars_rest[group_index]
                ns_other = self.X.shape[0] - ns_group

            if method == 't-test':
                ns_rest = ns_other
            elif method == 't-test_overestim_var':
                # hack for overestimating the variance for small groups
                ns_rest = ns_group
            else:
                raise ValueError('Method does not exist.')

            # TODO: Come up with better solution. Mask unexpressed genes?
            # See https://github.com/scipy/scipy/issues/10269
            with np.errstate(invalid="ignore"):
                scores, pvals = stats.ttest_ind_from_stats(
                    mean1=mean_group,
                    std1=np.sqrt(var_group),
                    nobs1=ns_group,
                    mean2=mean_rest,
                    std2=np.sqrt(var_rest),
                    nobs2=ns_rest,
                    equal_var=False,  # Welch's
                )

            # I think it's only nan when means are the same and vars are 0
            scores[np.isnan(scores)] = 0
            # This also has to happen for Benjamini Hochberg
            pvals[np.isnan(pvals)] = 1

            yield group_index, scores, pvals

    def wilcoxon(self, tie_correct):
        from scipy import stats

        self._basic_stats()

        n_genes = self.X.shape[1]
        # First loop: Loop over all genes
        if self.ireference is not None:
            # initialize space for z-scores
            scores = np.zeros(n_genes)
            # initialize space for tie correction coefficients
            if tie_correct:
                T = np.zeros(n_genes)
            else:
                T = 1

            for group_index, mask in enumerate(self.groups_masks):
                if group_index == self.ireference:
                    continue

                mask_rest = self.groups_masks[self.ireference]

                n_active = np.count_nonzero(mask)
                m_active = np.count_nonzero(mask_rest)

                if n_active <= 25 or m_active <= 25:
                    logg.hint(
                        'Few observations in a group for '
                        'normal approximation (<=25). Lower test accuracy.'
                    )

                # Calculate rank sums for each chunk for the current mask
                for ranks, left, right in _ranks(self.X, mask, mask_rest):
                    scores[left:right] = np.sum(ranks.iloc[0:n_active, :])
                    if tie_correct:
                        T[left:right] = _tiecorrect(ranks)

                std_dev = np.sqrt(
                    T * n_active * m_active * (n_active + m_active + 1) / 12.0
                )

                scores = (
                    scores - (n_active * ((n_active + m_active + 1) / 2.0))
                ) / std_dev
                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores))

                yield group_index, scores, pvals
        # If no reference group exists,
        # ranking needs only to be done once (full mask)
        else:
            n_groups = self.groups_masks.shape[0]
            scores = np.zeros((n_groups, n_genes))
            n_cells = self.X.shape[0]

            if tie_correct:
                T = np.zeros((n_groups, n_genes))

            for ranks, left, right in _ranks(self.X):
                # sum up adjusted_ranks to calculate W_m,n
                for imask, mask in enumerate(self.groups_masks):
                    scores[imask, left:right] = np.sum(ranks.iloc[mask, :])
                    if tie_correct:
                        T[imask, left:right] = _tiecorrect(ranks)

            for group_index, mask in enumerate(self.groups_masks):
                n_active = np.count_nonzero(mask)

                if tie_correct:
                    T_i = T[group_index]
                else:
                    T_i = 1

                std_dev = np.sqrt(
                    T_i * n_active * (n_cells - n_active) * (n_cells + 1) / 12.0
                )

                scores[group_index, :] = (
                    scores[group_index, :] - (n_active * (n_cells + 1) / 2.0)
                ) / std_dev
                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores[group_index, :]))

                yield group_index, scores[group_index], pvals

    def logreg(self, **kwds):
        # if reference is not set, then the groups listed will be compared to the rest
        # if reference is set, then the groups listed will be compared only to the other groups listed
        from sklearn.linear_model import LogisticRegression

        # Indexing with a series causes issues, possibly segfault
        X = self.X[self.grouping_mask.values, :]

        if len(self.groups_order) == 1:
            raise ValueError('Cannot perform logistic regression on a single cluster.')

        clf = LogisticRegression(**kwds)
        clf.fit(X, self.grouping.cat.codes)
        scores_all = clf.coef_
        for igroup, _ in enumerate(self.groups_order):
            if len(self.groups_order) <= 2:  # binary logistic regression
                scores = scores_all[0]
            else:
                scores = scores_all[igroup]

            yield igroup, scores, None

            if len(self.groups_order) <= 2:
                break

    def compute_statistics(
        self,
        method,
        corr_method='benjamini-hochberg',
        n_genes_user=None,
        rankby_abs=False,
        tie_correct=False,
        **kwds,
    ):

        if method in {'t-test', 't-test_overestim_var'}:
            generate_test_results = self.t_test(method)
        elif method == 'wilcoxon':
            generate_test_results = self.wilcoxon(tie_correct)
        elif method == 'logreg':
            generate_test_results = self.logreg(**kwds)

        self.stats = None

        n_genes = self.X.shape[1]

        for group_index, scores, pvals in generate_test_results:
            group_name = str(self.groups_order[group_index])

            if n_genes_user is not None:
                scores_sort = np.abs(scores) if rankby_abs else scores
                global_indices = _select_top_n(scores_sort, n_genes_user)
                first_col = 'names'
            else:
                global_indices = slice(None)
                first_col = 'scores'

            if self.stats is None:
                idx = pd.MultiIndex.from_tuples([(group_name, first_col)])
                self.stats = pd.DataFrame(columns=idx)

            if n_genes_user is not None:
                self.stats[group_name, 'names'] = self.var_names[global_indices]

            self.stats[group_name, 'scores'] = scores[global_indices]

            if pvals is not None:
                self.stats[group_name, 'pvals'] = pvals[global_indices]
                if corr_method == 'benjamini-hochberg':
                    from statsmodels.stats.multitest import multipletests

                    pvals[np.isnan(pvals)] = 1
                    _, pvals_adj, _, _ = multipletests(
                        pvals, alpha=0.05, method='fdr_bh'
                    )
                elif corr_method == 'bonferroni':
                    pvals_adj = np.minimum(pvals * n_genes, 1.0)
                self.stats[group_name, 'pvals_adj'] = pvals_adj[global_indices]

            if self.means is not None:
                mean_group = self.means[group_index]
                if self.ireference is None:
                    mean_rest = self.means_rest[group_index]
                else:
                    mean_rest = self.means[self.ireference]
                foldchanges = (self.expm1_func(mean_group) + 1e-9) / (
                    self.expm1_func(mean_rest) + 1e-9
                )  # add small value to remove 0's
                self.stats[group_name, 'logfoldchanges'] = np.log2(
                    foldchanges[global_indices]
                )

        if n_genes_user is None:
            self.stats.index = self.var_names


# TODO: Make arguments after groupby keyword only
def rank_genes_groups(
    adata: AnnData,
    groupby: str,
    use_raw: Optional[bool] = None,
    groups: Union[Literal['all'], Iterable[str]] = 'all',
    reference: str = 'rest',
    n_genes: Optional[int] = None,
    rankby_abs: bool = False,
    pts: bool = False,
    key_added: Optional[str] = None,
    copy: bool = False,
    method: _Method = None,
    corr_method: _CorrMethod = 'benjamini-hochberg',
    tie_correct: bool = False,
    layer: Optional[str] = None,
    **kwds,
) -> Optional[AnnData]:
    """\
    Rank genes for characterizing groups.

    Expects logarithmized data.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    groups
        Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
        shall be restricted, or `'all'` (default), for all groups.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    n_genes
        The number of genes that appear in the returned tables.
        Defaults to all genes.
    method
        The default method is `'t-test'`,
        `'t-test_overestim_var'` overestimates variance of each group,
        `'wilcoxon'` uses Wilcoxon rank-sum,
        `'logreg'` uses logistic regression. See [Ntranos18]_,
        `here <https://github.com/scverse/scanpy/issues/95>`__ and `here
        <http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,
        for why this is meaningful.
    corr_method
        p-value correction method.
        Used only for `'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'`.
    tie_correct
        Use tie correction for `'wilcoxon'` scores.
        Used only for `'wilcoxon'`.
    rankby_abs
        Rank genes by the absolute value of the score, not by the
        score. The returned scores are never the absolute values.
    pts
        Compute the fraction of cells expressing the genes.
    key_added
        The key in `adata.uns` information is saved to.
    **kwds
        Are passed to test methods. Currently this affects only parameters that
        are passed to :class:`sklearn.linear_model.LogisticRegression`.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).

    Returns
    -------
    **names** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    **scores** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the z-score
        underlying the computation of a p-value for each gene for each
        group. Ordered according to scores.
    **logfoldchanges** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the log2
        fold change for each gene for each group. Ordered according to
        scores. Only provided if method is 't-test' like.
        Note: this is an approximation calculated from mean-log values.
    **pvals** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        p-values.
    **pvals_adj** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Corrected p-values.
    **pts** : `pandas.DataFrame` (`.uns['rank_genes_groups']`)
        Fraction of cells expressing the genes for each group.
    **pts_rest** : `pandas.DataFrame` (`.uns['rank_genes_groups']`)
        Only if `reference` is set to `'rest'`.
        Fraction of cells from the union of the rest of each group
        expressing the genes.

    Notes
    -----
    There are slight inconsistencies depending on whether sparse
    or dense data are passed. See `here <https://github.com/scverse/scanpy/blob/master/scanpy/tests/test_rank_genes_groups.py>`__.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels', method='wilcoxon')
    >>> # to visualize the results
    >>> sc.pl.rank_genes_groups(adata)
    """
    if use_raw is None:
        use_raw = adata.raw is not None
    elif use_raw is True and adata.raw is None:
        raise ValueError("Received `use_raw=True`, but `adata.raw` is empty.")

    if method is None:
        logg.warning(
            "Default of the method has been changed to 't-test' from 't-test_overestim_var'"
        )
        method = 't-test'

    if 'only_positive' in kwds:
        rankby_abs = not kwds.pop('only_positive')  # backwards compat

    start = logg.info('ranking genes')
    avail_methods = {'t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'}
    if method not in avail_methods:
        raise ValueError(f'Method must be one of {avail_methods}.')

    avail_corr = {'benjamini-hochberg', 'bonferroni'}
    if corr_method not in avail_corr:
        raise ValueError(f'Correction method must be one of {avail_corr}.')

    adata = adata.copy() if copy else adata
    _utils.sanitize_anndata(adata)
    # for clarity, rename variable
    if groups == 'all':
        groups_order = 'all'
    elif isinstance(groups, (str, int)):
        raise ValueError('Specify a sequence of groups')
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if reference != 'rest' and reference not in set(groups_order):
            groups_order += [reference]
    if reference != 'rest' and reference not in adata.obs[groupby].cat.categories:
        cats = adata.obs[groupby].cat.categories.tolist()
        raise ValueError(
            f'reference = {reference} needs to be one of groupby = {cats}.'
        )

    if key_added is None:
        key_added = 'rank_genes_groups'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = dict(
        groupby=groupby,
        reference=reference,
        method=method,
        use_raw=use_raw,
        layer=layer,
        corr_method=corr_method,
    )

    test_obj = _RankGenes(adata, groups_order, groupby, reference, use_raw, layer, pts)

    if check_nonnegative_integers(test_obj.X) and method != 'logreg':
        logg.warning(
            "It seems you use rank_genes_groups on the raw count data. "
            "Please logarithmize your data before calling rank_genes_groups."
        )

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    # defaults to all genes
    if n_genes_user is None or n_genes_user > test_obj.X.shape[1]:
        n_genes_user = test_obj.X.shape[1]

    logg.debug(f'consider {groupby!r} groups:')
    logg.debug(f'with sizes: {np.count_nonzero(test_obj.groups_masks, axis=1)}')

    test_obj.compute_statistics(
        method, corr_method, n_genes_user, rankby_abs, tie_correct, **kwds
    )

    if test_obj.pts is not None:
        groups_names = [str(name) for name in test_obj.groups_order]
        adata.uns[key_added]['pts'] = pd.DataFrame(
            test_obj.pts.T, index=test_obj.var_names, columns=groups_names
        )
    if test_obj.pts_rest is not None:
        adata.uns[key_added]['pts_rest'] = pd.DataFrame(
            test_obj.pts_rest.T, index=test_obj.var_names, columns=groups_names
        )

    test_obj.stats.columns = test_obj.stats.columns.swaplevel()

    dtypes = {
        'names': 'O',
        'scores': 'float32',
        'logfoldchanges': 'float32',
        'pvals': 'float64',
        'pvals_adj': 'float64',
    }

    for col in test_obj.stats.columns.levels[0]:
        adata.uns[key_added][col] = test_obj.stats[col].to_records(
            index=False, column_dtypes=dtypes[col]
        )

    logg.info(
        '    finished',
        time=start,
        deep=(
            f'added to `.uns[{key_added!r}]`\n'
            "    'names', sorted np.recarray to be indexed by group ids\n"
            "    'scores', sorted np.recarray to be indexed by group ids\n"
            + (
                "    'logfoldchanges', sorted np.recarray to be indexed by group ids\n"
                "    'pvals', sorted np.recarray to be indexed by group ids\n"
                "    'pvals_adj', sorted np.recarray to be indexed by group ids"
                if method in {'t-test', 't-test_overestim_var', 'wilcoxon'}
                else ''
            )
        ),
    )
    return adata if copy else None


def _calc_frac(X):
    if issparse(X):
        n_nonzero = X.getnnz(axis=0)
    else:
        n_nonzero = np.count_nonzero(X, axis=0)
    return n_nonzero / X.shape[0]


def filter_rank_genes_groups(
    adata: AnnData,
    key=None,
    groupby=None,
    use_raw=None,
    key_added='rank_genes_groups_filtered',
    min_in_group_fraction=0.25,
    min_fold_change=1,
    max_out_group_fraction=0.5,
    compare_abs=False,
) -> None:
    """\
    Filters out genes based on log fold change and fraction of genes expressing the
    gene within and outside the `groupby` categories.

    See :func:`~scanpy.tl.rank_genes_groups`.

    Results are stored in `adata.uns[key_added]`
    (default: 'rank_genes_groups_filtered').

    To preserve the original structure of adata.uns['rank_genes_groups'],
    filtered genes are set to `NaN`.

    Parameters
    ----------
    adata
    key
    groupby
    use_raw
    key_added
    min_in_group_fraction
    min_fold_change
    max_out_group_fraction
    compare_abs
        If `True`, compare absolute values of log fold change with `min_fold_change`.

    Returns
    -------
    Same output as :func:`scanpy.tl.rank_genes_groups` but with filtered genes names set to
    `nan`

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels', method='wilcoxon')
    >>> sc.tl.filter_rank_genes_groups(adata, min_fold_change=3)
    >>> # visualize results
    >>> sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered')
    >>> # visualize results using dotplot
    >>> sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered')
    """
    if key is None:
        key = 'rank_genes_groups'

    if groupby is None:
        groupby = adata.uns[key]['params']['groupby']

    if use_raw is None:
        use_raw = adata.uns[key]['params']['use_raw']

    same_params = (
        adata.uns[key]['params']['groupby'] == groupby
        and adata.uns[key]['params']['reference'] == 'rest'
        and adata.uns[key]['params']['use_raw'] == use_raw
    )

    use_logfolds = same_params and 'logfoldchanges' in adata.uns[key]
    use_fraction = same_params and 'pts_rest' in adata.uns[key]

    # convert structured numpy array into DataFrame
    gene_names = pd.DataFrame(adata.uns[key]['names'])

    fraction_in_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )
    fraction_out_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )

    if use_logfolds:
        fold_change_matrix = pd.DataFrame(adata.uns[key]['logfoldchanges'])
    else:
        fold_change_matrix = pd.DataFrame(
            np.zeros(gene_names.shape),
            columns=gene_names.columns,
            index=gene_names.index,
        )

        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            expm1_func = lambda x: np.expm1(x * np.log(adata.uns['log1p']['base']))
        else:
            expm1_func = np.expm1

    logg.info(
        f"Filtering genes using: "
        f"min_in_group_fraction: {min_in_group_fraction} "
        f"min_fold_change: {min_fold_change}, "
        f"max_out_group_fraction: {max_out_group_fraction}"
    )

    for cluster in gene_names.columns:
        # iterate per column
        var_names = gene_names[cluster].values

        if not use_logfolds or not use_fraction:
            sub_X = adata.raw[:, var_names].X if use_raw else adata[:, var_names].X
            in_group = adata.obs[groupby] == cluster
            X_in = sub_X[in_group]
            X_out = sub_X[~in_group]

        if use_fraction:
            fraction_in_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]['pts'][cluster].loc[var_names].values
            )
            fraction_out_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]['pts_rest'][cluster].loc[var_names].values
            )
        else:
            fraction_in_cluster_matrix.loc[:, cluster] = _calc_frac(X_in)
            fraction_out_cluster_matrix.loc[:, cluster] = _calc_frac(X_out)

        if not use_logfolds:
            # compute mean value
            mean_in_cluster = np.ravel(X_in.mean(0))
            mean_out_cluster = np.ravel(X_out.mean(0))
            # compute fold change
            fold_change_matrix.loc[:, cluster] = np.log2(
                (expm1_func(mean_in_cluster) + 1e-9)
                / (expm1_func(mean_out_cluster) + 1e-9)
            )

    if compare_abs:
        fold_change_matrix = fold_change_matrix.abs()
    # filter original_matrix
    gene_names = gene_names[
        (fraction_in_cluster_matrix > min_in_group_fraction)
        & (fraction_out_cluster_matrix < max_out_group_fraction)
        & (fold_change_matrix > min_fold_change)
    ]
    # create new structured array using 'key_added'.
    adata.uns[key_added] = adata.uns[key].copy()
    adata.uns[key_added]['names'] = gene_names.to_records(index=False)

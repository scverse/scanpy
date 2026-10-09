from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from anndata import AnnData

from ... import logging as logg
from ..._settings import settings
from ..._utils import _numba_thread_limit
from ._kernels import _ranks, _tiecorrect

if TYPE_CHECKING:
    from collections.abc import Generator, Iterator

    from numpy.typing import NDArray

    from ._comparison import _Comparison
    from ._stats import _GroupStats


type _TestResult = tuple[int, NDArray[np.floating], NDArray[np.floating] | None]


def _t_test(
    cmp: _Comparison, stats: _GroupStats, *, overestimate_var: bool
) -> Generator[_TestResult]:
    from scipy import stats as st

    for group_index, (mask_obs, mean_group, var_group) in enumerate(
        zip(cmp.groups_masks_obs, stats.means, stats.vars, strict=True)
    ):
        if cmp.ireference is not None and group_index == cmp.ireference:
            continue

        ns_group = np.count_nonzero(mask_obs)

        if cmp.ireference is not None:
            mean_rest = stats.means[cmp.ireference]
            var_rest = stats.vars[cmp.ireference]
            ns_other = np.count_nonzero(cmp.groups_masks_obs[cmp.ireference])
        else:
            mean_rest = stats.means_rest[group_index]
            var_rest = stats.vars_rest[group_index]
            ns_other = cmp.x.shape[0] - ns_group

        # hack for overestimating the variance for small groups
        ns_rest = ns_group if overestimate_var else ns_other

        # TODO: Come up with better solution. Mask unexpressed genes?
        # See https://github.com/scipy/scipy/issues/10269
        with np.errstate(invalid="ignore"):
            scores, pvals = st.ttest_ind_from_stats(
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


def _wilcoxon(cmp: _Comparison, *, tie_correct: bool) -> Generator[_TestResult]:
    if cmp.ireference is None:
        yield from _wilcoxon_vs_rest(cmp, tie_correct=tie_correct)
    else:
        yield from _wilcoxon_vs_reference(cmp, tie_correct=tie_correct)


def _wilcoxon_vs_reference(
    cmp: _Comparison, *, tie_correct: bool
) -> Generator[_TestResult]:
    from scipy import stats

    n_genes = cmp.x.shape[1]
    # initialize space for tie correction coefficients
    tc_coef = np.zeros(n_genes) if tie_correct else 1
    mask_obs_rest = cmp.groups_masks_obs[cmp.ireference]
    m_active = np.count_nonzero(mask_obs_rest)

    for group_index, mask_obs in enumerate(cmp.groups_masks_obs):
        if group_index == cmp.ireference:
            continue

        n_active = np.count_nonzero(mask_obs)

        if n_active <= 25 or m_active <= 25:
            logg.hint(
                "Few observations in a group for "
                "normal approximation (<=25). Lower test accuracy."
            )

        scores = np.zeros(n_genes)
        # Calculate rank sums for each chunk for the current mask
        for ranks, left, right in _ranks(cmp.x, mask_obs, mask_obs_rest):
            scores[left:right] = ranks[0:n_active, :].sum(axis=0)
            if tie_correct:
                tc_coef[left:right] = _tiecorrect(ranks)

        std_dev = np.sqrt(
            tc_coef * n_active * m_active * (n_active + m_active + 1) / 12.0
        )

        scores = (scores - (n_active * ((n_active + m_active + 1) / 2.0))) / std_dev
        scores[np.isnan(scores)] = 0
        pvals = 2 * stats.distributions.norm.sf(np.abs(scores))

        yield group_index, scores, pvals


def _wilcoxon_vs_rest(cmp: _Comparison, *, tie_correct: bool) -> Generator[_TestResult]:
    from scipy import stats

    # ranking needs only to be done once (full mask)
    n_genes = cmp.x.shape[1]
    scores = np.zeros((cmp.n_groups, n_genes))
    n_cells = cmp.x.shape[0]

    if tie_correct:
        tc_coef = np.zeros((cmp.n_groups, n_genes))

    for ranks, left, right in _ranks(cmp.x):
        if tie_correct:
            tc_coef[:, left:right] = _tiecorrect(ranks)
        # sum up adjusted_ranks to calculate W_m,n
        for group_index, mask_obs in enumerate(cmp.groups_masks_obs):
            scores[group_index, left:right] = ranks[mask_obs, :].sum(axis=0)

    for group_index, mask_obs in enumerate(cmp.groups_masks_obs):
        n_active = np.count_nonzero(mask_obs)

        coef = tc_coef[group_index] if tie_correct else 1

        std_dev = np.sqrt(coef * n_active * (n_cells - n_active) * (n_cells + 1) / 12.0)

        scores[group_index, :] = (
            scores[group_index, :] - (n_active * (n_cells + 1) / 2.0)
        ) / std_dev
        scores[np.isnan(scores)] = 0
        pvals = 2 * stats.distributions.norm.sf(np.abs(scores[group_index, :]))

        yield group_index, scores[group_index], pvals


def _illico_results_to_iter(
    illico_df: pd.DataFrame,
    groups_order: NDArray,
    ireference: int | None,
) -> Iterator[_TestResult]:
    """Yield per-group ``(index, z, p)`` from illico's long-form output.

    illico's frame is group-major with features in ``var`` order, so its
    ``z_score``/``p_value`` columns reshape directly to ``(n_present_groups, n_genes)``.
    """
    reference = None if ireference is None else groups_order[ireference]
    group_names = illico_df.index.get_level_values("pert").unique()
    zscores = illico_df["z_score"].to_numpy().reshape(len(group_names), -1)
    pvals = illico_df["p_value"].to_numpy().reshape(len(group_names), -1)
    group_index = {name: i for i, name in enumerate(groups_order)}
    return (
        (group_index[name], zscores[row], pvals[row])
        for row, name in enumerate(group_names)
        if name != reference
    )


def _illico(cmp: _Comparison, *, tie_correct: bool) -> Iterator[_TestResult]:
    from illico import asymptotic_wilcoxon

    adata = AnnData(
        X=cmp.x,
        var=pd.DataFrame(index=cmp.var_names),
        obs=pd.DataFrame(
            index=pd.RangeIndex(cmp.x.shape[0]).astype("str"),
            data={"group": cmp.labels},
        ),
    )
    reference = None if cmp.ireference is None else cmp.groups_order[cmp.ireference]
    with _numba_thread_limit(settings.n_jobs) as n_threads:
        result = asymptotic_wilcoxon(
            adata,
            reference=reference,
            group_keys="group",
            is_log1p=True,
            tie_correct=tie_correct,
            use_continuity=False,
            alternative="two-sided",
            use_rust=False,
            n_threads=n_threads,
            groups=cmp.groups_order,
            return_as_scanpy=False,
        )
    return _illico_results_to_iter(result, cmp.groups_order, cmp.ireference)


def _fit_logreg(cmp: _Comparison, **kwds):
    from sklearn.linear_model import LogisticRegression

    if len(cmp.groups_order) == 1:
        msg = "Cannot perform logistic regression on a single cluster."
        raise ValueError(msg)

    # Indexing with a series causes issues, possibly segfault
    x = cmp.x[cmp.grouping_mask, :]
    clf = LogisticRegression(**kwds)
    clf.fit(x, cmp.grouping.codes)
    return clf


def _logreg(cmp: _Comparison, **kwds) -> Generator[_TestResult]:
    # if reference is not set, then the groups listed will be compared to the rest
    # if reference is set, then the groups listed will be compared only to the other groups listed
    scores_all = _fit_logreg(cmp, **kwds).coef_
    # not all codes necessarily appear in data
    existing_codes = np.unique(cmp.grouping.codes)
    for igroup, cat in enumerate(cmp.groups_order):
        if len(cmp.groups_order) <= 2:  # binary logistic regression
            scores = scores_all[0]
        else:
            # cat code is index of cat value in .categories
            cat_code: int = np.argmax(cmp.grouping.categories == cat)
            # index of scores row is index of cat code in array of existing codes
            scores_idx: int = np.argmax(existing_codes == cat_code)
            scores = scores_all[scores_idx]
        yield igroup, scores, None

        if len(cmp.groups_order) <= 2:
            break


def _logreg_signed(cmp: _Comparison, **kwds) -> Generator[_TestResult]:
    clf = _fit_logreg(cmp, **kwds)
    categories = pd.Index(cmp.grouping.categories)
    for igroup, cat in enumerate(cmp.groups_order):
        code = categories.get_loc(cat)
        if clf.coef_.shape[0] == 1:
            yield (
                igroup,
                clf.coef_[0] if code == clf.classes_[1] else -clf.coef_[0],
                None,
            )
        else:
            yield igroup, clf.coef_[np.flatnonzero(clf.classes_ == code)[0]], None

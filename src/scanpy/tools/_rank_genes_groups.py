"""Rank genes according to differential expression."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
import pandas as pd
from anndata import AnnData
from fast_array_utils.numba import njit
from scipy import sparse

from .. import _utils
from .. import logging as logg
from .._compat import CSBase, DaskArray
from .._settings import Default
from .._settings.presets import DETest
from .._utils import (
    _numba_thread_limit,
    check_nonnegative_integers,
    get_literal_vals,
    raise_not_implemented_error_if_backed_type,
)
from ..get import _check_mask, _get_obs_rep, aggregate
from ..get._aggregated import _chan_combine

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable
    from typing import Literal

    from numpy.typing import NDArray


type _CorrMethod = Literal["benjamini-hochberg", "bonferroni"]

_CONST_MAX_SIZE: int = 10_000_000


def _select_top_n(scores: NDArray, n_top: int):
    n_from = scores.shape[0]
    reference_indices = np.arange(n_from, dtype=int)
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]

    return global_indices


def _illico_results_to_iter(
    illico_df: pd.DataFrame,
    groups_order: NDArray,
    ireference: int | None,
):
    """Yield per-group ``(index, z, p)`` tuples from illico's long-form output.

    illico returns a DataFrame with a 2-level MultiIndex ``(pert, feature)``
    and columns including ``z_score`` and ``p_value``. We stream one group
    at a time via `pandas.Series.loc`, trusting illico_df groups are ordered
    by ``var_name``.
    """
    ref_label = None if ireference is None else groups_order[ireference]
    z_series = illico_df["z_score"]
    p_series = illico_df["p_value"]
    illico_groups = set(illico_df.index.unique(level="pert"))
    return (
        (
            group_index,
            z_series.loc[group_name].to_numpy(),
            p_series.loc[group_name].to_numpy(),
        )
        for group_index, group_name in enumerate(groups_order)
        if group_name != ref_label and group_name in illico_groups
    )


@njit
def rankdata(data: NDArray[np.number]) -> NDArray[np.float64]:
    """Parallelized version of scipy.stats.rankdata."""
    ranked = np.empty(data.shape, dtype=np.float64)
    for j in numba.prange(data.shape[1]):
        arr = np.ravel(data[:, j])
        sorter = np.argsort(arr)

        arr = arr[sorter]
        obs = np.concatenate((np.array([True]), arr[1:] != arr[:-1]))

        dense = np.empty(obs.size, dtype=np.int64)
        dense[sorter] = obs.cumsum()

        # cumulative counts of each unique value
        count = np.concatenate((np.flatnonzero(obs), np.array([len(obs)])))
        ranked[:, j] = 0.5 * (count[dense] + count[dense - 1] + 1)

    return ranked


@njit
def _tiecorrect(rankvals: NDArray[np.number]) -> NDArray[np.float64]:
    """Parallelized version of scipy.stats.tiecorrect."""
    tc = np.ones(rankvals.shape[1], dtype=np.float64)
    for j in numba.prange(rankvals.shape[1]):
        arr = np.sort(np.ravel(rankvals[:, j]))
        idx = np.flatnonzero(
            np.concatenate((np.array([True]), arr[1:] != arr[:-1], np.array([True])))
        )
        cnt = np.diff(idx).astype(np.float64)

        size = np.float64(arr.size)
        if size >= 2:
            tc[j] = 1.0 - (cnt**3 - cnt).sum() / (size**3 - size)

    return tc


def _ranks(
    x: NDArray[np.number] | CSBase,
    /,
    mask_obs: NDArray[np.bool] | None = None,
    mask_obs_rest: NDArray[np.bool] | None = None,
) -> Generator[tuple[NDArray[np.float64], int, int], None, None]:
    n_genes = x.shape[1]

    if isinstance(x, CSBase):
        merge = lambda tpl: sparse.vstack(tpl).toarray()
        adapt = lambda x: x.toarray()
    else:
        merge = np.vstack
        adapt = lambda x: x

    masked = mask_obs is not None and mask_obs_rest is not None

    if masked:
        n_cells = np.count_nonzero(mask_obs) + np.count_nonzero(mask_obs_rest)
        get_chunk = lambda x, left, right: merge((
            x[mask_obs, left:right],
            x[mask_obs_rest, left:right],
        ))
    else:
        n_cells = x.shape[0]
        get_chunk = lambda x, left, right: adapt(x[:, left:right])

    # Calculate chunk frames
    max_chunk = max(_CONST_MAX_SIZE // n_cells, 1)

    for left in range(0, n_genes, max_chunk):
        right = min(left + max_chunk, n_genes)

        ranks = rankdata(get_chunk(x, left, right))
        yield ranks, left, right


def _apply_expm1_preserving_sparsity(X, expm1_func):
    """Apply ``expm1`` to X, lazily and chunk-wise for dask, keeping sparse
    data sparse (``expm1(0) == 0``).
    """
    if isinstance(X, DaskArray):
        return X.map_blocks(
            _apply_expm1_preserving_sparsity,
            expm1_func,
            dtype=X.dtype,
            meta=X._meta,
        )
    if isinstance(X, CSBase):
        Xp = X.copy()
        Xp.data = expm1_func(Xp.data)
        return Xp
    return expm1_func(X)


@numba.njit
def _scan(
    part_n: NDArray[np.float64],
    mean: NDArray[np.float64],
    m2: NDArray[np.float64],
    j: int,
    step: int,
) -> NDArray[np.float64]:
    """Running Chan combine over the partitions for gene ``j``.

    ``acc[i]`` is the combined ``(count, mean, M2)`` of partition ``i`` together with
    every partition on the ``step`` side of it (``step=1`` forward / prefix,
    ``step=-1`` backward / suffix).
    """
    n_parts = part_n.shape[0]
    acc = np.empty((n_parts, 3))
    n = m = v = 0.0
    i = 0 if step == 1 else n_parts - 1
    for _ in range(n_parts):
        n, m, v = _chan_combine(n, m, v, part_n[i], mean[i, j], m2[i, j])
        acc[i, 0], acc[i, 1], acc[i, 2] = n, m, v
        i += step
    return acc


@njit
def _vars_rest(
    part_n: NDArray[np.float64],
    mean: NDArray[np.float64],
    m2: NDArray[np.float64],
    k: int,
) -> NDArray[np.float64]:
    """Leave-one-out variance for each selected group, parallel over genes.

    Group ``g``'s "rest" combines the forward scan up to ``g - 1`` with the backward
    scan from ``g + 1`` — every partition except ``g`` — so variances are never
    subtracted (Chan's cancellation-free combine).
    """
    n_parts, n_genes = mean.shape
    vars_rest = np.zeros((k, n_genes))
    for j in numba.prange(n_genes):
        prefix = _scan(part_n, mean, m2, j, 1)
        suffix = _scan(part_n, mean, m2, j, -1)
        for g in range(k):
            right = suffix[g + 1]  # everything after g; g + 1 < n_parts always holds
            if g >= 1:
                left = prefix[g - 1]
                n_r, _, m2_r = _chan_combine(
                    left[0], left[1], left[2], right[0], right[1], right[2]
                )
            else:
                n_r, m2_r = right[0], right[2]
            denom = n_r - 1.0 if n_r >= 2.0 else 1.0
            v = m2_r / denom
            vars_rest[g, j] = v if v > 0.0 else 0.0
    return vars_rest


class _RankGenes:
    def __init__(
        self,
        adata: AnnData,
        groups: Iterable[str] | Literal["all"],
        groupby: str,
        *,
        mask_var: NDArray[np.bool] | None = None,
        reference: Literal["rest"] | str = "rest",
        use_raw: bool = True,
        layer: str | None = None,
        comp_pts: bool = False,
    ) -> None:
        self.mask_var = mask_var
        if (base := adata.uns.get("log1p", {}).get("base")) is not None:
            self.expm1_func = lambda x: np.expm1(x * np.log(base))
        else:
            self.expm1_func = np.expm1
        self.group_col = adata.obs[groupby].array

        self.groups_order, self.groups_masks_obs = _utils.select_groups(
            adata, groups, groupby
        )

        # Singlet groups cause division by zero errors
        invalid_groups_selected = set(self.groups_order) & set(
            adata.obs[groupby].value_counts().loc[lambda x: x < 2].index
        )

        if len(invalid_groups_selected) > 0:
            msg = (
                f"Could not calculate statistics for groups {', '.join(invalid_groups_selected)} "
                "since they only contain one sample."
            )
            raise ValueError(msg)

        adata_comp = adata
        if layer is not None:
            if use_raw:
                msg = "Cannot specify `layer` and have `use_raw=True`."
                raise ValueError(msg)
            x = adata_comp.layers[layer]
        else:
            if use_raw and adata.raw is not None:
                adata_comp = adata.raw
            x = adata_comp.X
        raise_not_implemented_error_if_backed_type(x, "rank_genes_groups")

        # for correct getnnz calculation
        if isinstance(x, CSBase):
            x.eliminate_zeros()

        if self.mask_var is not None:
            self.X = x[:, self.mask_var]
            self.var_names = adata_comp.var_names[self.mask_var]

        else:
            self.X = x
            self.var_names = adata_comp.var_names

        self.ireference = None
        if reference != "rest":
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

    def _basic_stats(
        self, *, exponentiate_values: bool = False, need_var: bool = False
    ) -> None:
        """Populate per-group stats, and (in vs_rest mode) rest-group stats.

        ``need_var`` controls whether variance (per-group and per-rest) is
        computed; only the t-test family reads it. In vs_rest mode every cell
        is partitioned into its selected group or a single "remainder"
        partition (cells in no selected group), and each group's "rest" is the
        forward Chan-combine of all other partitions — a sum of non-negative
        terms, hence free of catastrophic cancellation for any group sizes.
        """
        X = (
            _apply_expm1_preserving_sparsity(self.X, self.expm1_func)
            if exponentiate_values
            else self.X
        )
        if self.ireference is None:
            self._stats_vs_rest(X, need_var=need_var)
        else:
            self._stats_vs_reference(X, need_var=need_var)

    def _aggregate_group_stats(
        self, X_used, codes: NDArray[np.int64], n_parts: int, *, need_var: bool
    ):
        """One batched :func:`scanpy.get.aggregate` over ``X_used`` grouped by
        ``codes`` (values ``0 .. n_parts-1``).

        Returns ``(mean, var, nnz)`` of shape ``(n_parts, n_genes)``,
        zero-filled for partitions with no cells. ``var`` is ``None`` unless
        ``need_var``; ``nnz`` is ``None`` unless ``self.comp_pts``.
        """
        n_genes = X_used.shape[1]
        mean = np.zeros((n_parts, n_genes))
        var = np.zeros((n_parts, n_genes)) if need_var else None
        nnz = np.zeros((n_parts, n_genes)) if self.comp_pts else None

        funcs = ["mean"]
        if need_var:
            funcs.append("var")
        if self.comp_pts:
            funcs.append("count_nonzero")
        agg_adata = AnnData(
            X=X_used,
            obs=pd.DataFrame(
                {"_g": pd.Categorical(codes, categories=range(n_parts))},
                index=pd.RangeIndex(len(codes)).astype(str),
            ),
        )
        out = aggregate(agg_adata, by="_g", func=funcs, dof=1)
        idx = out.obs_names.astype(int).to_numpy()
        mean[idx] = np.asarray(out.layers["mean"])
        if need_var:
            var[idx] = np.asarray(out.layers["var"])
        if self.comp_pts:
            nnz[idx] = np.asarray(out.layers["count_nonzero"])
        return mean, var, nnz

    def _stats_vs_reference(self, X, *, need_var: bool) -> None:
        """vs-reference: aggregate the selected-group cells only (the
        reference is itself one of the selected groups). No rest derivation.
        """
        mask = self.grouping_mask.to_numpy()
        X_used = X if mask.all() else X[mask]
        codes = pd.Categorical(self.grouping, categories=self.groups_order).codes
        k = self.groups_masks_obs.shape[0]

        self.means, self.vars, nnz = self._aggregate_group_stats(
            X_used, codes, k, need_var=need_var
        )
        if self.comp_pts:
            n_per_group = self.groups_masks_obs.sum(axis=1)
            self.pts = nnz / n_per_group[:, None]
        else:
            self.pts = None

    def _stats_vs_rest(self, X, *, need_var: bool) -> None:
        """vs-rest: partition *all* cells into the ``k`` selected-group
        partitions plus one "remainder" partition (cells in no selected group
        — non-selected groups and unassigned/NaN). Each group's "rest" is the
        forward Chan-combine of every other partition.
        """
        k = self.groups_masks_obs.shape[0]

        # Codes: each cell's selected-group index, or `k` (the remainder
        # partition) for cells in no selected group (non-selected / NaN).
        sel = pd.Categorical(self.group_col, categories=self.groups_order).codes
        codes = np.where(sel >= 0, sel, k).astype(np.int64)
        part_n = np.bincount(codes, minlength=k + 1)  # partition k == remainder
        n_sel = part_n[:k]

        mean, var, nnz = self._aggregate_group_stats(X, codes, k + 1, need_var=need_var)

        # Selected-group arm of the test (the remainder partition is excluded).
        self.means = mean[:k]
        self.vars = var[:k] if need_var else None
        self.pts = nnz[:k] / n_sel[:, None] if self.comp_pts else None

        # M2 = var * (n - 1); forced to 0 for partitions with <= 1 cell so a
        # singleton remainder (aggregate var undefined there) is harmless.
        if need_var:
            with np.errstate(invalid="ignore"):
                M2 = var * (part_n[:, None] - 1)
            M2[part_n <= 1] = 0.0
        else:
            M2 = None

        self._derive_rest_stats(part_n, mean, M2, nnz, k, need_var=need_var)

    def _derive_rest_stats(
        self, part_n, mean, M2, nnz, k: int, *, need_var: bool
    ) -> None:
        """Set ``means_rest``/``vars_rest``/``pts_rest`` for each selected group ``g``
        (statistics over every cell *not* in ``g``).

        ``means_rest`` and ``pts_rest`` are linear in the partitions, so they are the
        exact total-minus-group difference. Variance would lose precision under such
        subtraction, so ``vars_rest`` uses the cancellation-free Chan leave-one-out
        scan (:func:`_vars_rest`) — only the variance-based tests request it.
        """
        n_rest = (self.X.shape[0] - part_n[:k])[:, None]
        total = (part_n[:, None] * mean).sum(axis=0)
        self.means_rest = (total - part_n[:k, None] * mean[:k]) / n_rest
        self.vars_rest = (
            _vars_rest(np.ascontiguousarray(part_n, dtype=np.float64), mean, M2, k)
            if need_var
            else None
        )
        self.pts_rest = (nnz.sum(axis=0) - nnz[:k]) / n_rest if self.comp_pts else None

    def t_test(
        self, method: Literal["t-test", "t-test_overestim_var"]
    ) -> Generator[tuple[int, NDArray[np.floating], NDArray[np.floating]], None, None]:
        from scipy import stats

        for group_index, (mask_obs, mean_group, var_group) in enumerate(
            zip(self.groups_masks_obs, self.means, self.vars, strict=True)
        ):
            if self.ireference is not None and group_index == self.ireference:
                continue

            ns_group = np.count_nonzero(mask_obs)

            if self.ireference is not None:
                mean_rest = self.means[self.ireference]
                var_rest = self.vars[self.ireference]
                ns_other = np.count_nonzero(self.groups_masks_obs[self.ireference])
            else:
                mean_rest = self.means_rest[group_index]
                var_rest = self.vars_rest[group_index]
                ns_other = self.X.shape[0] - ns_group

            if method == "t-test":
                ns_rest = ns_other
            elif method == "t-test_overestim_var":
                # hack for overestimating the variance for small groups
                ns_rest = ns_group
            else:
                msg = "Method does not exist."
                raise ValueError(msg)

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

    def wilcoxon(
        self, *, tie_correct: bool
    ) -> Generator[tuple[int, NDArray[np.floating], NDArray[np.floating]], None, None]:
        from scipy import stats

        n_genes = self.X.shape[1]
        # First loop: Loop over all genes
        if self.ireference is not None:
            # initialize space for z-scores
            scores = np.zeros(n_genes)
            # initialize space for tie correction coefficients
            tc_coef = np.zeros(n_genes) if tie_correct else 1

            for group_index, mask_obs in enumerate(self.groups_masks_obs):
                if group_index == self.ireference:
                    continue

                mask_obs_rest = self.groups_masks_obs[self.ireference]

                n_active = np.count_nonzero(mask_obs)
                m_active = np.count_nonzero(mask_obs_rest)

                if n_active <= 25 or m_active <= 25:
                    logg.hint(
                        "Few observations in a group for "
                        "normal approximation (<=25). Lower test accuracy."
                    )

                # Calculate rank sums for each chunk for the current mask
                for ranks, left, right in _ranks(self.X, mask_obs, mask_obs_rest):
                    scores[left:right] = ranks[0:n_active, :].sum(axis=0)
                    if tie_correct:
                        tc_coef[left:right] = _tiecorrect(ranks)

                std_dev = np.sqrt(
                    tc_coef * n_active * m_active * (n_active + m_active + 1) / 12.0
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
            n_groups = self.groups_masks_obs.shape[0]
            scores = np.zeros((n_groups, n_genes))
            n_cells = self.X.shape[0]

            if tie_correct:
                tc_coef = np.zeros((n_groups, n_genes))

            for ranks, left, right in _ranks(self.X):
                if tie_correct:
                    tc_coef[:, left:right] = _tiecorrect(ranks)
                # sum up adjusted_ranks to calculate W_m,n
                for group_index, mask_obs in enumerate(self.groups_masks_obs):
                    scores[group_index, left:right] = ranks[mask_obs, :].sum(axis=0)

            for group_index, mask_obs in enumerate(self.groups_masks_obs):
                n_active = np.count_nonzero(mask_obs)

                coef = tc_coef[group_index] if tie_correct else 1

                std_dev = np.sqrt(
                    coef * n_active * (n_cells - n_active) * (n_cells + 1) / 12.0
                )

                scores[group_index, :] = (
                    scores[group_index, :] - (n_active * (n_cells + 1) / 2.0)
                ) / std_dev
                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores[group_index, :]))

                yield group_index, scores[group_index], pvals

    def logreg(
        self, **kwds
    ) -> Generator[tuple[int, NDArray[np.floating], None], None, None]:
        # if reference is not set, then the groups listed will be compared to the rest
        # if reference is set, then the groups listed will be compared only to the other groups listed
        from sklearn.linear_model import LogisticRegression

        # Indexing with a series causes issues, possibly segfault
        x = self.X[self.grouping_mask.to_numpy(), :]

        if len(self.groups_order) == 1:
            msg = "Cannot perform logistic regression on a single cluster."
            raise ValueError(msg)

        clf = LogisticRegression(**kwds)
        clf.fit(x, self.grouping.cat.codes)
        scores_all = clf.coef_
        # not all codes necessarily appear in data
        existing_codes = np.unique(self.grouping.cat.codes)
        for igroup, cat in enumerate(self.groups_order):
            if len(self.groups_order) <= 2:  # binary logistic regression
                scores = scores_all[0]
            else:
                # cat code is index of cat value in .categories
                cat_code: int = np.argmax(self.grouping.cat.categories == cat)
                # index of scores row is index of cat code in array of existing codes
                scores_idx: int = np.argmax(existing_codes == cat_code)
                scores = scores_all[scores_idx]
            yield igroup, scores, None

            if len(self.groups_order) <= 2:
                break

    def _run_illico(self, *, tie_correct: bool):
        """Invoke `illico.asymptotic_wilcoxon` on `self.X` / `self.group_col`."""
        from illico import asymptotic_wilcoxon

        return asymptotic_wilcoxon(
            AnnData(
                X=self.X,
                var=pd.DataFrame(index=self.var_names),
                obs=pd.DataFrame(
                    index=pd.RangeIndex(self.X.shape[0]).astype("str"),
                    # This self.group_col means illico will run tests against
                    # *all* data instead of what's in self.groups_order as
                    # controlled by the `groups` arg.
                    # TODO: Only run the subset once illico supports a `groups` argument
                    data={"group": self.group_col},
                ),
            ),
            reference=self.groups_order[self.ireference]
            if self.ireference is not None
            else None,
            group_keys="group",
            return_as_scanpy=False,
            is_log1p=True,
            tie_correct=tie_correct,
            use_continuity=False,
            alternative="two-sided",
            use_rust=False,
        )

    def compute_statistics(
        self,
        method: DETest,
        *,
        corr_method: _CorrMethod,
        n_genes_user: int | None,
        rankby_abs: bool,
        tie_correct: bool,
        mean_in_log_space: bool,
        **kwds,
    ) -> None:
        if method in {"t-test", "t-test_overestim_var"}:
            self._basic_stats(exponentiate_values=not mean_in_log_space, need_var=True)
            generate_test_results = self.t_test(method)
        elif "wilcoxon" in method:
            if "illico" in method:
                illico_df = self._run_illico(tie_correct=tie_correct)
                generate_test_results = _illico_results_to_iter(
                    illico_df,
                    self.groups_order,
                    self.ireference,
                )
            else:
                generate_test_results = self.wilcoxon(tie_correct=tie_correct)
            # Wilcoxon paths only need means (for fold-change); skip var.
            self._basic_stats(exponentiate_values=not mean_in_log_space, need_var=False)
        elif method == "logreg":
            generate_test_results = self.logreg(**kwds)

        self.stats = _build_stats_dataframe(
            self,
            generate_test_results,
            corr_method=corr_method,
            n_genes_user=n_genes_user,
            rankby_abs=rankby_abs,
            mean_in_log_space=mean_in_log_space,
        )


def _build_stats_dataframe(
    rg,
    results,
    *,
    corr_method: _CorrMethod,
    n_genes_user: int | None,
    rankby_abs: bool,
    mean_in_log_space: bool,
):
    """Drain the per-group ``(group_index, scores, pvals)`` iterator into a
    wide-form ``(group, statistic)`` MultiIndex DataFrame: top-N selection,
    multiple-testing correction, and (when ``rg.means`` is set) log2
    fold-change. Read-only on ``rg``.
    """
    from statsmodels.stats.multitest import multipletests

    n_genes_total = rg.X.shape[1]
    df: pd.DataFrame | None = None

    for group_index, scores, pvals in results:
        group_name = str(rg.groups_order[group_index])

        if n_genes_user is not None:
            scores_sort = np.abs(scores) if rankby_abs else scores
            global_indices = _select_top_n(scores_sort, n_genes_user)
            first_col = "names"
        else:
            global_indices = slice(None)
            first_col = "scores"

        if df is None:
            idx = pd.MultiIndex.from_tuples([(group_name, first_col)])
            df = pd.DataFrame(columns=idx)

        if n_genes_user is not None:
            df[group_name, "names"] = rg.var_names[global_indices]
        df[group_name, "scores"] = scores[global_indices]

        if pvals is not None:
            df[group_name, "pvals"] = pvals[global_indices]
            if corr_method == "benjamini-hochberg":
                pvals[np.isnan(pvals)] = 1
                _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
            elif corr_method == "bonferroni":
                pvals_adj = np.minimum(pvals * n_genes_total, 1.0)
            df[group_name, "pvals_adj"] = pvals_adj[global_indices]

        if rg.means is not None:
            mean_group = rg.means[group_index]
            mean_rest = (
                rg.means_rest[group_index]
                if rg.ireference is None
                else rg.means[rg.ireference]
            )
            foldchanges = (
                (rg.expm1_func(mean_group) + 1e-9) / (rg.expm1_func(mean_rest) + 1e-9)
                if mean_in_log_space
                else (mean_group + 1e-9) / (mean_rest + 1e-9)
            )  # add small value to avoid zeros
            df[group_name, "logfoldchanges"] = np.log2(foldchanges[global_indices])

    if df is not None and n_genes_user is None:
        df.index = rg.var_names
    return df


def rank_genes_groups(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    groupby: str,
    *,
    mask_var: NDArray[np.bool] | str | None | Default = Default(
        preset=("rank_genes_groups", "mask_var")
    ),
    use_raw: bool | None = None,
    groups: Literal["all"] | Iterable[str] = "all",
    reference: str = "rest",
    n_genes: int | None = None,
    rankby_abs: bool = False,
    pts: bool = False,
    key_added: str | None = None,
    copy: bool = False,
    method: DETest | Default = Default(preset=("rank_genes_groups", "method")),
    corr_method: _CorrMethod = "benjamini-hochberg",
    tie_correct: bool = False,
    layer: str | None = None,
    mean_in_log_space: bool | Default = Default(
        preset=("rank_genes_groups", "mean_in_log_space")
    ),
    **kwds,
) -> AnnData | None:
    r"""Rank genes for characterizing groups.

    Expects logarithmized data.

    .. array-support:: tl.rank_genes_groups

    ..  warning::

        Comparing between cells leads to highly inflated p-values,
        since cells are not independent observations :cite:p:`Squair2021`.
        Especially in single-cell data, consider instead to use more appropriate methods such as combining pseudobulking with :doc:`pydeseq2:index`.

        :func:`decoupler.pp.pseudobulk` or :func:`scanpy.get.aggregate` can be used to aggregate samples for pseudobulking.
        Ours is a bit more verbose, but supports :doc:`dask:index` arrays for improved performance.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    mask_var
        Select subset of genes to use in statistical tests.
    use_raw
        Use `raw` attribute of `adata` if present. The default behavior is to use `raw` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    groups
        Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
        shall be restricted, or `'all'` (default), for all groups. Note that if
        `reference='rest'` all groups will still be used as the reference, not
        just those specified in `groups`.
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
        `'logreg'` uses logistic regression. See :cite:t:`Ntranos2019`,
        `here <https://github.com/scverse/scanpy/issues/95>`__ and `here
        <https://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,
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
    copy
        Whether to copy `adata` or modify it inplace.
    mean_in_log_space
        Whether to do :math:`\log(\operatorname{mean}(e^x))` (`False`)
        or :math:`\log(e^{\operatorname{mean}(x)})` (`True`).
        The former is accurate, while the latter is a faster approximation
        that underestimates this accurate result in the presence of many outliers.
    kwds
        Are passed to test methods. Currently this affects only parameters that
        are passed to :class:`sklearn.linear_model.LogisticRegression`.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.uns['rank_genes_groups' | key_added]['names']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['scores']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the z-score
        underlying the computation of a p-value for each gene for each
        group. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['logfoldchanges']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the log2
        fold change for each gene for each group. Ordered according to
        scores. Only provided if method is 't-test' like.
        Note: if `mean_in_log_space=True`, this is an approximation calculated from mean-log values.
    `adata.uns['rank_genes_groups' | key_added]['pvals']` : structured :class:`numpy.ndarray` (dtype `float`)
        p-values.
    `adata.uns['rank_genes_groups' | key_added]['pvals_adj']` : structured :class:`numpy.ndarray` (dtype `float`)
        Corrected p-values.
    `adata.uns['rank_genes_groups' | key_added]['pts']` : :class:`pandas.DataFrame` (dtype `float`)
        Fraction of cells expressing the genes for each group.
    `adata.uns['rank_genes_groups' | key_added]['pts_rest']` : :class:`pandas.DataFrame` (dtype `float`)
        Only if `reference` is set to `'rest'`.
        Fraction of cells from the union of the rest of each group
        expressing the genes.

    Notes
    -----
    There are slight inconsistencies depending on whether sparse
    or dense data are passed. See `here <https://github.com/scverse/scanpy/blob/main/tests/test_rank_genes_groups.py>`__.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")
    >>> # to visualize the results
    >>> sc.pl.rank_genes_groups(adata)

    """
    from scanpy import settings

    if isinstance(mask_var, Default):
        mask_var = settings.preset.rank_genes_groups.mask_var
    if isinstance(mean_in_log_space, Default):
        mean_in_log_space = settings.preset.rank_genes_groups.mean_in_log_space
    if method is None or isinstance(method, Default):
        method = settings.preset.rank_genes_groups.method

    mask_var = _check_mask(adata, mask_var, "var")

    if use_raw is None:
        use_raw = adata.raw is not None
    elif use_raw is True and adata.raw is None:
        msg = "Received `use_raw=True`, but `adata.raw` is empty."
        raise ValueError(msg)

    if "only_positive" in kwds:
        rankby_abs = not kwds.pop("only_positive")  # backwards compat

    start = logg.info("ranking genes")
    if method not in (avail_methods := get_literal_vals(DETest)):
        msg = f"Method must be one of {avail_methods}."
        raise ValueError(msg)

    avail_corr = {"benjamini-hochberg", "bonferroni"}
    if corr_method not in avail_corr:
        msg = f"Correction method must be one of {avail_corr}."
        raise ValueError(msg)

    adata = adata.copy() if copy else adata
    _utils.sanitize_anndata(adata)
    # for clarity, rename variable
    if groups == "all":
        groups_order = "all"
    elif isinstance(groups, str | int):
        msg = "Specify a sequence of groups"
        raise ValueError(msg)
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if reference != "rest" and reference not in set(groups_order):
            groups_order += [reference]
    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        cats = adata.obs[groupby].cat.categories.tolist()
        msg = f"reference = {reference} needs to be one of groupby = {cats}."
        raise ValueError(msg)

    if key_added is None:
        key_added = "rank_genes_groups"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = dict(
        groupby=groupby,
        reference=reference,
        method=method,
        use_raw=use_raw,
        layer=layer,
        corr_method=corr_method,
    )

    test_obj = _RankGenes(
        adata,
        groups_order,
        groupby,
        mask_var=mask_var,
        reference=reference,
        use_raw=use_raw,
        layer=layer,
        comp_pts=pts,
    )

    if check_nonnegative_integers(test_obj.X) and method != "logreg":
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

    logg.debug(f"consider {groupby!r} groups:")
    logg.debug(f"with sizes: {np.count_nonzero(test_obj.groups_masks_obs, axis=1)}")

    with _numba_thread_limit(settings.n_jobs if method == "wilcoxon" else None):
        test_obj.compute_statistics(
            method,
            corr_method=corr_method,
            n_genes_user=n_genes_user,
            rankby_abs=rankby_abs,
            tie_correct=tie_correct,
            mean_in_log_space=mean_in_log_space,
            **kwds,
        )

    if test_obj.pts is not None:
        groups_names = [str(name) for name in test_obj.groups_order]
        adata.uns[key_added]["pts"] = pd.DataFrame(
            test_obj.pts.T, index=test_obj.var_names, columns=groups_names
        )
    if test_obj.pts_rest is not None:
        adata.uns[key_added]["pts_rest"] = pd.DataFrame(
            test_obj.pts_rest.T, index=test_obj.var_names, columns=groups_names
        )

    test_obj.stats.columns = test_obj.stats.columns.swaplevel()

    dtypes = {
        "names": "O",
        "scores": "float32",
        "logfoldchanges": "float32",
        "pvals": "float64",
        "pvals_adj": "float64",
    }

    for col in test_obj.stats.columns.levels[0]:
        adata.uns[key_added][col] = test_obj.stats[col].to_records(
            index=False, column_dtypes=dtypes[col]
        )

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added to `.uns[{key_added!r}]`\n"
            "    'names', sorted np.recarray to be indexed by group ids\n"
            "    'scores', sorted np.recarray to be indexed by group ids\n"
            + (
                "    'logfoldchanges', sorted np.recarray to be indexed by group ids\n"
                "    'pvals', sorted np.recarray to be indexed by group ids\n"
                "    'pvals_adj', sorted np.recarray to be indexed by group ids"
                if method in {"t-test", "t-test_overestim_var", "wilcoxon"}
                else ""
            )
        ),
    )
    return adata if copy else None


def _calc_frac(x: NDArray[np.number] | CSBase, /) -> NDArray[np.float64]:
    n_nonzero = (
        x.getnnz(axis=0) if isinstance(x, CSBase) else np.count_nonzero(x, axis=0)
    )
    return n_nonzero / x.shape[0]


def filter_rank_genes_groups(  # noqa: PLR0912
    adata: AnnData,
    *,
    key: str | None = None,
    groupby: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    key_added: str = "rank_genes_groups_filtered",
    min_in_group_fraction: float = 0.25,
    min_fold_change: float = 1,
    max_out_group_fraction: float = 0.5,
    compare_abs: bool = False,
) -> None:
    """Filter out genes based on two criteria.

    1. log fold change and
    2. fraction of genes expressing the
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
    layer
    use_raw
    key_added
    min_in_group_fraction
    min_fold_change
    max_out_group_fraction
    compare_abs
        If `True`, compare absolute values of log fold change with `min_fold_change`.

    Returns
    -------
    Same output as :func:`scanpy.tl.rank_genes_groups` but with filtered genes names set to `nan`.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")
    >>> sc.tl.filter_rank_genes_groups(adata, min_fold_change=3)
    >>> # visualize results
    >>> sc.pl.rank_genes_groups(adata, key="rank_genes_groups_filtered")
    >>> # visualize results using dotplot
    >>> sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered")

    """
    if key is None:
        key = "rank_genes_groups"

    if groupby is None:
        groupby = adata.uns[key]["params"]["groupby"]

    if use_raw is None:
        use_raw = adata.uns[key]["params"]["use_raw"] if layer is None else False

    x = _get_obs_rep(adata, use_raw=use_raw, layer=layer)

    same_params = (
        adata.uns[key]["params"]["groupby"] == groupby
        and adata.uns[key]["params"]["reference"] == "rest"
        and adata.uns[key]["params"]["use_raw"] == use_raw
    )

    use_logfolds = same_params and "logfoldchanges" in adata.uns[key]
    use_fraction = same_params and "pts_rest" in adata.uns[key]

    # convert structured numpy array into DataFrame
    gene_names = pd.DataFrame(adata.uns[key]["names"])

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
        fold_change_matrix = pd.DataFrame(adata.uns[key]["logfoldchanges"])
    else:
        fold_change_matrix = pd.DataFrame(
            np.zeros(gene_names.shape),
            columns=gene_names.columns,
            index=gene_names.index,
        )

        if (base := adata.uns.get("log1p", {}).get("base")) is not None:
            expm1_func = lambda x: np.expm1(x * np.log(base))
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
        var_names = gene_names[cluster].array

        if not use_logfolds or not use_fraction:
            var_idx = (adata.raw if use_raw else adata).var_names.get_indexer(var_names)
            sub_x = x[:, var_idx]
            in_group = (adata.obs[groupby] == cluster).to_numpy()
            x_in = sub_x[in_group]
            x_out = sub_x[~in_group]

        if use_fraction:
            fraction_in_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]["pts"][cluster].loc[var_names].array
            )
            fraction_out_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]["pts_rest"][cluster].loc[var_names].array
            )
        else:
            fraction_in_cluster_matrix.loc[:, cluster] = _calc_frac(x_in)
            fraction_out_cluster_matrix.loc[:, cluster] = _calc_frac(x_out)

        if not use_logfolds:
            # compute mean value
            mean_in_cluster = np.ravel(x_in.mean(0))
            mean_out_cluster = np.ravel(x_out.mean(0))
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
    adata.uns[key_added]["names"] = gene_names.to_records(index=False)

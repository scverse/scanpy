from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from typing import Literal

    from numpy.typing import NDArray

    from ._comparison import _Comparison
    from ._scorers import _TestResult
    from ._stats import _GroupStats


type _CorrMethod = Literal["benjamini-hochberg", "bonferroni"]


class _GroupResult(NamedTuple):
    group_index: int
    scores: NDArray[np.floating]
    pvals: NDArray[np.floating] | None
    pvals_adj: NDArray[np.floating] | None
    foldchanges: NDArray[np.floating] | None


def _select_top_n(scores: NDArray, n_top: int) -> NDArray[np.intp]:
    partition = np.argpartition(scores, -n_top)[-n_top:]
    return partition[np.argsort(scores[partition])[::-1]]


def _adjust_pvals(pvals: NDArray, corr_method: _CorrMethod) -> NDArray:
    from statsmodels.stats.multitest import multipletests

    if corr_method == "benjamini-hochberg":
        pvals_no_nan = np.where(np.isnan(pvals), 1.0, pvals)
        _, pvals_adj, _, _ = multipletests(pvals_no_nan, alpha=0.05, method="fdr_bh")
        return pvals_adj
    if corr_method == "bonferroni":
        return np.minimum(pvals * pvals.shape[0], 1.0)
    msg = f"Correction method must be one of {{'benjamini-hochberg', 'bonferroni'}}, not {corr_method!r}."
    raise ValueError(msg)


def _group_results(
    cmp: _Comparison,
    lfc_stats: _GroupStats | None,
    results: Iterable[_TestResult],
    *,
    corr_method: _CorrMethod,
    mean_in_log_space: bool,
) -> Iterator[_GroupResult]:
    """Add multiple-testing correction and fold changes to each group’s scores.

    Fold changes are returned as ratios (``log2`` is taken by the caller).
    """
    for group_index, scores, pvals in results:
        pvals_adj = None if pvals is None else _adjust_pvals(pvals, corr_method)
        foldchanges = None
        if lfc_stats is not None:
            mean_group = lfc_stats.means[group_index]
            mean_rest = (
                lfc_stats.means_rest[group_index]
                if cmp.ireference is None
                else lfc_stats.means[cmp.ireference]
            )
            foldchanges = (
                (cmp.expm1_func(mean_group) + 1e-9) / (cmp.expm1_func(mean_rest) + 1e-9)
                if mean_in_log_space
                else (mean_group + 1e-9) / (mean_rest + 1e-9)
            )  # add small value to avoid zeros
        yield _GroupResult(group_index, scores, pvals, pvals_adj, foldchanges)


def _legacy_stats_frame(
    cmp: _Comparison,
    group_results: Iterable[_GroupResult],
    *,
    n_genes_user: int | None,
    rankby_abs: bool,
) -> pd.DataFrame | None:
    """Build the wide ``(group, statistic)`` frame stored by `rank_genes_groups`."""
    cols: dict[tuple[str, str], NDArray] = {}

    for res in group_results:
        group_name = str(cmp.groups_order[res.group_index])

        if n_genes_user is not None:
            scores_sort = np.abs(res.scores) if rankby_abs else res.scores
            global_indices = _select_top_n(scores_sort, n_genes_user)
            cols[group_name, "names"] = cmp.var_names[global_indices]
        else:
            global_indices = slice(None)
        cols[group_name, "scores"] = res.scores[global_indices]

        if res.pvals is not None:
            cols[group_name, "pvals"] = res.pvals[global_indices]
            cols[group_name, "pvals_adj"] = res.pvals_adj[global_indices]

        if res.foldchanges is not None:
            cols[group_name, "logfoldchanges"] = np.log2(
                res.foldchanges[global_indices]
            )

    if not cols:
        return None
    df = pd.DataFrame(cols)
    if n_genes_user is None:
        df.index = cmp.var_names
    return df


def _marker_table(
    cmp: _Comparison,
    stats: _GroupStats,
    group_results: Iterable[_GroupResult],
    *,
    n_top: int | None,
) -> pd.DataFrame:
    """Build the long results table returned by the `scanpy.tl.markers` functions."""
    n_genes = cmp.x.shape[1]
    n_top = n_genes if n_top is None else min(n_top, n_genes)
    frames: list[pd.DataFrame] = []
    for res in group_results:
        idx = _select_top_n(res.scores, n_top)
        frac_reference = (
            stats.pts_rest[res.group_index]
            if cmp.ireference is None
            else stats.pts[cmp.ireference]
        )
        columns: dict[str, object] = dict(
            group=str(cmp.groups_order[res.group_index]),
            reference=cmp.reference,
            gene=cmp.var_names[idx],
            score=res.scores[idx],
            log_fc=np.log2(res.foldchanges[idx]),
        )
        if res.pvals is not None:
            columns.update(p_value=res.pvals[idx], adj_p_value=res.pvals_adj[idx])
        columns.update(
            frac_group=stats.pts[res.group_index][idx],
            frac_reference=frac_reference[idx],
        )
        frames.append(pd.DataFrame(columns))
    df = pd.concat(frames, ignore_index=True)
    groups = [str(g) for g in cmp.groups_order]
    df["group"] = pd.Categorical(
        df["group"], categories=[g for g in groups if g in set(df["group"])]
    )
    df["reference"] = df["reference"].astype("category")
    return df

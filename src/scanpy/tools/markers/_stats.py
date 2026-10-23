from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from anndata import AnnData

from ..._compat import CSBase, DaskArray
from ..._utils import dim_acc
from ...get import aggregate
from ._kernels import _vars_rest

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from ._comparison import _Comparison


@dataclass(frozen=True, eq=False)
class _GroupStats:
    """Per-group statistics and, in vs-rest mode, the statistics of each group’s rest."""

    means: NDArray[np.floating]
    vars: NDArray[np.floating] | None
    pts: NDArray[np.floating] | None
    means_rest: NDArray[np.floating] | None
    vars_rest: NDArray[np.floating] | None
    pts_rest: NDArray[np.floating] | None


def _apply_expm1_preserving_sparsity(x, expm1_func: Callable[[NDArray], NDArray]):
    """Apply ``expm1`` to ``x`` while keeping sparse data sparse.

    Applied lazily and chunk-wise for dask; ``expm1(0) == 0`` preserves sparsity.
    """
    if isinstance(x, DaskArray):
        return x.map_blocks(
            _apply_expm1_preserving_sparsity,
            expm1_func,
            dtype=x.dtype,
            meta=x._meta,
        )
    if isinstance(x, CSBase):
        xp = x.copy()
        xp.data = expm1_func(xp.data)
        return xp
    return expm1_func(x)


def _group_stats(
    cmp: _Comparison,
    *,
    exponentiate_values: bool = False,
    need_var: bool = False,
    comp_pts: bool = False,
) -> _GroupStats:
    """Compute per-group stats, and (in vs_rest mode) rest-group stats.

    ``need_var`` controls whether variance (per-group and per-rest) is
    computed; only the t-test family reads it. In vs_rest mode every cell
    is assigned to its selected group or a single "remainder" group (cells
    in no selected group), and each group's "rest" is the forward
    Chan-combine of all other groups — a sum of non-negative terms, hence
    free of catastrophic cancellation for any group sizes.
    """
    x = (
        _apply_expm1_preserving_sparsity(cmp.x, cmp.expm1_func)
        if exponentiate_values
        else cmp.x
    )
    if cmp.ireference is None:
        return _stats_vs_rest(cmp, x, need_var=need_var, comp_pts=comp_pts)
    return _stats_vs_reference(cmp, x, need_var=need_var, comp_pts=comp_pts)


def _aggregate_group_stats(
    x_used,
    codes: NDArray[np.int64],
    n_groups: int,
    *,
    need_var: bool,
    comp_pts: bool,
) -> tuple[NDArray, NDArray | None, NDArray | None]:
    """Aggregate ``x_used`` in one batched :func:`scanpy.get.aggregate`.

    Grouped by ``codes`` (values ``0 .. n_groups-1``). Returns ``(mean, var,
    nnz)`` of shape ``(n_groups, n_genes)``, zero-filled for groups with
    no cells. ``var`` is ``None`` unless ``need_var``; ``nnz`` is ``None``
    unless ``comp_pts``.
    """
    n_genes = x_used.shape[1]
    mean = np.zeros((n_groups, n_genes))
    var = np.zeros((n_groups, n_genes)) if need_var else None
    nnz = np.zeros((n_groups, n_genes)) if comp_pts else None

    funcs = ["mean"]
    if need_var:
        funcs.append("var")
    if comp_pts:
        funcs.append("count_nonzero")
    agg_adata = AnnData(
        X=x_used,
        obs=pd.DataFrame(
            {"_g": pd.Categorical(codes, categories=range(n_groups))},
            index=pd.RangeIndex(len(codes)).astype(str),
        ),
    )
    out = aggregate(agg_adata, by=dim_acc("_g", dim="obs"), func=funcs, dof=1)
    idx = out.obs_names.astype(int).to_numpy()
    mean[idx] = np.asarray(out.layers["mean"])
    if need_var:
        var[idx] = np.asarray(out.layers["var"])
    if comp_pts:
        nnz[idx] = np.asarray(out.layers["count_nonzero"])
    return mean, var, nnz


def _stats_vs_reference(
    cmp: _Comparison, x, *, need_var: bool, comp_pts: bool
) -> _GroupStats:
    """Aggregate the selected-group cells only (vs-reference; no rest derivation).

    The reference is itself one of the selected groups.
    """
    mask = cmp.grouping_mask
    x_used = x if mask.all() else x[mask]
    codes = pd.Index(cmp.groups_order).get_indexer(cmp.grouping)

    means, vars_, nnz = _aggregate_group_stats(
        x_used, codes, cmp.n_groups, need_var=need_var, comp_pts=comp_pts
    )
    pts = None
    if comp_pts:
        n_per_group = cmp.groups_masks_obs.sum(axis=1)
        pts = nnz / n_per_group[:, None]
    return _GroupStats(
        means=means,
        vars=vars_,
        pts=pts,
        means_rest=None,
        vars_rest=None,
        pts_rest=None,
    )


def _stats_vs_rest(
    cmp: _Comparison, x, *, need_var: bool, comp_pts: bool
) -> _GroupStats:
    """Assign every cell to one of the ``k`` selected groups or a remainder group.

    The remainder group holds cells in no selected group (non-selected
    groups and unassigned/NaN). Each group's "rest" is the forward
    Chan-combine of every other group.
    """
    k = cmp.n_groups

    # each cell's selected-group index, or `k` (the remainder group) for
    # cells in no selected group (non-selected / NaN)
    sel = pd.Index(cmp.groups_order).get_indexer(cmp.labels)
    codes = np.where(sel >= 0, sel, k).astype(np.int64)
    group_counts = np.bincount(codes, minlength=k + 1)  # group k == remainder
    n_sel = group_counts[:k]

    mean, var, nnz = _aggregate_group_stats(
        x, codes, k + 1, need_var=need_var, comp_pts=comp_pts
    )

    # m2 = var * (n - 1); forced to 0 for groups with <= 1 cell so a
    # singleton remainder (aggregate var undefined there) is harmless
    if need_var:
        with np.errstate(invalid="ignore"):
            m2 = var * (group_counts[:, None] - 1)
        m2[group_counts <= 1] = 0.0
    else:
        m2 = None

    n_rest = (x.shape[0] - group_counts[:k])[:, None]
    total = (group_counts[:, None] * mean).sum(axis=0)
    return _GroupStats(
        means=mean[:k],
        vars=var[:k] if need_var else None,
        pts=nnz[:k] / n_sel[:, None] if comp_pts else None,
        means_rest=(total - group_counts[:k, None] * mean[:k]) / n_rest,
        vars_rest=(
            _vars_rest(
                np.ascontiguousarray(group_counts, dtype=np.float64), mean, m2, k
            )
            if need_var
            else None
        ),
        pts_rest=(nnz.sum(axis=0) - nnz[:k]) / n_rest if comp_pts else None,
    )

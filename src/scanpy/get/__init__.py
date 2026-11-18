"""Get data from AnnData."""

from __future__ import annotations

from ._aggregated import aggregate
from ._out import _set_arr, _write_out
from ._types import _Rep
from ._use import _check_mask, _get_arr, _get_vec
from .get import obs_df, pca, rank_genes_groups_df, var_df

__all__ = [
    "_Rep",
    "_check_mask",
    "_get_arr",
    "_get_vec",
    "_set_arr",
    "_write_out",
    "aggregate",
    "obs_df",
    "pca",
    "rank_genes_groups_df",
    "var_df",
]

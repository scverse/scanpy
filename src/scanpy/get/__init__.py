"""Get data from AnnData."""

from __future__ import annotations

from ._aggregated import aggregate
from .get import (
    _check_mask,
    _get_arr,
    _get_vec,
    _get_vec_compat,
    _Rep,
    _set_arr,
    _write_out,
    obs_df,
    pca,
    rank_genes_groups_df,
    var_df,
)

__all__ = [
    "_Rep",
    "_check_mask",
    "_get_arr",
    "_get_vec",
    "_get_vec_compat",
    "_set_arr",
    "_write_out",
    "aggregate",
    "obs_df",
    "pca",
    "rank_genes_groups_df",
    "var_df",
]

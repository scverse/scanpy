"""Get data from AnnData."""

from __future__ import annotations

from ._aggregated import aggregate
from .get import (
    _check_mask,
    _get_arr,
    _get_vec,
    _Rep,
    _set_obs_rep,
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
    "_set_obs_rep",
    "aggregate",
    "obs_df",
    "pca",
    "rank_genes_groups_df",
    "var_df",
]

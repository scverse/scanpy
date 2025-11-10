"""Get data from AnnData."""

from __future__ import annotations

from ._aggregated import aggregate
from .get import (
    _check_mask,
    _get_obs_rep,
    _ObsRep,
    _set_obs_rep,
    obs_df,
    rank_genes_groups_df,
    var_df,
)

__all__ = [
    "_ObsRep",
    "_check_mask",
    "_get_obs_rep",
    "_set_obs_rep",
    "aggregate",
    "obs_df",
    "rank_genes_groups_df",
    "var_df",
]

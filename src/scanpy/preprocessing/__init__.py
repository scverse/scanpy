"""Preprocessing functions."""

from __future__ import annotations

from ..neighbors import neighbors
from ._combat import combat
from ._deprecated.highly_variable_genes import filter_genes_dispersion
from ._deprecated.sampling import subsample
from ._highly_variable_genes import highly_variable_genes
from ._normalization import normalize_total
from ._pca import pca
from ._qc import calculate_qc_metrics
from ._recipes import recipe_seurat, recipe_weinreb17, recipe_zheng17
from ._scale import scale
from ._scrublet import scrublet, scrublet_simulate_doublets
from ._simple import (
    downsample_counts,
    filter_cells,
    filter_genes,
    log1p,
    normalize_per_cell,
    regress_out,
    sample,
    sqrt,
)

__all__ = [
    "calculate_qc_metrics",
    "combat",
    "downsample_counts",
    "filter_cells",
    "filter_genes",
    "filter_genes_dispersion",
    "highly_variable_genes",
    "log1p",
    "neighbors",
    "normalize_per_cell",
    "normalize_total",
    "pca",
    "recipe_seurat",
    "recipe_weinreb17",
    "recipe_zheng17",
    "regress_out",
    "sample",
    "scale",
    "scrublet",
    "scrublet_simulate_doublets",
    "sqrt",
    "subsample",
]

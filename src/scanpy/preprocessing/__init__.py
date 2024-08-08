from __future__ import annotations

from ..neighbors import neighbors
from ._combat import combat
from ._deprecated.highly_variable_genes import filter_genes_dispersion
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
    sqrt,
    subsample,
)

__all__ = [
    "neighbors",
    "combat",
    "filter_genes_dispersion",
    "highly_variable_genes",
    "normalize_total",
    "pca",
    "calculate_qc_metrics",
    "recipe_seurat",
    "recipe_weinreb17",
    "recipe_zheng17",
    "scrublet",
    "scrublet_simulate_doublets",
    "downsample_counts",
    "filter_cells",
    "filter_genes",
    "log1p",
    "normalize_per_cell",
    "regress_out",
    "scale",
    "sqrt",
    "subsample",
]

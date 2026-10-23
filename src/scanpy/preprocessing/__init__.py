"""Preprocessing functions."""

from __future__ import annotations

from ..neighbors import neighbors
from ..neighbors._bbknn import bbknn
from ..neighbors._types import NeighborsOut
from ._combat import combat
from ._deprecated.sampling import subsample
from ._harmony import harmony_integrate
from ._hashsolo import HashsoloOut, hashsolo
from ._highly_variable_genes import highly_variable_genes
from ._normalization import NormalizeTotalOut, normalize_shifted_clr, normalize_total
from ._pca import PcaOut, pca
from ._qc import QcMetricsOut, calculate_qc_metrics
from ._recipes import recipe_seurat, recipe_weinreb17, recipe_zheng17
from ._scale import ScaleOut, scale
from ._scrublet import scrublet, scrublet_simulate_doublets
from ._simple import (
    downsample_counts,
    filter_cells,
    filter_genes,
    log1p,
    regress_out,
    sample,
    sqrt,
)

__all__ = [
    "HashsoloOut",
    "NeighborsOut",
    "NormalizeTotalOut",
    "PcaOut",
    "QcMetricsOut",
    "ScaleOut",
    "bbknn",
    "calculate_qc_metrics",
    "combat",
    "downsample_counts",
    "filter_cells",
    "filter_genes",
    "harmony_integrate",
    "hashsolo",
    "highly_variable_genes",
    "log1p",
    "neighbors",
    "normalize_shifted_clr",
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

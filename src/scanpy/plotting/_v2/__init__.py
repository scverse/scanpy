"""Scanpy plots."""

from __future__ import annotations

from ._api import hv_init
from ._core import (
    diffmap,
    dotplot,
    heatmap,
    matrixplot,
    pca,
    scatter,
    stacked_violin,
    tracksplot,
    tsne,
    umap,
    violin,
)
from ._pp import highest_expr_genes, highly_variable_genes, scrublet_score_distribution
from ._tl import draw_graph, embedding_density, ranking

__all__ = [
    "diffmap",
    "dotplot",
    "draw_graph",
    "embedding_density",
    "heatmap",
    "highest_expr_genes",
    "highly_variable_genes",
    "hv_init",
    "matrixplot",
    "pca",
    "ranking",
    "scatter",
    "scrublet_score_distribution",
    "stacked_violin",
    "tracksplot",
    "tsne",
    "umap",
    "violin",
]

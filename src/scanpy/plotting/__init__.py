"""Plotting functions and classes."""

from __future__ import annotations

from . import palettes
from ._anndata import (
    clustermap,
    correlation_matrix,
    dendrogram,
    heatmap,
    ranking,
    scatter,
    tracksplot,
    violin,
)
from ._dotplot import DotPlot, dotplot
from ._matrixplot import MatrixPlot, matrixplot
from ._preprocessing import filter_genes_dispersion, highly_variable_genes
from ._qc import highest_expr_genes
from ._rcmod import set_rcParams_defaults, set_rcParams_scanpy
from ._scrublet import scrublet_score_distribution
from ._stacked_violin import StackedViolin, stacked_violin
from ._tools import (
    dpt_groups_pseudotime,
    dpt_timeseries,
    embedding_density,
    pca_loadings,
    pca_overview,
    pca_scatter,
    pca_variance_ratio,
    rank_genes_groups,
    rank_genes_groups_dotplot,
    rank_genes_groups_heatmap,
    rank_genes_groups_matrixplot,
    rank_genes_groups_stacked_violin,
    rank_genes_groups_tracksplot,
    rank_genes_groups_violin,
    sim,
)
from ._tools.paga import (
    paga,
    paga_adjacency,  # noqa: F401
    paga_compare,
    paga_path,
)
from ._tools.scatterplots import (
    diffmap,
    draw_graph,
    embedding,
    pca,
    spatial,
    tsne,
    umap,
)
from ._utils import matrix, timeseries, timeseries_as_heatmap, timeseries_subplot

__all__ = [
    "DotPlot",
    "MatrixPlot",
    "StackedViolin",
    "clustermap",
    "correlation_matrix",
    "dendrogram",
    "diffmap",
    "dotplot",
    "dpt_groups_pseudotime",
    "dpt_timeseries",
    "draw_graph",
    "embedding",
    "embedding_density",
    "filter_genes_dispersion",
    "heatmap",
    "highest_expr_genes",
    "highly_variable_genes",
    "matrix",
    "matrixplot",
    "paga",
    "paga_compare",
    "paga_path",
    "palettes",
    "pca",
    "pca_loadings",
    "pca_overview",
    "pca_scatter",
    "pca_variance_ratio",
    "rank_genes_groups",
    "rank_genes_groups_dotplot",
    "rank_genes_groups_heatmap",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_tracksplot",
    "rank_genes_groups_violin",
    "ranking",
    "scatter",
    "scrublet_score_distribution",
    "set_rcParams_defaults",
    "set_rcParams_scanpy",
    "sim",
    "spatial",
    "stacked_violin",
    "timeseries",
    "timeseries_as_heatmap",
    "timeseries_subplot",
    "tracksplot",
    "tsne",
    "umap",
    "violin",
]

"""Plotting functions and classes."""

from __future__ import annotations

from scverse_misc import Deprecation, deprecated

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
from ._easter_egg import dogplot
from ._matrixplot import MatrixPlot, matrixplot
from ._preprocessing import highly_variable_genes
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
from ._utils import matrix

__all__ = [
    "DotPlot",
    "MatrixPlot",
    "StackedViolin",
    "clustermap",
    "correlation_matrix",
    "dendrogram",
    "diffmap",
    "dogplot",
    "dotplot",
    "dpt_groups_pseudotime",
    "dpt_timeseries",
    "draw_graph",
    "embedding",
    "embedding_density",
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
    "tracksplot",
    "tsne",
    "umap",
    "violin",
]


@deprecated(Deprecation("1.11.5", "Use :func:`scanpy.pl.dpt_timeseries` instead."))
def timeseries(*args, **kwargs):
    from ._utils import timeseries

    return timeseries(*args, **kwargs)


@deprecated(Deprecation("1.11.5", "Use :func:`scanpy.pl.dpt_timeseries` instead."))
def timeseries_as_heatmap(*args, **kwargs):
    from ._utils import timeseries_as_heatmap

    return timeseries_as_heatmap(*args, **kwargs)


@deprecated(Deprecation("1.11.5", "Use :func:`scanpy.pl.dpt_timeseries` instead."))
def timeseries_subplot(*args, **kwargs):
    from ._utils import timeseries_subplot

    return timeseries_subplot(*args, **kwargs)


@deprecated(Deprecation("1.13.0", "Use :func:`scanpy.pl.pca` instead."))
def pca_scatter(*args, **params):
    return pca(*args, **params)

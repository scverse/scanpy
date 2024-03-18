from __future__ import annotations

from typing import Any

from ._dendrogram import dendrogram
from ._diffmap import diffmap
from ._dpt import dpt
from ._draw_graph import draw_graph
from ._embedding_density import embedding_density
from ._ingest import (
    Ingest,  # noqa: F401
    ingest,
)
from ._leiden import leiden
from ._louvain import louvain
from ._marker_gene_overlap import marker_gene_overlap
from ._paga import (
    paga,
    paga_compare_paths,  # noqa: F401
    paga_degrees,  # noqa: F401
    paga_expression_entropies,  # noqa: F401
)
from ._rank_genes_groups import filter_rank_genes_groups, rank_genes_groups
from ._score_genes import score_genes, score_genes_cell_cycle
from ._sim import sim
from ._tsne import tsne
from ._umap import umap


def __getattr__(name: str) -> Any:
    if name == "pca":
        from ..preprocessing import pca

        return pca
    raise AttributeError(name)


__all__ = [
    "dendrogram",
    "diffmap",
    "dpt",
    "draw_graph",
    "embedding_density",
    "ingest",
    "leiden",
    "louvain",
    "marker_gene_overlap",
    "paga",
    "filter_rank_genes_groups",
    "rank_genes_groups",
    "score_genes",
    "score_genes_cell_cycle",
    "sim",
    "tsne",
    "umap",
]

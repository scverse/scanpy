"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from ._utils import pbmc68k_reduced

if TYPE_CHECKING:
    from anndata import AnnData

# setup variables

adata: AnnData


def setup():
    global adata  # noqa: PLW0603
    adata = pbmc68k_reduced()
    assert "X_pca" in adata.obsm


def time_umap():
    sc.tl.umap(adata)


def peakmem_umap():
    sc.tl.umap(adata)


def time_diffmap():
    sc.tl.diffmap(adata)


def peakmem_diffmap():
    sc.tl.diffmap(adata)


def time_leiden():
    sc.tl.leiden(adata, flavor="igraph")


def peakmem_leiden():
    sc.tl.leiden(adata, flavor="igraph")


def time_rank_genes_groups() -> None:
    sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")


def peakmem_rank_genes_groups() -> None:
    sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")

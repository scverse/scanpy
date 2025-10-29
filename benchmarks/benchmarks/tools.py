"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

import anndata as ad

import scanpy as sc

from ._utils import pbmc68k_reduced


class ToolsSuite:  # noqa: D101
    def setup_cache(self):
        adata = pbmc68k_reduced()
        assert "X_pca" in adata.obsm
        adata.write_h5ad("adata.h5ad")

    def setup(self):
        self.adata = ad.read_h5ad("adata.h5ad")

    def time_umap(self):
        sc.tl.umap(self.adata)

    def peakmem_umap(self):
        sc.tl.umap(self.adata)

    def time_diffmap(self):
        sc.tl.diffmap(self.adata)

    def peakmem_diffmap(self):
        sc.tl.diffmap(self.adata)

    def time_leiden(self):
        sc.tl.leiden(self.adata, flavor="igraph")

    def peakmem_leiden(self):
        sc.tl.leiden(self.adata, flavor="igraph")

    def time_rank_genes_groups(self) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

    def peakmem_rank_genes_groups(self) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

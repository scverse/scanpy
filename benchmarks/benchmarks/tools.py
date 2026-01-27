"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

import anndata as ad

import scanpy as sc

from ._utils import pbmc68k_reduced


class ToolsSuite:  # noqa: D101
    def setup_cache(self) -> None:
        adata = pbmc68k_reduced()
        assert "X_pca" in adata.obsm
        adata.write_h5ad("adata.h5ad")

    def setup(self) -> None:
        self.adata = ad.read_h5ad("adata.h5ad")

    def time_umap(self) -> None:
        sc.tl.umap(self.adata)

    def peakmem_umap(self) -> None:
        sc.tl.umap(self.adata)

    def time_diffmap(self) -> None:
        sc.tl.diffmap(self.adata)

    def peakmem_diffmap(self) -> None:
        sc.tl.diffmap(self.adata)

    def time_leiden(self) -> None:
        sc.tl.leiden(self.adata, flavor="igraph")

    def peakmem_leiden(self) -> None:
        sc.tl.leiden(self.adata, flavor="igraph")

    def time_rank_genes_groups(self) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

    def peakmem_rank_genes_groups(self) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

    def time_score_genes(self) -> None:
        sc.tl.score_genes(self.adata, self.adata.var_names[:500])

    def peakmem_score_genes(self) -> None:
        sc.tl.score_genes(self.adata, self.adata.var_names[:500])

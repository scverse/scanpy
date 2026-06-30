"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from ._utils import pbmc3k, pbmc68k_reduced, to_off_axis

if TYPE_CHECKING:
    import anndata as ad


class ToolsSuite:  # noqa: D101
    def setup_cache(self) -> ad.AnnData:
        adata = pbmc68k_reduced()
        assert "X_pca" in adata.obsm
        return adata

    def setup(self, adata: ad.AnnData) -> None:
        self.adata = adata.copy()

    def time_umap(self, *_) -> None:
        sc.tl.umap(self.adata, rng=None)

    def peakmem_umap(self, *_) -> None:
        sc.tl.umap(self.adata, rng=None)

    def time_diffmap(self, *_) -> None:
        sc.tl.diffmap(self.adata)

    def peakmem_diffmap(self, *_) -> None:
        sc.tl.diffmap(self.adata)

    def time_leiden(self, *_) -> None:
        sc.tl.leiden(self.adata, flavor="igraph")

    def peakmem_leiden(self, *_) -> None:
        sc.tl.leiden(self.adata, flavor="igraph")

    def time_rank_genes_groups(self, *_) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

    def peakmem_rank_genes_groups(self, *_) -> None:
        sc.tl.rank_genes_groups(self.adata, "bulk_labels", method="wilcoxon")

    def time_combat(self, *_) -> None:
        sc.pp.combat(self.adata, key="bulk_labels")

    def peakmem_combat(self, *_) -> None:
        sc.pp.combat(self.adata, key="bulk_labels")


class ScoreGenesSuite:
    """End-to-end benchmark for `sc.tl.score_genes` on sparse data.

    `score_genes` reduces the control- and target-gene blocks with
    `_sparse_nanmean` when `.X` is sparse, so this covers the public path that
    consumes the kernel for both the CSR and the off-axis (CSC) layout.
    """

    params: tuple[str, ...] = ("pbmc3k", "pbmc3k-off-axis")
    param_names = ("layout",)

    def setup_cache(self) -> dict[str, ad.AnnData]:
        adata = pbmc3k()
        adata_orig = adata.copy()
        adata.X = to_off_axis(adata.X)
        return {"pbmc3k": adata_orig, "pbmc3k-off-axis": adata}

    def setup(self, cache: dict[str, ad.AnnData], layout: str) -> None:
        self.adata = cache[layout].copy()
        self.gene_list = self.adata.var_names[:100].tolist()
        # warm up the numba JIT (score_genes -> _sparse_nanmean) so compilation
        # is excluded from the timing
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

    def time_score_genes(self, *_) -> None:
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

    def peakmem_score_genes(self, *_) -> None:
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

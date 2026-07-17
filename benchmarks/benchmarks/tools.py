"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

import anndata as ad

import scanpy as sc

from ._utils import pbmc3k, pbmc68k_reduced, to_off_axis


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

    def time_combat(self) -> None:
        sc.pp.combat(self.adata, key="bulk_labels")

    def peakmem_combat(self) -> None:
        sc.pp.combat(self.adata, key="bulk_labels")


class ScoreGenesSuite:
    """End-to-end benchmark for `sc.tl.score_genes` on sparse data.

    `score_genes` reduces the control- and target-gene blocks with
    `_sparse_nanmean` when `.X` is sparse, so this covers the public path that
    consumes the kernel for both the CSR and the off-axis (CSC) layout.
    """

    params: tuple[str, ...] = ("pbmc3k", "pbmc3k-off-axis")
    param_names = ("layout",)

    def setup_cache(self) -> None:
        adata = pbmc3k()
        adata.write_h5ad("pbmc3k.h5ad")
        adata.X = to_off_axis(adata.X)
        adata.write_h5ad("pbmc3k-off-axis.h5ad")

    def setup(self, layout: str) -> None:
        self.adata = ad.read_h5ad(f"{layout}.h5ad")
        self.gene_list = self.adata.var_names[:100].tolist()
        # warm up the numba JIT (score_genes -> _sparse_nanmean) so compilation
        # is excluded from the timing
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

    def time_score_genes(self, *_) -> None:
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

    def peakmem_score_genes(self, *_) -> None:
        sc.tl.score_genes(self.adata, self.gene_list, rng=0)

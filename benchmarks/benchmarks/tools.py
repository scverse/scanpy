"""Benchmark tool operations in Scanpy.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/tools.html>.
"""

from __future__ import annotations

import anndata as ad

import scanpy as sc

from ._utils import pbmc3k, pbmc68k_reduced


class ToolsSuite:  # noqa: D101
    def setup_cache(self) -> None:
        adata = pbmc68k_reduced()
        assert "X_pca" in adata.obsm
        adata.write_h5ad("adata.h5ad")

    def setup(self) -> None:
        self.adata = ad.read_h5ad("adata.h5ad")

    def time_umap(self) -> None:
        sc.tl.umap(self.adata, rng=None)

    def peakmem_umap(self) -> None:
        sc.tl.umap(self.adata, rng=None)

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


class CombatSuite:
    """Benchmark combat batch correction."""

    def setup_cache(self) -> None:
        import numpy as np

        adata = pbmc3k()
        sc.pp.highly_variable_genes(adata, n_top_genes=500)
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata, max_value=10)
        # assign cells to 3 batches deterministically
        np.random.seed(0)
        adata.obs["batch"] = np.random.choice(["A", "B", "C"], size=adata.n_obs)
        adata.write_h5ad("adata_combat.h5ad")

    def setup(self) -> None:
        self.adata = ad.read_h5ad("adata_combat.h5ad")

    def time_combat(self) -> None:
        sc.pp.combat(self.adata, key="batch")

    def peakmem_combat(self) -> None:
        sc.pp.combat(self.adata, key="batch")

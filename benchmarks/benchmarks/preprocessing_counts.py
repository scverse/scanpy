"""Benchmark preprocessing operations in Scanpy that run on counts.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/preprocessing.html>.
"""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING

import anndata as ad

import scanpy as sc

from ._utils import get_count_dataset

if TYPE_CHECKING:
    from ._utils import Dataset, KeyCount


# ASV suite
class PreprocessingCountsSuite:  # noqa: D101
    params: tuple[list[Dataset], list[KeyCount]] = (
        ["pbmc68k_reduced", "pbmc3k"],
        ["counts", "counts-off-axis"],
    )
    param_names = ("dataset", "layer")

    def setup_cache(self) -> None:
        """Without this caching, asv was running several processes which meant the data was repeatedly downloaded."""
        for dataset, layer in product(*self.params):
            adata, batch_key = get_count_dataset(dataset, layer=layer)
            assert "lop1p" not in adata.uns
            adata.uns["batch_key"] = batch_key
            adata.write_h5ad(f"{dataset}_{layer}.h5ad")

    def setup(self, dataset, layer) -> None:
        self.adata = ad.read_h5ad(f"{dataset}_{layer}.h5ad")

    def time_filter_cells(self, *_) -> None:
        sc.pp.filter_cells(self.adata, min_genes=100)

    def peakmem_filter_cells(self, *_) -> None:
        sc.pp.filter_cells(self.adata, min_genes=100)

    def time_filter_genes(self, *_) -> None:
        sc.pp.filter_genes(self.adata, min_cells=3)

    def peakmem_filter_genes(self, *_) -> None:
        sc.pp.filter_genes(self.adata, min_cells=3)

    def time_scrublet(self, *_) -> None:
        sc.pp.scrublet(self.adata, batch_key=self.adata.uns["batch_key"])

    def peakmem_scrublet(self, *_) -> None:
        sc.pp.scrublet(self.adata, batch_key=self.adata.uns["batch_key"])


# Can’t do seurat v3 yet: https://github.com/conda-forge/scikit-misc-feedstock/issues/17
"""
def time_hvg_seurat_v3(self):
    # seurat v3 runs on counts
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3_paper")


def peakmem_hvg_seurat_v3(self):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3_paper")
"""


class FastSuite:
    """Suite for fast preprocessing operations."""

    params: tuple[list[Dataset], list[KeyCount]] = (
        ["pbmc3k", "pbmc68k_reduced", "bmmc", "lung93k"],
        ["counts", "counts-off-axis"],
    )
    param_names = ("dataset", "layer")

    def setup_cache(self) -> None:
        """Without this caching, asv was running several processes which meant the data was repeatedly downloaded."""
        for dataset, layer in product(*self.params):
            adata, _ = get_count_dataset(dataset, layer=layer)
            assert "lop1p" not in adata.uns
            adata.write_h5ad(f"{dataset}_{layer}.h5ad")

    def setup(self, dataset, layer) -> None:
        self.adata = ad.read_h5ad(f"{dataset}_{layer}.h5ad")

    def time_calculate_qc_metrics(self, *_) -> None:
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def peakmem_calculate_qc_metrics(self, *_) -> None:
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def time_normalize_total(self, *_) -> None:
        sc.pp.normalize_total(self.adata, target_sum=1e4)

    def peakmem_normalize_total(self, *_) -> None:
        sc.pp.normalize_total(self.adata, target_sum=1e4)

    def time_log1p(self, *_) -> None:
        # TODO: This would fail: assert "log1p" not in self.adata.uns, "ASV bug?"
        # https://github.com/scverse/scanpy/issues/3052
        self.adata.uns.pop("log1p", None)
        sc.pp.log1p(self.adata)

    def peakmem_log1p(self, *_) -> None:
        self.adata.uns.pop("log1p", None)
        sc.pp.log1p(self.adata)

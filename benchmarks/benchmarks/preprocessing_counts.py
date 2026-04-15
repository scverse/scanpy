"""Benchmark preprocessing operations in Scanpy that run on counts.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/preprocessing.html>.
"""

from __future__ import annotations

from inspect import signature
from itertools import product
from typing import TYPE_CHECKING

import anndata as ad

import scanpy as sc
from scanpy._utils import get_literal_vals
from scanpy.get._aggregated import AggType

from ._utils import get_count_dataset, get_dataset

if TYPE_CHECKING:
    from typing import Any

    from ._utils import Dataset, KeyCount


def cache_adata(dataset: Dataset, layer: KeyCount) -> None:
    """Without this caching, asv was running several processes which meant the data was repeatedly downloaded."""
    adata, batch_key = get_count_dataset(dataset, layer=layer)
    assert "lop1p" not in adata.uns
    adata.uns["batch_key"] = batch_key
    adata.write_h5ad(f"{dataset}_{layer}.h5ad")


class PreprocessingCountsSuite:  # noqa: D101
    params: tuple[list[Dataset], list[KeyCount]] = (
        ["pbmc68k_reduced", "pbmc3k"],
        ["counts", "counts-off-axis"],
    )
    param_names = ("dataset", "layer")

    def setup_cache(self) -> None:
        for dataset, layer in product(*self.params):
            cache_adata(dataset, layer)

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

    # sciki-misc does not exit on osx-arm64
    # https://github.com/conda-forge/scikit-misc-feedstock/pull/29
    # def time_hvg_seurat_v3(self, *_):
    #     # seurat v3 runs on counts
    #     sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3_paper")

    # def peakmem_hvg_seurat_v3(self, *_):
    #     sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3_paper")


class PreprocessingCountsRngSuite:  # noqa: D101
    params: tuple[list[Dataset], list[str], list[str]] = (
        ["pbmc68k_reduced", "pbmc3k"],
        ["rng", "random_state"],
    )
    param_names = ("dataset", "layer")

    def setup_cache(self) -> None:
        for dataset in self.params[0]:
            cache_adata(dataset, "counts")

    def setup(self, dataset, rng_arg) -> None:
        if (
            rng_arg == "rng"
            and "rng" not in signature(sc.pp.downsample_counts).parameters
        ):
            raise NotImplementedError
        self.adata = ad.read_h5ad(f"{dataset}_counts.h5ad")
        self.rng_kw: Any = {rng_arg: 0}
        self.total = self.adata.X.sum() / 10

    def time_downsample_per_cell(self, *_) -> None:
        sc.pp.downsample_counts(self.adata, counts_per_cell=3, **self.rng_kw)

    def peakmem_downsample_per_cell(self, *_) -> None:
        sc.pp.downsample_counts(self.adata, counts_per_cell=3, **self.rng_kw)

    def time_downsample_total(self, *_) -> None:
        sc.pp.downsample_counts(self.adata, total_counts=self.total, **self.rng_kw)

    def peakmem_downsample_total(self, *_) -> None:
        sc.pp.downsample_counts(self.adata, total_counts=self.total, **self.rng_kw)


class FastSuite:
    """Suite for fast preprocessing operations."""

    params: tuple[list[Dataset], list[KeyCount]] = (
        ["pbmc3k", "pbmc68k_reduced", "bmmc", "lung93k"],
        ["counts", "counts-off-axis"],
    )
    param_names = ("dataset", "layer")

    def setup_cache(self) -> None:
        for dataset, layer in product(*self.params):
            cache_adata(dataset, layer)

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


class Agg:  # noqa: D101
    params: tuple[AggType] = tuple(get_literal_vals(AggType))
    param_names = ("agg_name",)

    def setup_cache(self) -> None:
        """Without this caching, asv was running several processes which meant the data was repeatedly downloaded."""
        adata, _ = get_dataset("lung93k")
        adata.write_h5ad("lung93k.h5ad")

    def setup(self, agg_name: AggType) -> None:
        self.adata = ad.read_h5ad("lung93k.h5ad")
        self.agg_name = agg_name

    def time_agg(self, *_) -> None:
        sc.get.aggregate(self.adata, by="PatientNumber", func=self.agg_name)

    def peakmem_agg(self, *_) -> None:
        sc.get.aggregate(self.adata, by="PatientNumber", func=self.agg_name)

"""Benchmark preprocessing operations in Scanpy that run on counts.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/preprocessing.html>.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from ._utils import get_count_dataset

if TYPE_CHECKING:
    from anndata import AnnData

    from ._utils import Dataset, KeyCount

# setup variables

adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, layer: KeyCount, *_):
    """Set up global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_count_dataset(dataset, layer=layer)
    assert "log1p" not in adata.uns


# ASV suite

params: tuple[list[Dataset], list[KeyCount]] = (
    ["pbmc68k_reduced", "pbmc3k"],
    ["counts", "counts-off-axis"],
)
param_names = ["dataset", "layer"]


def time_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def peakmem_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def time_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


def peakmem_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


def time_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


def peakmem_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


# Can’t do seurat v3 yet: https://github.com/conda-forge/scikit-misc-feedstock/issues/17
"""
def time_hvg_seurat_v3(*_):
    # seurat v3 runs on counts
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3_paper")


def peakmem_hvg_seurat_v3(*_):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3_paper")
"""


class FastSuite:
    """Suite for fast preprocessing operations."""

    params: tuple[list[Dataset], list[KeyCount]] = (
        ["pbmc3k", "pbmc68k_reduced", "bmmc", "lung93k"],
        ["counts", "counts-off-axis"],
    )
    param_names = ("dataset", "layer")

    def time_calculate_qc_metrics(self, *_):
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def peakmem_calculate_qc_metrics(self, *_):
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def time_normalize_total(self, *_):
        sc.pp.normalize_total(adata, target_sum=1e4)

    def peakmem_normalize_total(self, *_):
        sc.pp.normalize_total(adata, target_sum=1e4)

    def time_log1p(self, *_):
        # TODO: This would fail: assert "log1p" not in adata.uns, "ASV bug?"
        # https://github.com/scverse/scanpy/issues/3052
        adata.uns.pop("log1p", None)
        sc.pp.log1p(adata)

    def peakmem_log1p(self, *_):
        adata.uns.pop("log1p", None)
        sc.pp.log1p(adata)

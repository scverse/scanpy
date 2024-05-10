"""
This module will benchmark preprocessing operations in Scanpy
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from asv_runner.benchmarks.mark import skip_for_params

import scanpy as sc
from scanpy.preprocessing._utils import _get_mean_var

from .utils import get_dataset

if TYPE_CHECKING:
    from anndata import AnnData

    from .utils import Dataset


adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, *_):
    """Setup global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_dataset(dataset)


# The actual test suite begins here

params: list[Dataset] = ["pbmc68k_reduced", "pbmc3k"]
param_names = ["dataset"]


def time_calculate_qc_metrics(*_):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def peakmem_calculate_qc_metrics(*_):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def time_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def peakmem_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def time_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


def peakmem_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


# scublet doesn’t work with these datasets
@skip_for_params([("pbmc68k_reduced",)])
def time_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


@skip_for_params([("pbmc68k_reduced",)])
def peakmem_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


def time_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def peakmem_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def time_log1p(*_):
    sc.pp.log1p(adata, layer="counts" if "counts" in adata.layers else None)


def peakmem_time_log1p(*_):
    sc.pp.log1p(adata, layer="counts" if "counts" in adata.layers else None)


def time_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def peakmem_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def time_highly_variable_genes(*_):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def peakmem_highly_variable_genes(*_):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# regress_out is very slow for this dataset
@skip_for_params([("pbmc3k",)])
def time_regress_out(*_):
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


@skip_for_params([("pbmc3k",)])
def peakmem_regress_out(*_):
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


def time_scale(*_):
    sc.pp.scale(adata, max_value=10)


def peakmem_scale(*_):
    sc.pp.scale(adata, max_value=10)


class FastSuite:
    """Suite for fast preprocessing operations."""

    params: list[Dataset] = ["pbmc3k", "pbmc68k_reduced", "bmmc", "lung93k"]
    param_names = ["dataset"]

    def time_mean_var(self, *_):
        _get_mean_var(adata.X)

    def peakmem_mean_var(self, *_):
        _get_mean_var(adata.X)

"""
This module will benchmark preprocessing operations in Scanpy that run on log-transformed data
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from asv_runner.benchmarks.mark import skip_for_params

import scanpy as sc
from scanpy.preprocessing._utils import _get_mean_var

from ._utils import get_dataset

if TYPE_CHECKING:
    from anndata import AnnData

    from ._utils import Dataset

# setup variables


adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, *_):
    """Setup global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_dataset(dataset)


# ASV suite

params: list[Dataset] = ["pbmc68k_reduced", "pbmc3k"]
param_names = ["dataset"]


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

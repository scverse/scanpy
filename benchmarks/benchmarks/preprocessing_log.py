"""Benchmark preprocessing operations in Scanpy that run on log-transformed data.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/preprocessing.html>.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc
from scanpy.preprocessing._utils import _get_mean_var

from ._utils import get_dataset, param_skipper

if TYPE_CHECKING:
    from anndata import AnnData

    from ._utils import Dataset, KeyX

# setup variables


adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, layer: KeyX, *_):
    """Set up global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_dataset(dataset, layer=layer)


# ASV suite

params: tuple[list[Dataset], list[KeyX]] = (
    ["pbmc68k_reduced", "pbmc3k"],
    [None, "off-axis"],
)
param_names = ["dataset", "layer"]

skip_when = param_skipper(param_names, params)


def time_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def peakmem_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def time_highly_variable_genes(*_):
    # the default flavor runs on log-transformed data
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def peakmem_highly_variable_genes(*_):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# regress_out is very slow for this dataset
@skip_when(dataset={"pbmc3k"})
def time_regress_out(*_):
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


@skip_when(dataset={"pbmc3k"})
def peakmem_regress_out(*_):
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


def time_scale(*_):
    sc.pp.scale(adata, max_value=10)


def peakmem_scale(*_):
    sc.pp.scale(adata, max_value=10)


class FastSuite:
    """Suite for fast preprocessing operations."""

    params: tuple[list[Dataset], list[KeyX]] = (
        ["pbmc3k", "pbmc68k_reduced", "bmmc", "lung93k"],
        [None, "off-axis"],
    )
    param_names = ("dataset", "layer")

    def time_mean_var(self, *_):
        _get_mean_var(adata.X)

    def peakmem_mean_var(self, *_):
        _get_mean_var(adata.X)

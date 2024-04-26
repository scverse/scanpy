"""
This module will benchmark preprocessing operations in Scanpy
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from asv_runner.benchmarks.mark import skip_for_params

import scanpy as sc
from scanpy.preprocessing._utils import _get_mean_var

from .utils import bmmc8k, lung93k, pbmc3k, pbmc68k_reduced

if TYPE_CHECKING:
    from anndata import AnnData


adata: AnnData

params = ["pbmc68k_reduced", "bmmc8k"]
param_names = ["dataset"]


def setup(*params: str):
    """Setup all tests.

    The SparseDenseSuite below defines a parameter,
    the other tests none.
    """
    global adata

    if params[0] == "pbmc3k":
        adata = pbmc3k()
    elif params[0] == "pbmc68k_reduced":
        adata = pbmc68k_reduced()
    elif params[0] == "bmmc8k":
        adata = bmmc8k()
    elif params[0] == "lung93k":
        adata = lung93k()
    else:
        raise ValueError(f"Unknown dataset {params[0]}")


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
@skip_for_params([("pbmc3k",), ("pbmc68k_reduced",)])
def time_scrublet(*_):
    sc.pp.scrublet(adata, batch_key="sample")


@skip_for_params([("pbmc3k",), ("pbmc68k_reduced",)])
def peakmem_scrublet(*_):
    sc.pp.scrublet(adata, batch_key="sample")


def time_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def peakmem_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def time_log1p(*_):
    sc.pp.log1p(adata)


def peakmem_time_log1p(*_):
    sc.pp.log1p(adata)


def time_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def peakmem_pca(*_):
    sc.pp.pca(adata, svd_solver="arpack")


def time_highly_variable_genes(*_):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def peakmem_highly_variable_genes(*_):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def time_regress_out(*_):
    sc.pp.regress_out(adata, ["n_counts", "percent_mito"])


def peakmem_regress_out(*_):
    sc.pp.regress_out(adata, ["n_counts", "percent_mito"])


def time_scale(*_):
    sc.pp.scale(adata, max_value=10)


def peakmem_scale(*_):
    sc.pp.scale(adata, max_value=10)


class SparseDenseSuite:
    params = ["pbmc3k", "pbmc68k_reduced", "bmmc8k", "lung93k"]
    param_names = ["dataset"]

    def time_mean_var(self, *_):
        _get_mean_var(adata.X)

    def peakmem_mean_var(self, *_):
        _get_mean_var(adata.X)

"""
This module will benchmark preprocessing operations in Scanpy
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc
from scanpy.preprocessing._utils import _get_mean_var

from .utils import lung93k, pbmc3k, pbmc68k_reduced

if TYPE_CHECKING:
    from anndata import AnnData


adata: AnnData


def setup(*params: str):
    """Setup all tests.

    The SparseDenseSuite below defines a parameter,
    the other tests none.
    """
    global adata

    if len(params) == 0 or params[0] == "pbmc68k_reduced":
        adata = pbmc68k_reduced()
    elif params[0] == "pbmc3k":
        adata = pbmc3k()
    elif params[0] == "lung93k":
        adata = lung93k()
    else:
        raise ValueError(f"Unknown dataset {params[0]}")


def time_calculate_qc_metrics():
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def peakmem_calculate_qc_metrics():
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def time_filter_cells():
    sc.pp.filter_cells(adata, min_genes=200)


def peakmem_filter_cells():
    sc.pp.filter_cells(adata, min_genes=200)


def time_filter_genes():
    sc.pp.filter_genes(adata, min_cells=3)


def peakmem_filter_genes():
    sc.pp.filter_genes(adata, min_cells=3)


def time_normalize_total():
    sc.pp.normalize_total(adata, target_sum=1e4)


def peakmem_normalize_total():
    sc.pp.normalize_total(adata, target_sum=1e4)


def time_log1p():
    sc.pp.log1p(adata)


def peakmem_time_log1p():
    sc.pp.log1p(adata)


def time_pca():
    sc.pp.pca(adata, svd_solver="arpack")


def peakmem_pca():
    sc.pp.pca(adata, svd_solver="arpack")


def time_highly_variable_genes():
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def peakmem_highly_variable_genes():
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


def time_regress_out():
    sc.pp.regress_out(adata, ["n_counts", "percent_mito"])


def peakmem_regress_out():
    sc.pp.regress_out(adata, ["n_counts", "percent_mito"])


def time_scale():
    sc.pp.scale(adata, max_value=10)


def peakmem_scale():
    sc.pp.scale(adata, max_value=10)


class SparseDenseSuite:
    params = ["pbmc68k_reduced", "pbmc3k", "lung93k"]
    param_names = ["dataset"]

    def time_mean_var(self, *_):
        _get_mean_var(adata.X)

    def peakmem_mean_var(self, *_):
        _get_mean_var(adata.X)

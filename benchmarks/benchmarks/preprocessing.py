"""
This module will benchmark preprocessing operations in Scanpy
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from .utils import pbmc68k_reduced

if TYPE_CHECKING:
    from anndata import AnnData


adata: AnnData


def setup():
    global adata
    adata = pbmc68k_reduced()


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

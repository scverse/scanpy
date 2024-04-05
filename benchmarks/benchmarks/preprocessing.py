"""
This module will benchmark preprocessing operations in ScanPy
API documentation: https://scanpy.readthedocs.io/en/stable/api.html#module-scanpy.pp


"""

from __future__ import annotations

import scanpy as sc


class PreprocessingSuite:
    _data_dict = dict(pbmc68k_reduced=sc.datasets.pbmc68k_reduced())
    params = _data_dict.keys()
    param_names = ["input_data"]

    def setup(self, input_data):
        self.adata = self._data_dict[input_data]

    def time_calculate_qc_metrics(self, input_data):
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def peakmem_calculate_qc_metrics(self, input_data):
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

    def time_filter_cells(self, input_data):
        sc.pp.filter_cells(self.adata, min_genes=200)

    def peakmem_filter_cells(self, input_data):
        sc.pp.filter_cells(self.adata, min_genes=200)

    def time_filter_genes(self, input_data):
        sc.pp.filter_genes(self.adata, min_cells=3)

    def peakmem_filter_genes(self, input_data):
        sc.pp.filter_genes(self.adata, min_cells=3)

    def time_normalize_total(self, input_data):
        sc.pp.normalize_total(self.adata, target_sum=1e4)

    def peakmem_normalize_total(self, input_data):
        sc.pp.normalize_total(self.adata, target_sum=1e4)

    def time_log1p(self, input_data):
        sc.pp.log1p(self.adata)

    def peakmem_time_log1p(self, input_data):
        sc.pp.log1p(self.adata)

    def time_pca(self, input_data):
        sc.pp.pca(self.adata, svd_solver="arpack")

    def peakmem_pca(self, input_data):
        sc.pp.pca(self.adata, svd_solver="arpack")

    def time_highly_variable_genes(self, input_data):
        sc.pp.highly_variable_genes(
            self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )

    def peakmem_highly_variable_genes(self, input_data):
        sc.pp.highly_variable_genes(
            self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )

    def time_regress_out(self, input_data):
        sc.pp.regress_out(self.adata, ["n_counts", "percent_mito"])

    def peakmem_regress_out(self, input_data):
        sc.pp.regress_out(self.adata, ["n_counts", "percent_mito"])

    def time_scale(self, input_data):
        sc.pp.scale(self.adata, max_value=10)

    def peakmem_scale(self, input_data):
        sc.pp.scale(self.adata, max_value=10)

"""Benchmark preprocessing operations in Scanpy that run on log-transformed data.

API documentation: <https://scanpy.readthedocs.io/en/stable/api/preprocessing.html>.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import anndata as ad

import scanpy as sc

from ._utils import get_dataset, param_skipper

if TYPE_CHECKING:
    from ._utils import Dataset, KeyX


# ASV suite

params: tuple[list[Dataset], list[KeyX]] = (
    ["pbmc68k_reduced", "pbmc3k"],
    [None, "off-axis"],
)
param_names = ("dataset", "layer")
skip_when = param_skipper(param_names, params)


class PreprocessingSuite:  # noqa: D101
    params = params
    param_names = param_names

    def setup_cache(self):
        """Without this caching, asv was running several processes which meant the data was repeatedly downloaded."""
        for dataset in self.params[0]:
            for layer in self.params[1]:
                adata, _ = get_dataset(dataset, layer=layer)
                adata.write_h5ad(f"{dataset}_{layer}.h5ad")

    def setup(self, dataset, layer):
        self.adata = ad.read_h5ad(f"{dataset}_{layer}.h5ad")

    def time_pca(self, *_):
        sc.pp.pca(self.adata, svd_solver="arpack")

    def peakmem_pca(self, *_):
        sc.pp.pca(self.adata, svd_solver="arpack")

    def time_highly_variable_genes(self, *_):
        # the default flavor runs on log-transformed data
        sc.pp.highly_variable_genes(
            self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )

    def peakmem_highly_variable_genes(self, *_):
        sc.pp.highly_variable_genes(
            self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )

    # regress_out is very slow for this dataset
    @skip_when(dataset={"pbmc3k"})
    def time_regress_out(self, *_):
        sc.pp.regress_out(self.adata, ["total_counts", "pct_counts_mt"])

    @skip_when(dataset={"pbmc3k"})
    def peakmem_regress_out(self, *_):
        sc.pp.regress_out(self.adata, ["total_counts", "pct_counts_mt"])

    def time_scale(self, *_):
        sc.pp.scale(self.adata, max_value=10)

    def peakmem_scale(self, *_):
        sc.pp.scale(self.adata, max_value=10)

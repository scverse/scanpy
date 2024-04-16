"""
This module will benchmark tool operations in Scanpy
API documentation: https://scanpy.readthedocs.io/en/stable/api/tools.html
"""

from __future__ import annotations

import scanpy as sc


class ToolsSuite:
    _data_dict = dict(pbmc68k_reduced=sc.datasets.pbmc68k_reduced())
    params = _data_dict.keys()
    param_names = ["input_data"]

    def setup(self, input_data):
        self.adata = self._data_dict[input_data].copy()
        assert "X_pca" in self.adata.obsm

    def time_umap(self, *_):
        sc.tl.umap(self.adata)

    def peakmem_umap(self, *_):
        sc.tl.umap(self.adata)

    def time_diffmap(self, *_):
        sc.tl.diffmap(self.adata)

    def peakmem_diffmap(self, *_):
        sc.tl.diffmap(self.adata)

    def time_leiden(self, *_):
        sc.tl.leiden(self.adata, flavor="igraph")

    def peakmem_leiden(self, *_):
        sc.tl.leiden(self.adata, flavor="igraph")

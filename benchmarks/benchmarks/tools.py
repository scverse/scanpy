"""
This module will benchmark preprocessing operations in ScanPy
API documentation: https://scanpy.readthedocs.io/en/stable/api.html#module-scanpy.pp


"""

from __future__ import annotations

import scanpy as sc


class ToolsSuite:
    _data_dict = dict(pbmc68k_reduced=sc.datasets.pbmc68k_reduced())
    params = _data_dict.keys()
    param_names = ["input_data"]

    def setup(self, input_data):
        self.adata = self._data_dict[input_data]

    def time_pca(self, input_data):
        sc.tl.pca(self.adata)

    def peakmem_pca(self, input_data):
        sc.tl.pca(self.adata)

    def time_umap(self, input_data):
        sc.tl.umap(self.adata)

    def peakmem_umap(self, input_data):
        sc.tl.umap(self.adata)

    def time_diffmap(self, input_data):
        sc.tl.diffmap(self.adata)

    def peakmem_diffmap(self, input_data):
        sc.tl.diffmap(self.adata)

    # def time_embedding_density(self, input_data):
    #     sc.tl.embedding_density(self.adata, groupby="batch")

    # def peakmem_embedding_density(self, input_data):
    #     sc.tl.embedding_density(self.adata, groupby="batch")

    # def time_louvain(self, input_data):
    #     sc.tl.louvain(self.adata)

    # def peakmem_louvain(self, input_data):
    #     sc.tl.louvain(self.adata)

    # def time_paga(self, input_data):
    #     sc.tl.paga(self.adata, groups='louvain')

    # def peakmem_paga(self, input_data):
    #     sc.tl.paga(self.adata, groups='louvain')

    # def time_leiden(self, input_data):
    #     sc.tl.leiden(self.adata)

    # def peakmem_leiden(self, input_data):
    #     sc.tl.leiden(self.adata)

    # def time_draw_graph(self, input_data):
    #     sc.tl.draw_graph(self.adata)

    # def peakmem_draw_graph(self, input_data):
    #     sc.tl.draw_graph(self.adata)

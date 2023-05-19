"""
This module will benchmark preprocessing operations in ScanPy
API documentation: https://scanpy.readthedocs.io/en/stable/api.html#module-scanpy.pp


"""
import tempfile
from pathlib import Path
import sys

from memory_profiler import memory_usage
import numpy as np
import scanpy as sc

from .utils import get_anndata_memsize, sedate, get_peak_mem, get_actualsize


class ToolsSuite:
    params = [sc.datasets.pbmc3k_processed()]
    param_names = ['input_data']

    def setup(self, input_data):
        self.adata = input_data
        
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
    
    

    

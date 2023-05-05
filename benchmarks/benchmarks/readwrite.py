"""
This module will benchmark io of AnnData objects

Things to test:

* Read time, write time
* Peak memory during io
* File sizes

Parameterized by:

* What method is being used
* What data is being included
* Size of data being used

Also interesting:

* io for views
* io for backed objects
* Reading dense as sparse, writing sparse as dense
"""
import tempfile
from pathlib import Path
import sys

from memory_profiler import memory_usage
import numpy as np
import pooch
import scanpy as sc

from .utils import get_anndata_memsize, sedate, get_peak_mem, get_actualsize


PBMC_3K_URL = "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"


class H5ADInMemorySizeSuite:
    params = [PBMC_3K_URL]
    param_names = ["input_url"]

    def setup(self, input_url):
        self.filepath = pooch.retrieve(url=input_url, known_hash=None)
        
    def track_in_memory_size(self, input_url):
        adata = sc.read_10x_mtx(self.filepath)
        adata_size = sys.getsizeof(adata)

        return adata_size

    def track_actual_in_memory_size(self, input_url):
        adata = sc.read_10x_mtx(self.filepath)
        adata_size = get_actualsize(adata)

        return adata_size


class H5ADReadSuite:
    # params = [PBMC_REDUCED_PATH, PBMC_3K_PATH, BM_43K_CSR_PATH]
    params = [PBMC_3K_URL]
    param_names = ["input_url"]

    def setup(self, input_url):
        self.filepath = pooch.retrieve(url=input_url, known_hash=None)

    def time_read_full(self, input_url):
        sc.read_10x_mtx(self.filepath)

    def peakmem_read_full(self, input_url):
        sc.read_10x_mtx(self.filepath)

    def mem_readfull_object(self, input_url):
        return sc.read_10x_mtx(self.filepath)

    def track_read_full_memratio(self, input_url):
        mem_recording = memory_usage(
            (sedate(sc.read_10x_mtx, 0.005), (self.filepath,)), interval=0.001
        )

        base_size = mem_recording[-1] - mem_recording[0]
        print(np.max(mem_recording) - np.min(mem_recording))
        print(base_size)
        return (np.max(mem_recording) - np.min(mem_recording)) / base_size

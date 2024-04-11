"""
This module will benchmark io of Scanpy readwrite operations

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

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING

import anndata
import numpy as np
from memory_profiler import memory_usage

import scanpy as sc

from .utils import get_actualsize, sedate

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path


@dataclass
class Dataset:
    path: Path
    get: Callable[[], anndata.AnnData]


pbmc3k = Dataset(
    path=sc.settings.datasetdir / "pbmc3k_raw.h5ad", get=sc.datasets.pbmc3k
)


class H5ADInMemorySizeSuite:
    _data_dict = dict(pbmc3k=pbmc3k)
    params = _data_dict.keys()
    param_names = ["input_data"]

    def setup(self, input_data: str):
        self.path = self._data_dict[input_data].path
        self.data = self._data_dict[input_data].get()

    def track_in_memory_size(self, *_):
        adata = anndata.read_h5ad(self.path)
        adata_size = sys.getsizeof(adata)

        return adata_size

    def track_actual_in_memory_size(self, *_):
        adata = anndata.read_h5ad(self.path)
        adata_size = get_actualsize(adata)

        return adata_size


class H5ADReadSuite:
    _data_dict = dict(pbmc3k=pbmc3k)
    params = _data_dict.keys()
    param_names = ["input_data"]

    def setup(self, input_data: str):
        self.path = self._data_dict[input_data].path
        self.data = self._data_dict[input_data].get()

    def time_read_full(self, *_):
        anndata.read_h5ad(self.path)

    def peakmem_read_full(self, *_):
        anndata.read_h5ad(self.path)

    def mem_read_full(self, *_):
        anndata.read_h5ad(self.path)

    def track_read_full_memratio(self, *_):
        mem_recording = memory_usage(
            (sedate(anndata.read_h5ad, 0.005), (self.path,)), interval=0.001
        )

        base_size = mem_recording[-1] - mem_recording[0]
        print(np.max(mem_recording) - np.min(mem_recording))
        print(base_size)
        return (np.max(mem_recording) - np.min(mem_recording)) / base_size

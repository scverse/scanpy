"""
Functions returning copies of datasets as cheaply as possible,
i.e. without having to hit the disk or (in case of ``_pbmc3k_normalized``) recomputing normalization.
"""

from __future__ import annotations

import warnings

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    def cache(func):
        return lru_cache(maxsize=None)(func)


from typing import TYPE_CHECKING

import dask.array as da
from dask import delayed
from scipy import sparse

import scanpy as sc

if TYPE_CHECKING:
    from anndata import AnnData
    from anndata._core.sparse_dataset import SparseDataset
# Functions returning the same objects (easy to misuse)


_pbmc3k = cache(sc.datasets.pbmc3k)
_pbmc3k_processed = cache(sc.datasets.pbmc3k_processed)
_pbmc68k_reduced = cache(sc.datasets.pbmc68k_reduced)
_krumsiek11 = cache(sc.datasets.krumsiek11)
_paul15 = cache(sc.datasets.paul15)


# Functions returning copies


def pbmc3k() -> AnnData:
    return _pbmc3k().copy()


def pbmc3k_processed() -> AnnData:
    return _pbmc3k_processed().copy()


def pbmc68k_reduced() -> AnnData:
    return _pbmc68k_reduced().copy()


def krumsiek11() -> AnnData:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "Observation names are not unique", module="anndata"
        )
        return _krumsiek11().copy()


def paul15() -> AnnData:
    return _paul15().copy()


# Derived datasets


@cache
def _pbmc3k_normalized() -> AnnData:
    pbmc = pbmc3k()
    pbmc.X = pbmc.X.astype("float64")  # For better accuracy
    sc.pp.filter_genes(pbmc, min_counts=1)
    sc.pp.log1p(pbmc)
    sc.pp.normalize_total(pbmc)
    sc.pp.highly_variable_genes(pbmc)
    return pbmc


def pbmc3k_normalized() -> AnnData:
    return _pbmc3k_normalized().copy()


class CSRCallable:
    """Dummy class to bypass dask checks"""

    def __new__(cls, shape, dtype):
        return csr_callable(shape, dtype)


def csr_callable(shape: tuple[int, int], dtype) -> sparse.csr_matrix:
    if len(shape) == 0:
        shape = (0, 0)
    if len(shape) == 1:
        shape = (shape[0], 0)
    elif len(shape) == 2:
        pass
    else:
        raise ValueError(shape)

    return sparse.csr_matrix(shape, dtype=dtype)


def make_dask_chunk(x: SparseDataset, start: int, end: int) -> da.Array:
    def take_slice(x, idx):
        return x[idx]

    return da.from_delayed(
        delayed(take_slice)(x, slice(start, end)),
        dtype=x.dtype,
        shape=(end - start, x.shape[1]),
        meta=CSRCallable,
    )


def sparse_dataset_as_dask(x: SparseDataset, stride: int):
    n_chunks, rem = divmod(x.shape[0], stride)

    chunks = []
    cur_pos = 0
    for i in range(n_chunks):
        chunks.append(make_dask_chunk(x, cur_pos, cur_pos + stride))
        cur_pos += stride
    if rem:
        chunks.append(make_dask_chunk(x, cur_pos, x.shape[0]))

    return da.concatenate(chunks, axis=0)

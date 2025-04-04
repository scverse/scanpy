from __future__ import annotations

from typing import TYPE_CHECKING

import dask.array as dd
from dask.distributed import Client, LocalCluster

import scanpy as sc

from ._utils import get_count_dataset

if TYPE_CHECKING:
    from anndata import AnnData

    from ._utils import Dataset, KeyCount

# Setup global variables
adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, layer: KeyCount, *_):
    """Setup global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_count_dataset(dataset, layer=layer)
    assert "log1p" not in adata.uns


# def setup_dask_cluster():
#    """Set up a local Dask cluster for benchmarking."""
#    cluster = LocalCluster(n_workers=4, threads_per_worker=2)
#    client = Client(cluster)
#    return client


def setup_dask_cluster():
    """Set up a local Dask cluster for benchmarking."""
    cluster = LocalCluster(
        n_workers=5, threads_per_worker=2, memory_limit="60GB", timeout="1200s"
    )
    client = Client(cluster)
    return client


# ASV suite
params: tuple[list[Dataset], list[KeyCount]] = (
    ["musmus_11m"],
    ["counts", "counts-off-axis"],
)
param_names = ["dataset", "layer"]


### Dask-Based Benchmarks ###
def time_filter_cells_dask(*_):
    client = setup_dask_cluster()
    try:
        optimal_chunks = (
            adata.X.shape[0] // (4 * len(client.nthreads())),
            adata.X.shape[1],
        )
        adata.X = dd.from_array(adata.X, chunks=optimal_chunks).persist()
        sc.pp.filter_cells(adata, min_genes=100)
        assert adata.n_obs > 0  # Ensure cells are retained
    finally:
        client.close()


# def time_filter_cells_dask(*_):
#    client = setup_dask_cluster()
#    try:
#        # Compute optimal chunks based on Dask cluster
#        optimal_chunks = (adata.X.shape[0] // (4 * len(client.nthreads())), adata.X.shape[1])
#        adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
#        adata.X = adata.X.persist()  # Persist to avoid recomputation
#        sc.pp.filter_cells(adata, min_genes=100)
#    finally:
#        client.close()


def peakmem_filter_cells_dask(*_):
    client = setup_dask_cluster()
    try:
        optimal_chunks = (
            adata.X.shape[0] // (4 * len(client.nthreads())),
            adata.X.shape[1],
        )
        adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
        sc.pp.filter_cells(adata, min_genes=100)
    finally:
        client.close()


def time_filter_genes_dask(*_):
    client = setup_dask_cluster()
    try:
        optimal_chunks = (
            adata.X.shape[0] // (4 * len(client.nthreads())),
            adata.X.shape[1],
        )
        adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
        adata.X = adata.X.persist()
        sc.pp.filter_genes(adata, min_cells=3)
    finally:
        client.close()


def peakmem_filter_genes_dask(*_):
    client = setup_dask_cluster()
    try:
        optimal_chunks = (
            adata.X.shape[0] // (4 * len(client.nthreads())),
            adata.X.shape[1],
        )
        adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
        sc.pp.filter_genes(adata, min_cells=3)
    finally:
        client.close()


### General Dask and Non-Dask Preprocessing Benchmarks ###


class FastSuite:
    """Suite for benchmarking preprocessing operations with Dask."""

    params: tuple[list[Dataset], list[KeyCount]] = (
        ["musmus_11m"],
        ["counts", "counts-off-axis"],
    )
    param_names = ["dataset", "layer"]

    def time_calculate_qc_metrics_dask(self, *_):
        client = setup_dask_cluster()
        try:
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            adata.X = adata.X.persist()
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )
        finally:
            client.close()

    def peakmem_calculate_qc_metrics_dask(self, *_):
        client = setup_dask_cluster()
        try:
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )
        finally:
            client.close()

    def time_normalize_total_dask(self, *_):
        client = setup_dask_cluster()
        try:
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            adata.X = adata.X.map_blocks(
                lambda x: x / x.sum(axis=1), dtype=float
            )  # Optimize normalization
            adata.X = adata.X.persist()
        finally:
            client.close()

    def peakmem_normalize_total_dask(self, *_):
        client = setup_dask_cluster()
        try:
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            sc.pp.normalize_total(adata, target_sum=1e4)
        finally:
            client.close()

    def time_log1p_dask(self, *_):
        client = setup_dask_cluster()
        try:
            adata.uns.pop("log1p", None)
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            adata.X = adata.X.persist()
            sc.pp.log1p(adata)
        finally:
            client.close()

    def peakmem_log1p_dask(self, *_):
        client = setup_dask_cluster()
        try:
            adata.uns.pop("log1p", None)
            optimal_chunks = (
                adata.X.shape[0] // (4 * len(client.nthreads())),
                adata.X.shape[1],
            )
            adata.X = dd.from_array(adata.X, chunks=optimal_chunks)
            sc.pp.log1p(adata)
        finally:
            client.close()

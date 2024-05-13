"""
This module will benchmark preprocessing operations in Scanpy that run on counts
API documentation: https://scanpy.readthedocs.io/en/stable/api/preprocessing.html
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from ._utils import get_count_dataset

if TYPE_CHECKING:
    from anndata import AnnData

    from ._utils import Dataset

# setup variables

adata: AnnData
batch_key: str | None


def setup(dataset: Dataset, *_):
    """Setup global variables before each benchmark."""
    global adata, batch_key
    adata, batch_key = get_count_dataset(dataset)
    assert "log1p" not in adata.uns


# ASV suite

params: list[Dataset] = ["pbmc68k_reduced", "pbmc3k"]
param_names = ["dataset"]


def time_calculate_qc_metrics(*_):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def peakmem_calculate_qc_metrics(*_):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def time_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def peakmem_filter_cells(*_):
    sc.pp.filter_cells(adata, min_genes=100)


def time_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


def peakmem_filter_genes(*_):
    sc.pp.filter_genes(adata, min_cells=3)


def time_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


def peakmem_scrublet(*_):
    sc.pp.scrublet(adata, batch_key=batch_key)


def time_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def peakmem_normalize_total(*_):
    sc.pp.normalize_total(adata, target_sum=1e4)


def time_log1p(*_):
    # TODO: This would fail: assert "log1p" not in adata.uns, "ASV bug?"
    # https://github.com/scverse/scanpy/issues/3052
    sc.pp.log1p(adata)


def peakmem_log1p(*_):
    sc.pp.log1p(adata)

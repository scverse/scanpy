from __future__ import annotations

import warnings
from functools import cache
from typing import TYPE_CHECKING

import numpy as np
import pooch
from anndata import concat

import scanpy as sc

if TYPE_CHECKING:
    from anndata import AnnData


@cache
def _pbmc68k_reduced() -> AnnData:
    return sc.datasets.pbmc68k_reduced()


def pbmc68k_reduced() -> AnnData:
    return _pbmc68k_reduced().copy()


@cache
def _pbmc3k() -> AnnData:
    return sc.datasets.pbmc3k()


def pbmc3k() -> AnnData:
    return _pbmc3k().copy()


@cache
def _bmmc4k() -> AnnData:
    registry = pooch.create(
        path=pooch.os_cache("pooch"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    registry.load_registry_from_doi()
    samples = {smp: f"{smp}_filtered_feature_bc_matrix.h5" for smp in ("s1d1", "s1d3")}
    adatas = {}

    for sample_id, filename in samples.items():
        path = registry.fetch(filename)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"Variable names are not unique")
            sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        sc.pp.subsample(sample_adata, n_obs=4000 // len(samples))
        adatas[sample_id] = sample_adata

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Observation names are not unique")
        adata = concat(adatas, label="sample")
    adata.obs_names_make_unique()

    adata.var["is_mito"] = adata.var_names.str.startswith("MT-").sum()
    adata.obs["percent_mito"] = (
        np.sum(adata[:, adata.var["is_mito"]].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    )
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1

    return adata


def bmmc4k() -> AnnData:
    return _bmmc4k().copy()


@cache
def _lung93k() -> AnnData:
    path = pooch.retrieve(
        url="https://figshare.com/ndownloader/files/45788454",
        known_hash="md5:4f28af5ff226052443e7e0b39f3f9212",
    )
    return sc.read_h5ad(path)


def lung93k() -> AnnData:
    return _lung93k().copy()

from __future__ import annotations

import warnings
from functools import cache
from typing import TYPE_CHECKING

import numpy as np
import pooch
from anndata import concat

import scanpy as sc

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

    Dataset = Literal["pbmc68k_reduced", "pbmc3k", "bmmc", "lung93k"]


@cache
def _pbmc68k_reduced() -> AnnData:
    adata = sc.datasets.pbmc68k_reduced()
    # raw has the same number of genes, so we can use it for counts
    # it doesn’t actually contain counts for some reason, but close enough
    adata.layers["counts"] = adata.raw.X.copy()
    mapper = dict(
        percent_mito="pct_counts_mt",
        n_counts="total_counts",
    )
    adata.obs.rename(columns=mapper, inplace=True)
    return adata


def pbmc68k_reduced() -> AnnData:
    return _pbmc68k_reduced().copy()


@cache
def _pbmc3k() -> AnnData:
    adata = sc.datasets.pbmc3k()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata.layers["counts"] = adata.X.astype(np.int32, copy=True)
    sc.pp.log1p(adata)
    return adata


def pbmc3k() -> AnnData:
    return _pbmc3k().copy()


@cache
def _bmmc(n_obs: int = 4000) -> AnnData:
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
        sc.pp.subsample(sample_adata, n_obs=n_obs // len(samples))
        adatas[sample_id] = sample_adata

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Observation names are not unique")
        adata = concat(adatas, label="sample")
    adata.obs_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1

    return adata


def bmmc(n_obs: int = 400) -> AnnData:
    return _bmmc(n_obs).copy()


@cache
def _lung93k() -> AnnData:
    path = pooch.retrieve(
        url="https://figshare.com/ndownloader/files/45788454",
        known_hash="md5:4f28af5ff226052443e7e0b39f3f9212",
    )
    return sc.read_h5ad(path)


def lung93k() -> AnnData:
    return _lung93k().copy()


def get_dataset(dataset: Dataset) -> tuple[AnnData, str | None]:
    if dataset == "pbmc68k_reduced":
        return pbmc68k_reduced(), None
    if dataset == "pbmc3k":
        return pbmc3k(), None  # can’t use this with batches
    if dataset == "bmmc":
        # TODO: allow specifying bigger variant
        return bmmc(400), "sample"
    if dataset == "lung93k":
        return lung93k(), "PatientNumber"

    msg = f"Unknown dataset {dataset}"
    raise AssertionError(msg)


def get_count_dataset(dataset: Dataset) -> tuple[AnnData, str | None]:
    adata, batch_key = get_dataset(dataset)

    adata.X = adata.layers.pop("counts")
    # remove indicators that X was transformed
    adata.uns.pop("log1p", None)

    return adata, batch_key

from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING

import pooch

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
def _lung93k() -> AnnData:
    path = pooch.retrieve(
        url="https://figshare.com/ndownloader/files/45788454",
        known_hash="md5:4f28af5ff226052443e7e0b39f3f9212",
    )
    return sc.read_h5ad(path)


def lung93k() -> AnnData:
    return _lung93k().copy()

"""
Functions returning copies of datasets as cheaply as possible,
i.e. without having to hit the disk or (in case of ``_pbmc3k_normalized``) recomputing normalization.
"""

from __future__ import annotations

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    def cache(func):
        return lru_cache(maxsize=None)(func)


from anndata import AnnData
import scanpy as sc


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

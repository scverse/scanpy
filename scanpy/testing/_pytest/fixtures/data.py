"""
These fixtures provide a per test new copy of the dataset, without
having to hit the disk or (in case of ``_pbmc3k_normalized``) recomputing normalization.
The private fixtures create the object while the public ones return deep copies.
"""

import pytest
from anndata import AnnData

import scanpy as sc


@pytest.fixture(scope='session')
def _pbmc3k() -> AnnData:
    return sc.datasets.pbmc3k()


@pytest.fixture(scope="session")
def _pbmc3k_normalized(_pbmc3k) -> AnnData:
    pbmc = _pbmc3k.copy()
    pbmc.X = pbmc.X.astype("float64")  # For better accuracy
    sc.pp.filter_genes(pbmc, min_counts=1)
    sc.pp.log1p(pbmc)
    sc.pp.normalize_total(pbmc)
    sc.pp.highly_variable_genes(pbmc)
    return pbmc


@pytest.fixture(scope='session')
def _pbmc3k_processed() -> AnnData:
    return sc.datasets.pbmc3k_processed()


@pytest.fixture(scope='session')
def _pbmc68k_reduced() -> AnnData:
    return sc.datasets.pbmc68k_reduced()


@pytest.fixture(scope='session')
def _krumsiek11() -> AnnData:
    return sc.datasets.krumsiek11()


@pytest.fixture(scope='session')
def _paul15() -> AnnData:
    return sc.datasets.paul15()


@pytest.fixture
def pbmc3k(_pbmc3k) -> AnnData:
    return _pbmc3k.copy()


@pytest.fixture
def pbmc3k_normalized(_pbmc3k_normalized) -> AnnData:
    return _pbmc3k_normalized.copy()


@pytest.fixture
def pbmc3k_processed(_pbmc3k_processed) -> AnnData:
    return _pbmc3k_processed.copy()


@pytest.fixture
def pbmc68k_reduced(_pbmc68k_reduced) -> AnnData:
    return _pbmc68k_reduced.copy()


@pytest.fixture
def krumsiek11(_krumsiek11) -> AnnData:
    return _krumsiek11.copy()


@pytest.fixture
def paul15(_paul15) -> AnnData:
    return _paul15.copy()

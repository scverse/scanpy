"""
These fixtures provide a per test new copy of the dataset, without
having to hit the disk or (in case of ``_pbmc3k_normalized``) recomputing normalization.
The private fixtures create the object while the public ones return deep copies.
"""

from itertools import product
import pytest
import numpy as np
from scipy import sparse
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


@pytest.fixture(
    scope='session',
    params=list(
        product([sparse.csr_matrix.toarray, sparse.csr_matrix], ['float32', 'int64'])
    ),
    ids=lambda x: f'{x[0].__name__}-{x[1]}',
)
def _pbmc3ks_parametrized(request, _pbmc3k):
    sparsity_func, dtype = request.param
    return {
        small: _prepare_pbmc_testdata(_pbmc3k.copy(), sparsity_func, dtype, small=small)
        for small in [True, False]
    }


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
def pbmc3k_parametrized(_pbmc3ks_parametrized):
    return _pbmc3ks_parametrized[False].copy()


@pytest.fixture
def pbmc3k_parametrized_small(_pbmc3ks_parametrized):
    return _pbmc3ks_parametrized[True].copy()


@pytest.fixture
def pbmc68k_reduced(_pbmc68k_reduced) -> AnnData:
    return _pbmc68k_reduced.copy()


@pytest.fixture
def krumsiek11(_krumsiek11) -> AnnData:
    return _krumsiek11.copy()


@pytest.fixture
def paul15(_paul15) -> AnnData:
    return _paul15.copy()


def _prepare_pbmc_testdata(adata, sparsity_func, dtype, *, small: bool):
    """Prepares 3k PBMC dataset with batch key `batch` and defined datatype/sparsity.

    Params
    ------
    sparsity_func
        sparsity function applied to adata.X (e.g. csr_matrix.toarray for dense or csr_matrix for sparse)
    dtype
        numpy dtype applied to adata.X (e.g. 'float32' or 'int64')
    small
        False (default) returns full data, True returns small subset of the data.
    """

    import scanpy as sc

    if small:
        adata = adata[:1000, :500]
        sc.pp.filter_cells(adata, min_genes=1)
    np.random.seed(42)
    adata.obs['batch'] = np.random.randint(0, 3, size=adata.shape[0])
    sc.pp.filter_genes(adata, min_cells=1)
    adata.X = sparsity_func(adata.X.astype(dtype))
    return adata

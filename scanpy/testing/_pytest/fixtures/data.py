"""Fixtures for parametrized datasets."""

from __future__ import annotations

import tempfile
from itertools import product
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import DTypeLike


@pytest.fixture(
    scope="session",
    params=list(
        product([sparse.csr_matrix.toarray, sparse.csr_matrix], ["float32", "int64"])
    ),
    ids=lambda x: f"{x[0].__name__}-{x[1]}",
)
def _pbmc3ks_parametrized_session(request) -> dict[bool, AnnData]:
    from ..._helpers.data import pbmc3k

    sparsity_func, dtype = request.param
    return {
        small: _prepare_pbmc_testdata(pbmc3k(), sparsity_func, dtype, small=small)
        for small in [True, False]
    }


@pytest.fixture
def pbmc3k_parametrized(_pbmc3ks_parametrized_session) -> Callable[[], AnnData]:
    return _pbmc3ks_parametrized_session[False].copy


@pytest.fixture
def pbmc3k_parametrized_small(_pbmc3ks_parametrized_session) -> Callable[[], AnnData]:
    return _pbmc3ks_parametrized_session[True].copy


@pytest.fixture(scope="session")
def backed_adata_path():
    tmp_path = f"{tempfile.gettempdir()}/test.h5ad"
    X = sparse.random(200, 10, format="csr").astype(np.float32)
    cat = np.random.randint(0, 3, (200,))
    adata = AnnData(X, obs={"cat": cat})
    adata.obs["cat"] = adata.obs["cat"].astype("category")
    adata.write_h5ad(tmp_path)
    return tmp_path


def _prepare_pbmc_testdata(
    adata: AnnData,
    sparsity_func: Callable[
        [np.ndarray | sparse.spmatrix], np.ndarray | sparse.spmatrix
    ],
    dtype: DTypeLike,
    *,
    small: bool,
) -> AnnData:
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

    if small:
        adata = adata[:1000, :500].copy()
        sc.pp.filter_cells(adata, min_genes=1)
    np.random.seed(42)
    adata.obs["batch"] = np.random.randint(0, 3, size=adata.shape[0])
    sc.pp.filter_genes(adata, min_cells=1)
    adata.X = sparsity_func(adata.X.astype(dtype))
    return adata

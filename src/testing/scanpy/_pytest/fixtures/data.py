"""Fixtures for parametrized datasets."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING, cast

import numpy as np
import pytest
from anndata import AnnData, read_h5ad
from anndata import __version__ as anndata_version
from packaging.version import Version
from scipy import sparse

if Version(anndata_version) >= Version("0.10.0"):
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )

    if Version(anndata_version) >= Version("0.11.0rc2"):
        from anndata.io import sparse_dataset
    else:
        from anndata.experimental import sparse_dataset

    def make_sparse(x):
        return sparse_dataset(x)
else:
    from anndata._core.sparse_dataset import SparseDataset

    def make_sparse(x):
        return SparseDataset(x)


if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from numpy.typing import DTypeLike

    from scanpy._compat import CSBase, CSRBase


@pytest.fixture(
    scope="session",
    params=list(
        product([sparse.csr_matrix.toarray, sparse.csr_matrix], ["float32", "int64"])  # noqa: TID251
    ),
    ids=lambda x: f"{x[0].__name__}-{x[1]}",
)
def pbmc3ks_parametrized_session(request) -> dict[bool, AnnData]:
    from ..._helpers.data import pbmc3k

    sparsity_func, dtype = request.param
    return {
        small: _prepare_pbmc_testdata(pbmc3k(), sparsity_func, dtype, small=small)
        for small in [True, False]
    }


@pytest.fixture
def pbmc3k_parametrized(pbmc3ks_parametrized_session) -> Callable[[], AnnData]:
    return pbmc3ks_parametrized_session[False].copy


@pytest.fixture
def pbmc3k_parametrized_small(pbmc3ks_parametrized_session) -> Callable[[], AnnData]:
    return pbmc3ks_parametrized_session[True].copy


def random_csr(m: int, n: int) -> CSRBase:
    return sparse.random(m, n, format="csr")


@pytest.fixture(params=[np.random.randn, random_csr], ids=["sparse", "dense"])
def backed_adata(request: pytest.FixtureRequest, tmp_path: Path) -> AnnData:
    rand_func = cast("Callable[[int, int], np.ndarray | CSRBase]", request.param)
    X = rand_func(200, 10).astype(np.float32)
    cat = np.random.randint(0, 3, (X.shape[0],)).ravel()
    adata = AnnData(X, obs={"cat": cat})
    adata.obs["percent_mito"] = np.random.rand(X.shape[0])
    adata.obs["n_counts"] = X.sum(axis=1)
    adata.obs["cat"] = adata.obs["cat"].astype("category")
    adata.layers["X_copy"] = adata.X[...]
    adata.write_h5ad(tmp_path / "test.h5ad")
    adata = read_h5ad(tmp_path / "test.h5ad", backed="r")
    adata.layers["X_copy"] = (
        make_sparse(adata.file["X"])
        if isinstance(adata.X, SparseDataset)
        else adata.file["X"]
    )
    return adata


def _prepare_pbmc_testdata(
    adata: AnnData,
    sparsity_func: Callable[[np.ndarray | CSBase], np.ndarray | CSBase],
    dtype: DTypeLike,
    *,
    small: bool,
) -> AnnData:
    """Prepare 3k PBMC dataset with batch key `batch` and defined datatype/sparsity.

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
        adata = adata[:1000, :500].copy()
        sc.pp.filter_cells(adata, min_genes=1)
    np.random.seed(42)
    adata.obs["batch"] = np.random.randint(0, 3, size=adata.shape[0])
    sc.pp.filter_genes(adata, min_cells=1)
    adata.X = sparsity_func(adata.X.astype(dtype))
    return adata

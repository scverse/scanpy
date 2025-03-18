from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy.testing as npt
import pytest
from anndata import read_zarr

from scanpy._compat import DaskArray, ZappyArray
from scanpy.datasets._utils import filter_oldformatwarning
from scanpy.preprocessing import (
    filter_cells,
    filter_genes,
    log1p,
    normalize_per_cell,
    normalize_total,
)
from scanpy.preprocessing._distributed import materialize_as_ndarray
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from anndata import AnnData

HERE = Path(__file__).parent / Path("_data/")
input_file = Path(HERE, "10x-10k-subset.zarr")

DIST_TYPES = (DaskArray, ZappyArray)


pytestmark = [needs.zarr]


@pytest.fixture
@filter_oldformatwarning
def adata() -> AnnData:
    a = read_zarr(input_file)
    a.var_names_make_unique()
    a.X = a.X[:]  # convert to numpy array
    return a


@pytest.fixture(
    params=[
        pytest.param("direct", marks=[needs.zappy]),
        pytest.param("dask", marks=[needs.dask, pytest.mark.anndata_dask_support]),
    ]
)
@filter_oldformatwarning
def adata_dist(request: pytest.FixtureRequest) -> AnnData:
    # regular anndata except for X, which we replace on the next line
    a = read_zarr(input_file)
    a.var_names_make_unique()
    a.uns["dist-mode"] = request.param
    input_file_X = f"{input_file}/X"
    if request.param == "direct":
        import zappy.direct

        a.X = zappy.direct.from_zarr(input_file_X)
        return a

    assert request.param == "dask"
    import dask.array as da

    a.X = da.from_zarr(input_file_X)
    return a


def test_log1p(adata: AnnData, adata_dist: AnnData):
    log1p(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    log1p(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(result, adata.X)


@pytest.mark.filterwarnings("ignore:Use sc.pp.normalize_total instead:FutureWarning")
def test_normalize_per_cell(
    request: pytest.FixtureRequest, adata: AnnData, adata_dist: AnnData
):
    if isinstance(adata_dist.X, DaskArray):
        reason = "normalize_per_cell deprecated and broken for Dask"
        request.applymarker(pytest.mark.xfail(reason=reason))
    normalize_per_cell(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    normalize_per_cell(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(result, adata.X)


def test_normalize_total(adata: AnnData, adata_dist: AnnData):
    normalize_total(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    normalize_total(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(result, adata.X)


def test_filter_cells_array(adata: AnnData, adata_dist: AnnData):
    cell_subset_dist, number_per_cell_dist = filter_cells(adata_dist.X, min_genes=3)
    assert isinstance(cell_subset_dist, DIST_TYPES)
    assert isinstance(number_per_cell_dist, DIST_TYPES)

    cell_subset, number_per_cell = filter_cells(adata.X, min_genes=3)
    npt.assert_allclose(materialize_as_ndarray(cell_subset_dist), cell_subset)
    npt.assert_allclose(materialize_as_ndarray(number_per_cell_dist), number_per_cell)


def test_filter_cells(adata: AnnData, adata_dist: AnnData):
    filter_cells(adata_dist, min_genes=3)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    filter_cells(adata, min_genes=3)

    assert result.shape == adata.shape
    npt.assert_array_equal(adata_dist.obs["n_genes"], adata.obs["n_genes"])
    npt.assert_allclose(result, adata.X)


def test_filter_genes_array(adata: AnnData, adata_dist: AnnData):
    gene_subset_dist, number_per_gene_dist = filter_genes(adata_dist.X, min_cells=2)
    assert isinstance(gene_subset_dist, DIST_TYPES)
    assert isinstance(number_per_gene_dist, DIST_TYPES)

    gene_subset, number_per_gene = filter_genes(adata.X, min_cells=2)
    npt.assert_allclose(materialize_as_ndarray(gene_subset_dist), gene_subset)
    npt.assert_allclose(materialize_as_ndarray(number_per_gene_dist), number_per_gene)


def test_filter_genes(adata: AnnData, adata_dist: AnnData):
    filter_genes(adata_dist, min_cells=2)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    filter_genes(adata, min_cells=2)
    assert result.shape == adata.shape
    npt.assert_allclose(result, adata.X)


@filter_oldformatwarning
def test_write_zarr(adata: AnnData, adata_dist: AnnData):
    import zarr

    log1p(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    temp_store = zarr.TempStore()
    chunks = adata_dist.X.chunks
    if isinstance(chunks[0], tuple):
        chunks = (chunks[0][0],) + chunks[1]

    # write metadata using regular anndata
    adata.write_zarr(temp_store, chunks=chunks)
    if adata_dist.uns["dist-mode"] == "dask":
        adata_dist.X.to_zarr(temp_store.dir_path("X"), overwrite=True)
    elif adata_dist.uns["dist-mode"] == "direct":
        adata_dist.X.to_zarr(temp_store.dir_path("X"), chunks=chunks)
    else:
        pytest.fail("add branch for new dist-mode")

    # read back as zarr directly and check it is the same as adata.X
    adata_log1p = read_zarr(temp_store)

    log1p(adata)
    npt.assert_allclose(adata_log1p.X, adata.X)

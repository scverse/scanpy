from __future__ import annotations

from pathlib import Path

import numpy.testing as npt
import pytest
import zarr
from anndata import AnnData, read_zarr
from anndata.experimental import read_elem, sparse_dataset
from scipy import sparse as sp

from scanpy._compat import DaskArray, ZappyArray
from scanpy.datasets._utils import filter_oldformatwarning
from scanpy.preprocessing import (
    filter_cells,
    filter_genes,
    log1p,
    normalize_per_cell,
    normalize_total,
    scale,
)
from scanpy.preprocessing._distributed import materialize_as_ndarray
from scanpy.testing._helpers.data import sparse_dataset_as_dask
from scanpy.testing._pytest.marks import needs

HERE = Path(__file__).parent / Path("_data/")
input_file = Path(HERE, "10x-10k-subset.zarr")

DIST_TYPES = (DaskArray, ZappyArray)


pytestmark = [needs.zarr]


def read_w_sparse_dask(group: zarr.Group, obs_chunk: int = 100) -> AnnData:
    return AnnData(
        X=sparse_dataset_as_dask(sparse_dataset(group["layers/CSR_X"]), obs_chunk),
        **{
            k: read_elem(group[k]) if k in group else {}
            for k in ["obs", "var", "obsm", "varm", "uns", "obsp", "varp"]
        },
    )


@pytest.fixture()
@filter_oldformatwarning
def adata() -> AnnData:
    a = read_zarr(input_file)
    a.var_names_make_unique()
    a.X = a.X[:]  # convert to numpy array
    return a


@filter_oldformatwarning
@pytest.fixture(
    params=[
        pytest.param("direct", marks=[needs.zappy]),
        pytest.param("dask", marks=[needs.dask, pytest.mark.anndata_dask_support]),
        pytest.param(
            "sparse_chunks_in_dask",
            marks=[needs.dask, pytest.mark.anndata_dask_support],
        ),
    ]
)
def adata_dist(request: pytest.FixtureRequest) -> AnnData:
    # regular anndata except for X, which we replace on the next line
    if request.param == "sparse_chunks_in_dask":
        return read_w_sparse_dask(zarr.open(input_file))
    a = read_zarr(input_file)
    a.var_names_make_unique()
    input_file_X = f"{input_file}/X"
    if request.param == "direct":
        import zappy.direct

        a.X = zappy.direct.from_zarr(input_file_X)
        return a

    assert request.param == "dask"
    import dask.array as da

    a.X = da.from_zarr(input_file_X, chunks=(100, 1000))
    return a


def test_log1p(adata: AnnData, adata_dist: AnnData):
    log1p(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    log1p(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(
        result.toarray() if isinstance(result, sp.spmatrix) else result, adata.X
    )


def test_normalize_per_cell(
    request: pytest.FixtureRequest, adata: AnnData, adata_dist: AnnData
):
    if isinstance(adata_dist.X, DaskArray):
        msg = "normalize_per_cell deprecated and broken for Dask"
        request.node.add_marker(pytest.mark.xfail(reason=msg))
    normalize_per_cell(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    normalize_per_cell(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(
        result.toarray() if isinstance(result, sp.spmatrix) else result, adata.X
    )


def test_normalize_total(adata: AnnData, adata_dist: AnnData):
    normalize_total(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    normalize_total(adata)
    assert result.shape == adata.shape
    npt.assert_allclose(
        result.toarray() if isinstance(result, sp.spmatrix) else result,
        adata.X,
        atol=1e-6,
    )


def test_filter_cells_array(adata: AnnData, adata_dist: AnnData):
    cell_subset_dist, number_per_cell_dist = filter_cells(adata_dist.X, min_genes=3)
    assert isinstance(cell_subset_dist, DIST_TYPES)
    assert isinstance(number_per_cell_dist, DIST_TYPES)
    cell_subset_dist, number_per_cell_dist = materialize_as_ndarray(
        (cell_subset_dist, number_per_cell_dist)
    )

    cell_subset, number_per_cell = filter_cells(adata.X, min_genes=3)
    npt.assert_allclose(
        cell_subset_dist.toarray()
        if isinstance(cell_subset_dist, sp.spmatrix)
        else cell_subset_dist,
        cell_subset,
    )
    npt.assert_allclose(
        number_per_cell_dist.toarray()
        if isinstance(number_per_cell_dist, sp.spmatrix)
        else number_per_cell_dist,
        number_per_cell,
    )


def test_filter_cells(adata: AnnData, adata_dist: AnnData):
    filter_cells(adata_dist, min_genes=3)
    assert isinstance(adata_dist.X, DIST_TYPES)
    result = materialize_as_ndarray(adata_dist.X)
    filter_cells(adata, min_genes=3)

    assert result.shape == adata.shape
    npt.assert_array_equal(adata_dist.obs["n_genes"], adata.obs["n_genes"])
    npt.assert_allclose(
        result.toarray() if isinstance(result, sp.spmatrix) else result, adata.X
    )


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
    npt.assert_allclose(
        result.toarray() if isinstance(result, sp.spmatrix) else result, adata.X
    )


@filter_oldformatwarning
def test_write_zarr(adata: AnnData, adata_dist: AnnData):
    import zarr

    log1p(adata_dist)
    assert isinstance(adata_dist.X, DIST_TYPES)
    temp_store = zarr.TempStore()  # write_elem needs a path
    chunks = adata_dist.X.chunks
    if isinstance(chunks[0], tuple):
        chunks = (chunks[0][0],) + chunks[1]

    adata_dist.write_zarr(temp_store)
    # read back as zarr directly and check it is the same as adata.X
    adata_log1p = read_zarr(temp_store)

    log1p(adata)
    expected = adata.X
    actual = adata_log1p.X
    npt.assert_allclose(
        expected.toarray() if isinstance(expected, sp.spmatrix) else expected,
        actual.toarray() if isinstance(actual, sp.spmatrix) else actual,
    )


def test_scale(adata: AnnData, adata_dist: AnnData):
    adata_zero_centered = scale(adata, copy=True)
    if sp.issparse(adata_dist.X._meta):
        with pytest.warns(
            UserWarning, match="zero-center being used with `DaskArray` sparse chunks"
        ):
            adata_dist_zero_centered = scale(adata_dist, copy=True)
    else:
        adata_dist_zero_centered = scale(adata_dist, copy=True)
    npt.assert_allclose(
        adata_zero_centered.X, adata_dist_zero_centered.X, rtol=1e-6, atol=1e-6
    )

    adata_non_zero_centered = scale(adata, copy=True, zero_center=False)
    adata_dist_non_zero_centered = scale(adata_dist, copy=True, zero_center=False)
    computed_X = adata_dist_non_zero_centered.X.compute()
    if sp.issparse(computed_X):
        computed_X = (
            computed_X.todense()
        )  # this is a COO matrix now because of the non-zero-centering
    npt.assert_allclose(adata_non_zero_centered.X, computed_X, rtol=1e-6, atol=1e-6)

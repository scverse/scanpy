from importlib.util import find_spec
from pathlib import Path

import anndata as ad
import numpy.testing as npt
import pytest

from scanpy.preprocessing import normalize_total, filter_genes
from scanpy.preprocessing import log1p, normalize_per_cell, filter_cells
from scanpy.preprocessing._distributed import materialize_as_ndarray
from scanpy.testing._pytest.marks import needs


HERE = Path(__file__).parent / Path('_data/')
input_file = str(Path(HERE, "10x-10k-subset.zarr"))


@needs("zarr")
class TestPreprocessingDistributed:
    @pytest.fixture()
    def adata(self):
        a = ad.read_zarr(input_file)  # regular anndata
        a.X = a.X[:]  # convert to numpy array
        return a

    @pytest.fixture(
        params=[
            pytest.param("direct", marks=[needs("zappy")]),
            pytest.param("dask", marks=[needs("dask")]),
        ]
    )
    def adata_dist(self, request):
        # regular anndata except for X, which we replace on the next line
        a = ad.read_zarr(input_file)
        a.uns["dist-mode"] = request.param
        input_file_X = f"{input_file}/X"
        if request.param == "direct":
            import zappy.direct

            a.X = zappy.direct.from_zarr(input_file_X)
            yield a
        elif request.param == "dask":
            import dask.array as da

            a.X = da.from_zarr(input_file_X)
            yield a

    def test_log1p(self, adata, adata_dist):
        log1p(adata_dist)
        result = materialize_as_ndarray(adata_dist.X)
        log1p(adata)
        assert result.shape == adata.shape
        assert result.shape == (adata.n_obs, adata.n_vars)
        npt.assert_allclose(result, adata.X)

    def test_normalize_per_cell(self, adata, adata_dist):
        if adata_dist.uns["dist-mode"] == "dask":
            pytest.xfail("TODO: Test broken for dask")
        normalize_per_cell(adata_dist)
        result = materialize_as_ndarray(adata_dist.X)
        normalize_per_cell(adata)
        assert result.shape == adata.shape
        assert result.shape == (adata.n_obs, adata.n_vars)
        npt.assert_allclose(result, adata.X)

    def test_normalize_total(self, adata, adata_dist):
        normalize_total(adata_dist)
        result = materialize_as_ndarray(adata_dist.X)
        normalize_total(adata)
        assert result.shape == adata.shape
        assert result.shape == (adata.n_obs, adata.n_vars)
        npt.assert_allclose(result, adata.X)

    def test_filter_cells(self, adata, adata_dist):
        filter_cells(adata_dist, min_genes=3)
        result = materialize_as_ndarray(adata_dist.X)
        filter_cells(adata, min_genes=3)
        assert result.shape == adata.shape
        assert result.shape == (adata.n_obs, adata.n_vars)
        npt.assert_allclose(result, adata.X)

    def test_filter_genes(self, adata, adata_dist):
        filter_genes(adata_dist, min_cells=2)
        result = materialize_as_ndarray(adata_dist.X)
        filter_genes(adata, min_cells=2)
        assert result.shape == adata.shape
        assert result.shape == (adata.n_obs, adata.n_vars)
        npt.assert_allclose(result, adata.X)

    def test_write_zarr(self, adata, adata_dist):
        import zarr

        log1p(adata_dist)
        temp_store = zarr.TempStore()
        chunks = adata_dist.X.chunks
        if isinstance(chunks[0], tuple):
            chunks = (chunks[0][0],) + chunks[1]
        # write metadata using regular anndata
        adata.write_zarr(temp_store, chunks)
        if adata_dist["dist-mode"] == "dask":
            adata_dist.X.to_zarr(temp_store.dir_path("X"), overwrite=True)
        elif adata_dist["dist-mode"] == "direct":
            adata_dist.X.to_zarr(temp_store.dir_path("X"), chunks)
        else:
            assert False, "add branch for new dist-mode"
        # read back as zarr directly and check it is the same as adata.X
        adata_log1p = ad.read_zarr(temp_store)
        log1p(adata)
        npt.assert_allclose(adata_log1p.X, adata.X)

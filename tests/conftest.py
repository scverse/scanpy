import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
from anndata import AnnData
from scipy import sparse
from anndata.tests.helpers import asarray, assert_equal

mpl.use('agg')
from matplotlib import pyplot
from matplotlib.testing.compare import compare_images
import pytest

import scanpy

scanpy.settings.verbosity = "hint"

# define this after importing scanpy but before running tests
IMPORTED = frozenset(sys.modules.keys())


@pytest.fixture
def imported_modules():
    return IMPORTED


def make_comparer(path_expected: Path, path_actual: Path, *, tol: int):
    def save_and_compare(basename, tolerance=None):
        path_actual.mkdir(parents=True, exist_ok=True)
        out_path = path_actual / f'{basename}.png'
        pyplot.savefig(out_path, dpi=40)
        pyplot.close()
        if tolerance is None:
            tolerance = tol
        res = compare_images(
            str(path_expected / f'{basename}.png'), str(out_path), tolerance
        )
        assert res is None, res

    return save_and_compare


@pytest.fixture
def image_comparer():
    return make_comparer


@pytest.fixture
def plt():
    return pyplot


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, asarray],
    ids=["scipy-csr", "scipy-csc", "np-ndarray"],
)
def array_type(request):
    """Function which converts passed array to one of the common array types."""
    return request.param


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param


###########################
# Representation choice
###########################
# These functions can be used to check that functions are correctly using arugments
# like `layers`, `obsm`, etc.


@pytest.fixture(scope='session')
def check_rep_mutation():
    def check_rep_mutation(func, X, **kwargs):
        """Check that only the array meant to be modified is modified."""
        adata = AnnData(
            X=X.copy(),
            layers={"layer": X.copy()},
            obsm={"obsm": X.copy()},
            dtype=X.dtype,
        )
        adata_X = func(adata, copy=True, **kwargs)
        adata_layer = func(adata, layer="layer", copy=True, **kwargs)
        adata_obsm = func(adata, obsm="obsm", copy=True, **kwargs)

        assert np.array_equal(asarray(adata_X.X), asarray(adata_layer.layers["layer"]))
        assert np.array_equal(asarray(adata_X.X), asarray(adata_obsm.obsm["obsm"]))

        assert np.array_equal(asarray(adata_layer.X), asarray(adata_layer.obsm["obsm"]))
        assert np.array_equal(
            asarray(adata_obsm.X), asarray(adata_obsm.layers["layer"])
        )
        assert np.array_equal(
            asarray(adata_X.layers["layer"]), asarray(adata_X.obsm["obsm"])
        )

    return check_rep_mutation


@pytest.fixture(scope='session')
def check_rep_results():
    def check_rep_results(func, X, **kwargs):
        """Checks that the results of a computation add values/
        mutate the anndata object in a consistent way.
        """
        # Gen data
        adata_X = AnnData(
            X=X.copy(),
            layers={"layer": np.zeros(shape=X.shape, dtype=X.dtype)},
            obsm={"obsm": np.zeros(shape=X.shape, dtype=X.dtype)},
        )
        adata_layer = AnnData(
            X=np.zeros(shape=X.shape, dtype=X.dtype),
            layers={"layer": X.copy()},
            obsm={"obsm": np.zeros(shape=X.shape, dtype=X.dtype)},
        )
        adata_obsm = AnnData(
            X=np.zeros(shape=X.shape, dtype=X.dtype),
            layers={"layer": np.zeros(shape=X.shape, dtype=X.dtype)},
            obsm={"obsm": X.copy()},
        )

        # Apply function
        func(adata_X, **kwargs)
        func(adata_layer, layer="layer", **kwargs)
        func(adata_obsm, obsm="obsm", **kwargs)

        # Reset X
        adata_X.X = np.zeros(shape=X.shape, dtype=X.dtype)
        adata_layer.layers["layer"] = np.zeros(shape=X.shape, dtype=X.dtype)
        adata_obsm.obsm["obsm"] = np.zeros(shape=X.shape, dtype=X.dtype)

        # Check equality
        assert_equal(adata_X, adata_layer)
        assert_equal(adata_X, adata_obsm)

    return check_rep_results

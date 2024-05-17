from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import pytest
from anndata import AnnData
from anndata.tests.helpers import (
    asarray,
    assert_equal,
)
from packaging.version import Version
from scipy import sparse
from scipy.sparse import issparse

import scanpy as sc
from testing.scanpy._helpers import as_dense_dask_array, as_sparse_dask_array
from testing.scanpy._helpers.data import pbmc3k_normalized
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import (
    ARRAY_TYPES,
    ARRAY_TYPES_SPARSE_DASK_UNSUPPORTED,
    param_with,
)

if TYPE_CHECKING:
    from typing import Literal

A_list = np.array(
    [
        [0, 0, 7, 0, 0],
        [8, 5, 0, 2, 0],
        [6, 0, 0, 2, 5],
        [0, 0, 0, 1, 0],
        [8, 8, 2, 1, 0],
        [0, 0, 0, 4, 5],
    ]
)

A_pca = np.array(
    [
        [-4.4783009, 5.55508466, 1.73111572, -0.06029139, 0.17292555],
        [5.4855141, -0.42651191, -0.74776055, -0.74532146, 0.74633582],
        [0.01161428, -4.0156662, 2.37252748, -1.33122372, -0.29044446],
        [-3.61934397, 0.48525412, -2.96861931, -1.16312545, -0.33230607],
        [7.14050048, 1.86330409, -0.05786325, 1.25045782, -0.50213107],
        [-4.53998399, -3.46146476, -0.32940009, 2.04950419, 0.20562023],
    ]
)

A_svd = np.array(
    [
        [-0.77034038, -2.00750922, 6.64603489, -0.39669256, -0.22212097],
        [-9.47135856, -0.6326006, -1.33787112, -0.24894361, -1.02044665],
        [-5.90007339, 4.99658727, 0.70712592, -2.15188849, 0.30430008],
        [-0.19132409, 0.42172251, 0.11169531, 0.50977966, -0.71637566],
        [-11.1286238, -2.73045559, 0.08040596, 1.06850585, 0.74173764],
        [-1.50180389, 5.56886849, 1.64034442, 2.24476032, -0.05109001],
    ]
)


# If one uses dask for PCA it will always require dask-ml
@pytest.fixture(
    params=[
        param_with(at, marks=[needs.dask_ml]) if "dask" in at.id else at
        for at in ARRAY_TYPES_SPARSE_DASK_UNSUPPORTED
    ]
)
def array_type(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture(params=[None, "valid", "invalid"])
def svd_solver_type(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture(params=[True, False], ids=["zero_center", "no_zero_center"])
def zero_center(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture
def pca_params(
    array_type, svd_solver_type: Literal[None, "valid", "invalid"], zero_center
):
    all_svd_solvers = {"auto", "full", "arpack", "randomized", "tsqr", "lobpcg"}

    expected_warning = None
    svd_solver = None
    if svd_solver_type is not None:
        # TODO: are these right for sparse?
        if array_type in {as_dense_dask_array, as_sparse_dask_array}:
            svd_solver = (
                {"auto", "full", "tsqr", "randomized"}
                if zero_center
                else {"tsqr", "randomized"}
            )
        elif array_type in {sparse.csr_matrix, sparse.csc_matrix}:
            svd_solver = (
                {"lobpcg", "arpack"} if zero_center else {"arpack", "randomized"}
            )
        elif array_type is asarray:
            svd_solver = (
                {"auto", "full", "arpack", "randomized"}
                if zero_center
                else {"arpack", "randomized"}
            )
        else:
            assert False, f"Unknown array type {array_type}"
        if svd_solver_type == "invalid":
            svd_solver = all_svd_solvers - svd_solver
            expected_warning = "Ignoring"

        svd_solver = np.random.choice(list(svd_solver))
    # explicit check for special case
    if (
        svd_solver == "randomized"
        and zero_center
        and array_type in [sparse.csr_matrix, sparse.csc_matrix]
    ):
        expected_warning = "not work with sparse input"

    return (svd_solver, expected_warning)


def test_pca_warnings(array_type, zero_center, pca_params):
    svd_solver, expected_warning = pca_params
    A = array_type(A_list).astype("float32")
    adata = AnnData(A)

    if expected_warning is not None:
        with pytest.warns(UserWarning, match=expected_warning):
            sc.pp.pca(adata, svd_solver=svd_solver, zero_center=zero_center)
        return

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            warnings.filterwarnings(
                "ignore",
                "pkg_resources is deprecated as an API",
                DeprecationWarning,
            )
            sc.pp.pca(adata, svd_solver=svd_solver, zero_center=zero_center)
    except UserWarning:
        # TODO: Fix this case, maybe by increasing test data size.
        # https://github.com/scverse/scanpy/issues/2744
        if svd_solver == "lobpcg":
            pytest.xfail(reason="lobpcg doesn’t work with this small test data")
        raise


# This warning test is out of the fixture because it is a special case in the logic of the function
def test_pca_warnings_sparse():
    for array_type in (sparse.csr_matrix, sparse.csc_matrix):
        A = array_type(A_list).astype("float32")
        adata = AnnData(A)
        with pytest.warns(UserWarning, match="not work with sparse input"):
            sc.pp.pca(adata, svd_solver="randomized", zero_center=True)


def test_pca_transform(array_type):
    A = array_type(A_list).astype("float32")
    A_pca_abs = np.abs(A_pca)
    A_svd_abs = np.abs(A_svd)

    adata = AnnData(A)

    with warnings.catch_warnings(record=True) as record:
        sc.pp.pca(adata, n_comps=4, zero_center=True, dtype="float64")
    assert len(record) == 0, record

    assert np.linalg.norm(A_pca_abs[:, :4] - np.abs(adata.obsm["X_pca"])) < 2e-05

    with warnings.catch_warnings(record=True) as record:
        sc.pp.pca(
            adata,
            n_comps=5,
            zero_center=True,
            svd_solver="randomized",
            dtype="float64",
            random_state=14,
        )
    if sparse.issparse(A):
        assert any(
            isinstance(r.message, UserWarning)
            and "svd_solver 'randomized' does not work with sparse input"
            in str(r.message)
            for r in record
        )
    else:
        assert len(record) == 0

    assert np.linalg.norm(A_pca_abs - np.abs(adata.obsm["X_pca"])) < 2e-05

    with warnings.catch_warnings(record=True) as record:
        sc.pp.pca(adata, n_comps=4, zero_center=False, dtype="float64", random_state=14)
    assert len(record) == 0

    assert np.linalg.norm(A_svd_abs[:, :4] - np.abs(adata.obsm["X_pca"])) < 2e-05


def test_pca_shapes():
    """
    Tests that n_comps behaves correctly
    See https://github.com/scverse/scanpy/issues/1051
    """
    adata = AnnData(np.random.randn(30, 20))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (30, 19)

    adata = AnnData(np.random.randn(20, 30))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (20, 19)

    with pytest.raises(ValueError):
        sc.pp.pca(adata, n_comps=100)


def test_pca_sparse():
    """
    Tests that implicitly centered pca on sparse arrays returns equivalent results to
    explicit centering on dense arrays.
    """
    pbmc = pbmc3k_normalized()

    pbmc_dense = pbmc.copy()
    pbmc_dense.X = pbmc_dense.X.toarray()

    implicit = sc.pp.pca(pbmc, dtype=np.float64, copy=True)
    explicit = sc.pp.pca(pbmc_dense, dtype=np.float64, copy=True)

    np.testing.assert_allclose(
        implicit.uns["pca"]["variance"], explicit.uns["pca"]["variance"]
    )
    np.testing.assert_allclose(
        implicit.uns["pca"]["variance_ratio"], explicit.uns["pca"]["variance_ratio"]
    )
    np.testing.assert_allclose(implicit.obsm["X_pca"], explicit.obsm["X_pca"])
    np.testing.assert_allclose(implicit.varm["PCs"], explicit.varm["PCs"])


def test_pca_reproducible(array_type):
    pbmc = pbmc3k_normalized()
    pbmc.X = array_type(pbmc.X)

    a = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=42)
    b = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=42)
    c = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=0)

    assert_equal(a, b)
    # Test that changing random seed changes result
    # Does not show up reliably with 32 bit computation
    assert not np.array_equal(a.obsm["X_pca"], c.obsm["X_pca"])


def test_pca_chunked():
    """
    See https://github.com/scverse/scanpy/issues/1590
    But this is also a more general test
    """

    # Subsetting for speed of test
    pbmc_full = pbmc3k_normalized()
    pbmc = pbmc_full[::6].copy()
    pbmc.X = pbmc.X.astype(np.float64)
    chunked = sc.pp.pca(pbmc_full, chunked=True, copy=True)
    default = sc.pp.pca(pbmc_full, copy=True)

    # Taking absolute value since sometimes dimensions are flipped
    np.testing.assert_allclose(
        np.abs(chunked.obsm["X_pca"]), np.abs(default.obsm["X_pca"])
    )
    np.testing.assert_allclose(np.abs(chunked.varm["PCs"]), np.abs(default.varm["PCs"]))
    np.testing.assert_allclose(
        np.abs(chunked.uns["pca"]["variance"]), np.abs(default.uns["pca"]["variance"])
    )
    np.testing.assert_allclose(
        np.abs(chunked.uns["pca"]["variance_ratio"]),
        np.abs(default.uns["pca"]["variance_ratio"]),
    )


def test_pca_n_pcs():
    """
    Tests that the n_pcs parameter also works for
    representations not called "X_pca"
    """
    pbmc = pbmc3k_normalized()
    sc.pp.pca(pbmc, dtype=np.float64)
    pbmc.obsm["X_pca_test"] = pbmc.obsm["X_pca"]
    original = sc.pp.neighbors(pbmc, n_pcs=5, use_rep="X_pca", copy=True)
    renamed = sc.pp.neighbors(pbmc, n_pcs=5, use_rep="X_pca_test", copy=True)

    assert np.allclose(original.obsm["X_pca"], renamed.obsm["X_pca_test"])
    assert np.allclose(
        original.obsp["distances"].toarray(), renamed.obsp["distances"].toarray()
    )


# We use all ARRAY_TYPES here since this error should be raised before
# PCA can realize that it got a Dask array
@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_mask_highly_var_error(array_type):
    """Check if use_highly_variable=True throws an error if the annotation is missing."""
    adata = AnnData(array_type(A_list).astype("float32"))
    with (
        pytest.warns(
            FutureWarning,
            match=r"Argument `use_highly_variable` is deprecated, consider using the mask argument\.",
        ),
        pytest.raises(
            ValueError,
            match=r"Did not find `adata\.var\['highly_variable'\]`\.",
        ),
    ):
        sc.pp.pca(adata, use_highly_variable=True)


def test_mask_length_error():
    """Check error for n_obs / mask length mismatch."""
    adata = AnnData(A_list)
    mask_var = np.random.choice([True, False], adata.shape[1] + 1)
    with pytest.raises(
        ValueError, match=r"The shape of the mask do not match the data\."
    ):
        sc.pp.pca(adata, mask_var=mask_var, copy=True)


def test_mask_var_argument_equivalence(float_dtype, array_type):
    """Test if pca result is equal when given mask as boolarray vs string"""

    adata_base = AnnData(array_type(np.random.random((100, 10))).astype(float_dtype))
    mask_var = np.random.choice([True, False], adata_base.shape[1])

    adata = adata_base.copy()
    sc.pp.pca(adata, mask_var=mask_var, dtype=float_dtype)

    adata_w_mask = adata_base.copy()
    adata_w_mask.var["mask"] = mask_var
    sc.pp.pca(adata_w_mask, mask_var="mask", dtype=float_dtype)

    assert np.allclose(
        adata.X.toarray() if issparse(adata.X) else adata.X,
        adata_w_mask.X.toarray() if issparse(adata_w_mask.X) else adata_w_mask.X,
    )


def test_mask(array_type, request):
    if array_type is as_dense_dask_array:
        pytest.xfail("TODO: Dask arrays are not supported")
    adata = sc.datasets.blobs(n_variables=10, n_centers=3, n_observations=100)
    adata.X = array_type(adata.X)

    if isinstance(adata.X, np.ndarray) and Version(ad.__version__) < Version("0.9"):
        request.node.add_marker(
            pytest.mark.xfail(
                reason="TODO: Previous version of anndata would return an F ordered array for one"
                " case here, which suprisingly considerably changes the results of PCA. "
            )
        )
    mask_var = np.random.choice([True, False], adata.shape[1])

    adata_masked = adata[:, mask_var].copy()
    sc.pp.pca(adata, mask_var=mask_var)
    sc.pp.pca(adata_masked)

    masked_var_loadings = adata.varm["PCs"][~mask_var]
    np.testing.assert_equal(masked_var_loadings, np.zeros_like(masked_var_loadings))

    np.testing.assert_equal(adata.obsm["X_pca"], adata_masked.obsm["X_pca"])
    # There are slight difference based on whether the matrix was column or row major
    np.testing.assert_allclose(
        adata.varm["PCs"][mask_var], adata_masked.varm["PCs"], rtol=1e-11
    )


def test_mask_order_warning(request):
    if Version(ad.__version__) >= Version("0.9"):
        request.node.add_marker(
            pytest.mark.xfail(
                reason="Not expected to warn in later versions of anndata"
            )
        )

    adata = ad.AnnData(X=np.random.randn(50, 5))
    mask = np.array([True, False, True, False, True])

    with pytest.warns(
        UserWarning,
        match="When using a mask parameter with anndata<0.9 on a dense array",
    ):
        sc.pp.pca(adata, mask_var=mask)


def test_mask_defaults(array_type, float_dtype):
    """
    Test if pca result is equal without highly variable and with-but mask is None
    and if pca takes highly variable as mask as default
    """
    A = array_type(A_list).astype("float64")
    adata = AnnData(A)

    without_var = sc.pp.pca(adata, copy=True, dtype=float_dtype)

    rng = np.random.default_rng(8)
    mask = rng.choice([True, False], adata.shape[1])
    adata.var["highly_variable"] = mask
    with_var = sc.pp.pca(adata, copy=True, dtype=float_dtype)
    assert without_var.uns["pca"]["params"]["mask_var"] is None
    assert with_var.uns["pca"]["params"]["mask_var"] == "highly_variable"
    assert not np.array_equal(without_var.obsm["X_pca"], with_var.obsm["X_pca"])
    with_no_mask = sc.pp.pca(adata, mask_var=None, copy=True, dtype=float_dtype)
    assert np.array_equal(without_var.obsm["X_pca"], with_no_mask.obsm["X_pca"])


def test_pca_layer():
    """
    Tests that layers works the same way as .X
    """
    X_adata = pbmc3k_normalized()

    layer_adata = X_adata.copy()
    layer_adata.layers["counts"] = X_adata.X.copy()
    del layer_adata.X

    sc.pp.pca(X_adata)
    sc.pp.pca(layer_adata, layer="counts")

    assert layer_adata.uns["pca"]["params"]["layer"] == "counts"
    assert "layer" not in X_adata.uns["pca"]["params"]

    np.testing.assert_equal(
        X_adata.uns["pca"]["variance"], layer_adata.uns["pca"]["variance"]
    )
    np.testing.assert_equal(
        X_adata.uns["pca"]["variance_ratio"], layer_adata.uns["pca"]["variance_ratio"]
    )
    np.testing.assert_equal(X_adata.obsm["X_pca"], layer_adata.obsm["X_pca"])
    np.testing.assert_equal(X_adata.varm["PCs"], layer_adata.varm["PCs"])

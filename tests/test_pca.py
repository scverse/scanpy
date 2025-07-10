from __future__ import annotations

import warnings
from contextlib import nullcontext
from functools import wraps
from typing import TYPE_CHECKING, Literal

import anndata as ad
import numpy as np
import pytest
from anndata import AnnData
from anndata.tests import helpers
from anndata.tests.helpers import assert_equal
from packaging.version import Version
from scipy import sparse

import scanpy as sc
from scanpy._compat import CSBase, DaskArray, pkg_version
from scanpy._utils import get_literal_vals
from scanpy.preprocessing._pca import SvdSolver as SvdSolverSupported
from scanpy.preprocessing._pca._dask import _cov_sparse_dask
from testing.scanpy import _helpers
from testing.scanpy._helpers.data import pbmc3k_normalized
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES as ARRAY_TYPES_ALL
from testing.scanpy._pytest.params import param_with

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    from anndata.typing import ArrayDataStructureType

    ArrayType = Callable[[np.ndarray], ArrayDataStructureType]


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


if pkg_version("anndata") < Version("0.9"):

    def to_memory(self: AnnData, *, copy: bool = False) -> AnnData:
        """Compatibility version of AnnData.to_memory() that works with old AnnData versions."""
        adata = self
        if adata.isbacked:
            adata = adata.to_memory()
        return adata.copy() if copy else adata
else:
    to_memory = AnnData.to_memory


def _chunked_1d(
    f: Callable[[np.ndarray], DaskArray],
) -> Callable[[np.ndarray], DaskArray]:
    @wraps(f)
    def wrapper(a: np.ndarray) -> DaskArray:
        da = f(a)
        return da.rechunk((da.chunksize[0], -1))

    return wrapper


DASK_CONVERTERS = {
    f: _chunked_1d(f)
    for f in (_helpers.as_dense_dask_array, _helpers.as_sparse_dask_array)
}


def maybe_convert_array_to_dask(array_type):
    # If one uses dask for PCA it will always require dask-ml.
    # dask-ml can’t do 2D-chunked arrays, so rechunk them.
    if as_dask_array := DASK_CONVERTERS.get(array_type):
        return (as_dask_array,)

    # When not using dask, just return the array type
    assert "dask" not in array_type.__name__, "add more branches or refactor"
    return (array_type,)


ARRAY_TYPES = [
    param_with(
        at,
        maybe_convert_array_to_dask,
        marks=[needs.dask_ml] if at.id == "dask_array_dense" else [],
    )
    for at in ARRAY_TYPES_ALL
]


@pytest.fixture(params=ARRAY_TYPES)
def array_type(request: pytest.FixtureRequest) -> ArrayType:
    return request.param


SVDSolverDeprecated = Literal["lobpcg"]
SVDSolver = SvdSolverSupported | SVDSolverDeprecated

SKLEARN_ADDITIONAL: frozenset[SvdSolverSupported] = frozenset(
    {"covariance_eigh"} if pkg_version("scikit-learn") >= Version("1.5") else ()
)


def gen_pca_params(
    *,
    array_type: ArrayType,
    svd_solver_type: Literal[None, "valid", "invalid"],
    zero_center: bool,
) -> Generator[tuple[SVDSolver | None, str | None, str | None], None, None]:
    if array_type is DASK_CONVERTERS[_helpers.as_sparse_dask_array] and not zero_center:
        xfail_reason = "Sparse-in-dask with zero_center=False not implemented yet"
        yield None, None, xfail_reason
        return
    if svd_solver_type is None:
        yield None, None, None
        return

    svd_solvers, warn_pat_expected = possible_solvers(
        array_type=array_type, svd_solver_type=svd_solver_type, zero_center=zero_center
    )

    # sorted to prevent https://github.com/pytest-dev/pytest-xdist/issues/432
    for svd_solver in sorted(svd_solvers):
        # explicit check for special case
        if (
            isinstance(array_type, type)
            and issubclass(array_type, CSBase)
            and zero_center
            and svd_solver == "lobpcg"
        ):
            pat = r"legacy code"
        else:
            pat = warn_pat_expected
        yield (svd_solver, pat, None)


def possible_solvers(
    *,
    array_type: ArrayType,
    svd_solver_type: Literal["valid", "invalid"],
    zero_center: bool,
) -> tuple[set[SVDSolver], str | None]:
    all_svd_solvers = get_literal_vals(SVDSolver)
    svd_solvers: set[SVDSolver]
    match array_type, zero_center:
        case (dc, True) if dc is DASK_CONVERTERS[_helpers.as_dense_dask_array]:
            svd_solvers = {"auto", "full", "tsqr", "randomized", "covariance_eigh"}
        case (dc, False) if dc is DASK_CONVERTERS[_helpers.as_dense_dask_array]:
            svd_solvers = {"tsqr", "randomized"}
        case (dc, True) if dc is DASK_CONVERTERS[_helpers.as_sparse_dask_array]:
            svd_solvers = {"covariance_eigh"}
        case (type() as dc, True) if issubclass(dc, CSBase):
            svd_solvers = {"arpack"} | SKLEARN_ADDITIONAL
        case (type() as dc, False) if issubclass(dc, CSBase):
            svd_solvers = {"arpack", "randomized"}
        case (helpers.asarray, True):
            svd_solvers = {"auto", "full", "arpack", "randomized"} | SKLEARN_ADDITIONAL
        case (helpers.asarray, False):
            svd_solvers = {"arpack", "randomized"}
        case _:
            pytest.fail(f"Unknown {array_type=} ({zero_center=})")

    if svd_solver_type == "invalid":
        svd_solvers = all_svd_solvers - svd_solvers
        warn_pat_expected = r"Ignoring svd_solver"
    elif svd_solver_type == "valid":
        warn_pat_expected = None
    else:
        pytest.fail(f"Unknown {svd_solver_type=}")
    return svd_solvers, warn_pat_expected


@pytest.mark.parametrize(
    ("array_type", "zero_center", "svd_solver", "warn_pat_expected"),
    [
        pytest.param(
            array_type.values[0],
            zero_center,
            svd_solver,
            warn_pat_expected,
            marks=(
                array_type.marks
                if xfail_reason is None
                else [pytest.mark.xfail(reason=xfail_reason)]
            ),
            id=(
                f"{array_type.id}-{'zero_center' if zero_center else 'no_zero_center'}-"
                f"{svd_solver or svd_solver_type}-{'xfail' if xfail_reason else warn_pat_expected}"
            ),
        )
        for array_type in ARRAY_TYPES
        for zero_center in [True, False]
        for svd_solver_type in [None, "valid", "invalid"]
        for svd_solver, warn_pat_expected, xfail_reason in gen_pca_params(
            array_type=array_type.values[0],
            zero_center=zero_center,
            svd_solver_type=svd_solver_type,
        )
    ],
)
def test_pca_warnings(
    *,
    array_type: ArrayType,
    zero_center: bool,
    svd_solver: SVDSolver,
    warn_pat_expected: str | None,
):
    A = array_type(A_list).astype("float32")
    adata = AnnData(A)

    if warn_pat_expected is not None:
        with pytest.warns((UserWarning, FutureWarning), match=warn_pat_expected):  # noqa: PT031
            warnings.filterwarnings(
                "ignore", r".*Using a dense eigensolver instead of LOBPCG", UserWarning
            )
            sc.pp.pca(adata, svd_solver=svd_solver, zero_center=zero_center)
        return

    warnings.simplefilter("error")
    sc.pp.pca(adata, svd_solver=svd_solver, zero_center=zero_center)


def test_pca_transform(array_type):
    adata = AnnData(array_type(A_list).astype("float32"))
    A_pca_abs = np.abs(A_pca)

    warnings.filterwarnings("error")
    sc.pp.pca(adata, n_comps=4, zero_center=True, dtype="float64")

    adata = to_memory(adata)
    assert np.linalg.norm(A_pca_abs[:, :4] - np.abs(adata.obsm["X_pca"])) < 2e-05


def test_pca_transform_randomized(array_type):
    adata = AnnData(array_type(A_list).astype("float32"))
    A_pca_abs = np.abs(A_pca)

    warnings.filterwarnings("error")
    if isinstance(adata.X, DaskArray) and isinstance(adata.X._meta, CSBase):
        patterns = (
            r"Ignoring random_state=14 when using a sparse dask array",
            r"Ignoring svd_solver='randomized' when using a sparse dask array",
        )
        ctx = _helpers.MultiContext(
            *(pytest.warns(UserWarning, match=pattern) for pattern in patterns)
        )
    elif isinstance(adata.X, CSBase):
        ctx = pytest.warns(UserWarning, match=r"Ignoring.*'randomized")
    else:
        ctx = nullcontext()

    with ctx:
        sc.pp.pca(
            adata,
            n_comps=4,
            zero_center=True,
            svd_solver="randomized",
            dtype="float64",
            random_state=14,
        )

    assert np.linalg.norm(A_pca_abs[:, :4] - np.abs(adata.obsm["X_pca"])) < 2e-05


def test_pca_transform_no_zero_center(request: pytest.FixtureRequest, array_type):
    adata = AnnData(array_type(A_list).astype("float32"))
    A_svd_abs = np.abs(A_svd)
    if isinstance(adata.X, DaskArray) and isinstance(adata.X._meta, CSBase):
        reason = "TruncatedSVD is not supported for sparse Dask yet"
        request.applymarker(pytest.mark.xfail(reason=reason))

    warnings.filterwarnings("error")
    sc.pp.pca(adata, n_comps=4, zero_center=False, dtype="float64", random_state=14)

    assert np.linalg.norm(A_svd_abs[:, :4] - np.abs(adata.obsm["X_pca"])) < 2e-05


def test_pca_shapes():
    """Tests that n_comps behaves correctly.

    See <https://github.com/scverse/scanpy/issues/1051>
    """
    adata = AnnData(np.random.randn(30, 20))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (30, 19)

    adata = AnnData(np.random.randn(20, 30))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (20, 19)

    with pytest.raises(
        ValueError,
        match=r"n_components=100 must be between 1 and.*20 with svd_solver='arpack'",
    ):
        sc.pp.pca(adata, n_comps=100)


@pytest.mark.parametrize(
    ("key_added", "keys_expected"),
    [
        pytest.param(None, ("X_pca", "PCs", "pca"), id="None"),
        pytest.param("custom_key", ("custom_key",) * 3, id="custom_key"),
    ],
)
def test_pca_sparse(key_added: str | None, keys_expected: tuple[str, str, str]):
    """Tests implicitly centered pca on sparse arrays.

    Checks if it returns equivalent results to explicit centering on dense arrays.
    """
    pbmc = pbmc3k_normalized()[:200].copy()

    pbmc_dense = pbmc.copy()
    pbmc_dense.X = pbmc_dense.X.toarray()

    implicit = sc.pp.pca(pbmc, dtype=np.float64, copy=True)
    explicit = sc.pp.pca(pbmc_dense, dtype=np.float64, key_added=key_added, copy=True)

    key_obsm, key_varm, key_uns = keys_expected

    np.testing.assert_allclose(
        implicit.uns["pca"]["variance"], explicit.uns[key_uns]["variance"]
    )
    np.testing.assert_allclose(
        implicit.uns["pca"]["variance_ratio"], explicit.uns[key_uns]["variance_ratio"]
    )
    np.testing.assert_allclose(implicit.obsm["X_pca"], explicit.obsm[key_obsm])
    np.testing.assert_allclose(implicit.varm["PCs"], explicit.varm[key_varm])


def test_pca_reproducible(array_type):
    pbmc = pbmc3k_normalized()
    pbmc.X = array_type(pbmc.X)

    with (
        pytest.warns(UserWarning, match=r"Ignoring random_state.*sparse dask array")
        if isinstance(pbmc.X, DaskArray) and isinstance(pbmc.X._meta, CSBase)
        else nullcontext()
    ):
        a = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=42)
        b = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=42)
        c = sc.pp.pca(pbmc, copy=True, dtype=np.float64, random_state=0)

    assert_equal(a, b)

    # Test that changing random seed changes result
    # Does not show up reliably with 32 bit computation
    # sparse-in-dask doesn’t use a random seed, so it also doesn’t work there.
    if not (isinstance(pbmc.X, DaskArray) and isinstance(pbmc.X._meta, CSBase)):
        a, c = map(to_memory, [a, c])
        assert not np.array_equal(a.obsm["X_pca"], c.obsm["X_pca"])


def test_pca_chunked():
    """Tests that chunked PCA is equivalent to default PCA.

    See also <https://github.com/scverse/scanpy/issues/1590>
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
    """Tests that the n_pcs parameter also works for representations not called "X_pca"."""
    pbmc = pbmc3k_normalized()
    sc.pp.pca(pbmc, dtype=np.float64)
    pbmc.obsm["X_pca_test"] = pbmc.obsm["X_pca"]
    original = sc.pp.neighbors(pbmc, n_pcs=5, use_rep="X_pca", copy=True)
    renamed = sc.pp.neighbors(pbmc, n_pcs=5, use_rep="X_pca_test", copy=True)

    assert np.allclose(original.obsm["X_pca"], renamed.obsm["X_pca_test"])
    assert np.allclose(
        original.obsp["distances"].toarray(), renamed.obsp["distances"].toarray()
    )


# We use all possible array types here since this error should be raised before
# PCA can realize that it got a Dask array
@pytest.mark.parametrize("array_type", ARRAY_TYPES_ALL)
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
    mask_var = _helpers.random_mask(adata.shape[1] + 1)
    with pytest.raises(
        ValueError, match=r"The shape of the mask do not match the data\."
    ):
        sc.pp.pca(adata, mask_var=mask_var, copy=True)


def test_mask_var_argument_equivalence(float_dtype, array_type):
    """Test if pca result is equal when given mask as boolarray vs string."""
    adata_base = AnnData(array_type(np.random.random((100, 10))).astype(float_dtype))
    mask_var = _helpers.random_mask(adata_base.shape[1])

    adata = adata_base.copy()
    sc.pp.pca(adata, mask_var=mask_var, dtype=float_dtype)

    adata_w_mask = adata_base.copy()
    adata_w_mask.var["mask"] = mask_var
    sc.pp.pca(adata_w_mask, mask_var="mask", dtype=float_dtype)

    adata, adata_w_mask = map(to_memory, [adata, adata_w_mask])
    assert np.allclose(
        adata.X.toarray() if isinstance(adata.X, CSBase) else adata.X,
        adata_w_mask.X.toarray()
        if isinstance(adata_w_mask.X, CSBase)
        else adata_w_mask.X,
    )


def test_mask(request: pytest.FixtureRequest, array_type):
    if array_type in DASK_CONVERTERS.values():
        reason = "TODO: Dask arrays are not supported"
        request.applymarker(pytest.mark.xfail(reason=reason))
    adata = sc.datasets.blobs(n_variables=10, n_centers=3, n_observations=100)
    adata.X = array_type(adata.X)

    if isinstance(adata.X, np.ndarray) and Version(ad.__version__) < Version("0.9"):
        reason = (
            "TODO: Previous version of anndata would return an F ordered array for one"
            " case here, which surprisingly considerably changes the results of PCA."
        )
        request.applymarker(pytest.mark.xfail(reason=reason))
    mask_var = _helpers.random_mask(adata.shape[1])

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


def test_mask_order_warning(request: pytest.FixtureRequest):
    if Version(ad.__version__) >= Version("0.9"):
        reason = "Not expected to warn in later versions of anndata"
        request.applymarker(pytest.mark.xfail(reason=reason))

    adata = ad.AnnData(X=np.random.randn(50, 5))
    mask = np.array([True, False, True, False, True])

    with pytest.warns(
        UserWarning,
        match="When using a mask parameter with anndata<0.9 on a dense array",
    ):
        sc.pp.pca(adata, mask_var=mask)


def test_mask_defaults(array_type, float_dtype):
    """Test if PCA behavior in relation to highly variable genes.

    1. That it’s equal withwithout and with – but mask is None
    2. If pca takes highly variable as mask as default
    """
    A = array_type(A_list).astype("float64")
    adata = AnnData(A)

    without_var = sc.pp.pca(adata, copy=True, dtype=float_dtype)

    rng = np.random.default_rng(8)
    mask = _helpers.random_mask(adata.shape[1], rng=rng)
    adata.var["highly_variable"] = mask
    with_var = sc.pp.pca(adata, copy=True, dtype=float_dtype)
    assert without_var.uns["pca"]["params"]["mask_var"] is None
    assert with_var.uns["pca"]["params"]["mask_var"] == "highly_variable"
    without_var, with_var = map(to_memory, [without_var, with_var])
    assert not np.array_equal(without_var.obsm["X_pca"], with_var.obsm["X_pca"])

    with_no_mask = sc.pp.pca(adata, mask_var=None, copy=True, dtype=float_dtype)
    with_no_mask = to_memory(with_no_mask)
    assert np.array_equal(without_var.obsm["X_pca"], with_no_mask.obsm["X_pca"])


def test_pca_layer():
    """Tests that layers works the same way as `X`."""
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


# Skipping these tests during min-deps testing shouldn't be an issue because the sparse-in-dask feature is not available on anndata<0.10 anyway
needs_anndata_dask = pytest.mark.skipif(
    pkg_version("anndata") < Version("0.10"),
    reason="Old AnnData doesn’t have dask test helpers",
)


@needs.dask
@needs_anndata_dask
@pytest.mark.parametrize(
    "other_array_type",
    [
        lambda x: x.toarray(),
        DASK_CONVERTERS[_helpers.as_sparse_dask_array],
        DASK_CONVERTERS[_helpers.as_dense_dask_array],
    ],
    ids=["dense-mem", "sparse-dask", "dense-dask"],
)
def test_covariance_eigh_impls(other_array_type):
    warnings.filterwarnings("error")

    adata_sparse_mem = pbmc3k_normalized()[:200, :100].copy()
    adata_other = adata_sparse_mem.copy()
    adata_other.X = other_array_type(adata_other.X)

    sc.pp.pca(adata_sparse_mem, svd_solver="covariance_eigh")
    sc.pp.pca(adata_other, svd_solver="covariance_eigh")

    to_memory(adata_other)
    np.testing.assert_allclose(
        np.abs(adata_sparse_mem.obsm["X_pca"]), np.abs(adata_other.obsm["X_pca"])
    )


@needs.dask
@needs_anndata_dask
@pytest.mark.parametrize(
    ("msg_re", "op"),
    [
        (
            r"Only sparse dask arrays with CSR-meta",
            lambda a: a.map_blocks(
                sparse.csc_matrix,  # noqa: TID251
                meta=sparse.csc_matrix(np.array([])),  # noqa: TID251
            ),
        ),
        (r"Only dask arrays with chunking", lambda a: a.rechunk((a.shape[0], 100))),
        (
            r"Only dask arrays with chunking",
            lambda a: a.map_blocks(np.array, meta=np.array([])).rechunk(
                (a.shape[0], 100)
            ),
        ),
    ],
    ids=["as-csc", "bad-chunking", "bad-chunking-dense"],
)
def test_sparse_dask_input_errors(msg_re: str, op: Callable[[DaskArray], DaskArray]):
    adata_sparse = pbmc3k_normalized()
    adata_sparse.X = op(DASK_CONVERTERS[_helpers.as_sparse_dask_array](adata_sparse.X))

    with pytest.raises(ValueError, match=msg_re):
        sc.pp.pca(adata_sparse, svd_solver="covariance_eigh")


@needs.dask
@needs_anndata_dask
@pytest.mark.parametrize(
    ("dtype", "dtype_arg", "rtol"),
    [
        pytest.param(np.float32, None, 1e-5, id="float32"),
        pytest.param(np.float32, np.float64, None, id="float32-float64"),
        pytest.param(np.float64, None, None, id="float64"),
        pytest.param(np.int64, None, None, id="int64"),
    ],
)
def test_cov_sparse_dask(dtype, dtype_arg, rtol):
    x_arr = A_list.astype(dtype)
    x = DASK_CONVERTERS[_helpers.as_sparse_dask_array](x_arr)
    cov, gram, mean = _cov_sparse_dask(x, return_gram=True, dtype=dtype_arg)
    np.testing.assert_allclose(mean, np.mean(x_arr, axis=0))
    np.testing.assert_allclose(gram, (x_arr.T @ x_arr) / x.shape[0])
    tol_args = dict(rtol=rtol) if rtol is not None else {}
    np.testing.assert_allclose(cov, np.cov(x_arr, rowvar=False, bias=True), **tol_args)

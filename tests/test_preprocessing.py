from __future__ import annotations

import warnings
from contextlib import nullcontext
from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, ImplicitModificationWarning
from anndata.tests.helpers import asarray, assert_equal
from fast_array_utils import conv
from numpy.testing import assert_allclose
from packaging.version import Version
from scipy import sparse

import scanpy as sc
from scanpy._compat import CSBase, pkg_version
from testing.scanpy._helpers import (
    check_rep_mutation,
    check_rep_results,
    maybe_dask_process_context,
)
from testing.scanpy._helpers.data import pbmc3k, pbmc68k_reduced
from testing.scanpy._pytest.params import ARRAY_TYPES, ARRAY_TYPES_SPARSE

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    from numpy.typing import DTypeLike, NDArray


class _MatrixFormat(NamedTuple):
    callback: Callable[[NDArray], CSBase | NDArray]
    is_sparse: bool

    def __call__(self, x: NDArray) -> CSBase | NDArray:
        return self.callback(x)


HERE = Path(__file__).parent
DATA_PATH = HERE / "_data"


@pytest.fixture(params=[np.asarray, sparse.csr_matrix, sparse.csc_matrix])  # noqa: TID251
def count_matrix_format(request: pytest.FixtureRequest) -> _MatrixFormat:
    is_sparse = isinstance(request.param, type) and issubclass(request.param, CSBase)
    return _MatrixFormat(request.param, is_sparse=is_sparse)


@pytest.fixture(params=[np.int64, np.float32, np.float64])
def dtype(request: pytest.FixtureRequest) -> DTypeLike:
    return request.param


@pytest.fixture(params=[True, False], ids=["zero_center", "no_zero_center"])
def zero_center(request: pytest.FixtureRequest) -> bool:
    return request.param


def test_log1p(tmp_path):
    A = np.random.rand(200, 10).astype(np.float32)
    A_l = np.log1p(A)
    ad = AnnData(A.copy())
    ad2 = AnnData(A.copy())
    ad3 = AnnData(A.copy())
    ad3.filename = tmp_path / "test.h5ad"
    sc.pp.log1p(ad)
    assert np.allclose(ad.X, A_l)
    sc.pp.log1p(ad2, chunked=True)
    assert np.allclose(ad2.X, ad.X)
    sc.pp.log1p(ad3, chunked=True)
    assert np.allclose(ad3.X, ad.X)

    # Test base
    ad4 = AnnData(A)
    sc.pp.log1p(ad4, base=2)
    assert np.allclose(ad4.X, A_l / np.log(2))


def test_log1p_deprecated_arg():
    A = np.random.rand(200, 10).astype(np.float32)
    with pytest.warns(FutureWarning, match=r".*`X` was renamed to `data`"):
        sc.pp.log1p(X=A)


@pytest.fixture(params=[None, 2])
def base(request):
    return request.param


def test_log1p_rep(count_matrix_format: _MatrixFormat, base, dtype: DTypeLike) -> None:
    X = count_matrix_format(
        np.abs(sparse.random(100, 200, density=0.3, dtype=dtype)).toarray()
    )
    check_rep_mutation(sc.pp.log1p, X, base=base)
    check_rep_results(sc.pp.log1p, X, base=base)


def test_normalize_per_cell() -> None:
    A = np.array([[1, 0], [3, 0], [5, 6]], dtype=np.float32)
    adata = AnnData(A.copy())
    with pytest.warns(FutureWarning, match=r"sc\.pp\.normalize_total"):
        sc.pp.normalize_per_cell(
            adata, counts_per_cell_after=1, key_n_counts="n_counts2"
        )
    assert adata.X.sum(axis=1).tolist() == [1.0, 1.0, 1.0]
    # now with copy option
    adata = AnnData(A.copy())
    # note that sc.pp.normalize_per_cell is also used in
    # pl.highest_expr_genes with parameter counts_per_cell_after=100
    with pytest.warns(FutureWarning, match=r"sc\.pp\.normalize_total"):
        adata_copy = sc.pp.normalize_per_cell(adata, counts_per_cell_after=1, copy=True)
    assert adata_copy.X.sum(axis=1).tolist() == [1.0, 1.0, 1.0]
    # now sparse
    adata = AnnData(A.copy())
    adata_sparse = AnnData(sparse.csr_matrix(A.copy()))  # noqa: TID251
    with pytest.warns(FutureWarning, match=r"sc\.pp\.normalize_total"):
        sc.pp.normalize_per_cell(adata)
    with pytest.warns(FutureWarning, match=r"sc\.pp\.normalize_total"):
        sc.pp.normalize_per_cell(adata_sparse)
    assert adata.X.sum(axis=1).tolist() == adata_sparse.X.sum(axis=1).A1.tolist()


def _random_probs(n: int, frac_zero: float) -> NDArray[np.float64]:
    """Generate a random probability distribution of `n` values between 0 and 1."""
    probs = np.random.randint(0, 10000, n).astype(np.float64)
    probs[probs < np.quantile(probs, frac_zero)] = 0
    probs /= probs.sum()
    np.testing.assert_almost_equal(probs.sum(), 1)
    return probs


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("which", ["copy", "inplace", "array"])
@pytest.mark.parametrize(
    ("axis", "f_or_n", "replace"),
    [
        pytest.param(0, 40, False, id="obs-40-no_replace"),
        pytest.param(0, 0.1, False, id="obs-0.1-no_replace"),
        pytest.param(0, 201, True, id="obs-201-replace"),
        pytest.param(0, 1, True, id="obs-1-replace"),
        pytest.param(1, 10, False, id="var-10-no_replace"),
        pytest.param(1, 11, True, id="var-11-replace"),
        pytest.param(1, 2.0, True, id="var-2.0-replace"),
    ],
)
@pytest.mark.parametrize(
    "ps",
    [
        dict(obs=None, var=None),
        dict(obs=np.tile([True, False], 100), var=np.tile([True, False], 5)),
        dict(obs=_random_probs(200, 0.3), var=_random_probs(10, 0.7)),
    ],
    ids=["all", "mask", "p"],
)
def test_sample(
    *,
    request: pytest.FixtureRequest,
    array_type: Callable[[np.ndarray], np.ndarray | CSBase],
    which: Literal["copy", "inplace", "array"],
    axis: Literal[0, 1],
    f_or_n: float | int,  # noqa: PYI041
    replace: bool,
    ps: dict[Literal["obs", "var"], NDArray[np.bool_] | None],
):
    adata = AnnData(array_type(np.ones((200, 10))))
    p = ps["obs" if axis == 0 else "var"]
    expected = int(adata.shape[axis] * f_or_n) if isinstance(f_or_n, float) else f_or_n
    if p is not None and not replace and expected > (n_possible := (p != 0).sum()):
        request.applymarker(pytest.xfail(f"Can’t draw {expected} out of {n_possible}"))

    # ignoring this warning declaratively is a pain so do it here
    if find_spec("dask"):
        import dask.array as da

        warnings.filterwarnings("ignore", category=da.PerformanceWarning)
    # can’t guarantee that duplicates are drawn when `replace=True`,
    # so we just ignore the warning instead using `with pytest.warns(...)`
    warnings.filterwarnings(
        "ignore" if replace else "error", r".*names are not unique", UserWarning
    )
    rv = sc.pp.sample(
        adata.X if which == "array" else adata,
        f_or_n if isinstance(f_or_n, float) else None,
        n=f_or_n if isinstance(f_or_n, int) else None,
        replace=replace,
        axis=axis,
        # `copy` only effects AnnData inputs
        copy=dict(copy=True, inplace=False, array=False)[which],
        p=p,
    )

    match which:
        case "copy":
            subset = rv
            assert rv is not adata
            assert adata.shape == (200, 10)
        case "inplace":
            subset = adata
            assert rv is None
        case "array":
            subset, indices = rv
            assert len(indices) == expected
            assert adata.shape == (200, 10)
        case _:
            pytest.fail(f"Unknown `{which=}`")

    assert subset.shape == ((expected, 10) if axis == 0 else (200, expected))


@pytest.mark.parametrize(
    ("args", "exc", "pattern"),
    [
        pytest.param(
            dict(), TypeError, r"Either `fraction` or `n` must be set", id="empty"
        ),
        pytest.param(
            dict(n=10, fraction=0.2),
            TypeError,
            r"Providing both `fraction` and `n` is not allowed",
            id="both",
        ),
        pytest.param(
            dict(fraction=2),
            ValueError,
            r"If `replace=False`, `fraction=2` needs to be",
            id="frac>1",
        ),
        pytest.param(
            dict(fraction=-0.3),
            ValueError,
            r"`fraction=-0\.3` needs to be nonnegative",
            id="frac<0",
        ),
        pytest.param(
            dict(n=3, p=np.ones(200, dtype=np.int32)),
            ValueError,
            r"mask/probabilities array must be boolean or floating point",
            id="type(p)",
        ),
    ],
)
def test_sample_error(args: dict[str, Any], exc: type[Exception], pattern: str):
    adata = AnnData(np.ones((200, 10)))
    with pytest.raises(exc, match=pattern):
        sc.pp.sample(adata, **args)


def test_sample_backwards_compat():
    expected = np.array(
        [26, 86, 2, 55, 75, 93, 16, 73, 54, 95, 53, 92, 78, 13, 7, 30, 22, 24, 33, 8]
    )
    legacy_result, indices = sc.pp.subsample(np.arange(100), n_obs=20)
    assert np.array_equal(indices, legacy_result), "arange choices should match indices"
    assert np.array_equal(legacy_result, expected)


def test_sample_copy_backed(tmp_path):
    adata_m = AnnData(np.random.rand(200, 10).astype(np.float32))
    adata_d = adata_m.copy()
    adata_d.filename = tmp_path / "test.h5ad"

    assert sc.pp.sample(adata_d, n=40, copy=True).shape == (40, 10)
    np.testing.assert_array_equal(
        sc.pp.sample(adata_m, n=40, copy=True, rng=0).X,
        sc.pp.sample(adata_d, n=40, copy=True, rng=0).X,
    )


def test_sample_copy_backed_error(tmp_path):
    adata_d = AnnData(np.random.rand(200, 10).astype(np.float32))
    adata_d.filename = tmp_path / "test.h5ad"
    with pytest.raises(NotImplementedError):
        sc.pp.sample(adata_d, n=40, copy=False)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("max_value", [None, 1.0], ids=["no_clip", "clip"])
def test_scale_matrix_types(array_type, zero_center, max_value):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    warn_ctx = pytest.warns(UserWarning, match=r"zero-centering.*densifies")
    with warn_ctx if zero_center else nullcontext():
        sc.pp.scale(adata, zero_center=zero_center, max_value=max_value)
        adata.X = conv.to_dense(adata.X)
    with (
        (
            warn_ctx
            if zero_center
            and any(pat in array_type.__name__ for pat in ("sparse", "csc", "csr"))
            else nullcontext()
        ),
        maybe_dask_process_context(),
    ):
        sc.pp.scale(adata_casted, zero_center=zero_center, max_value=max_value)
        if is_dask := ("dask" in array_type.__name__):
            assert not isinstance(adata_casted.X._meta, np.matrix)
        adata_casted.X = conv.to_dense(adata_casted.X, to_cpu_memory=True)
    assert_allclose(
        adata_casted.X,
        adata.X,
        rtol=1e-1 if is_dask else 1e-5,
        atol=1e-1 if is_dask else 1e-5,
    )


@pytest.mark.parametrize("array_type", ARRAY_TYPES_SPARSE)
def test_scale_zero_center_warns_dask_sparse(array_type):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    with pytest.warns(UserWarning, match=r"zero-center.*densifies"):
        sc.pp.scale(adata_casted)
    with pytest.warns(UserWarning, match=r"zero-center.*densifies"):
        sc.pp.scale(adata)
    assert_allclose(adata_casted.X, adata.X, rtol=1e-5, atol=1e-5)


def test_scale():
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    v = adata[:, 0 : adata.shape[1] // 2]
    # Should turn view to copy https://github.com/scverse/anndata/issues/171#issuecomment-508689965
    assert v.is_view
    with (
        pytest.warns(UserWarning, match=r"zero-center.*densifies"),
        pytest.warns(UserWarning, match=r"Received a view"),
    ):
        sc.pp.scale(v)
    assert not v.is_view
    assert_allclose(v.X.var(axis=0), np.ones(v.shape[1]), atol=0.01)
    assert_allclose(v.X.mean(axis=0), np.zeros(v.shape[1]), atol=0.00001)


def test_scale_rep(*, count_matrix_format: _MatrixFormat, zero_center: bool) -> None:
    """Test that it doesn't matter where the array being scaled is in the anndata object."""
    X = count_matrix_format(sparse.random(100, 200, density=0.3).toarray())
    ctx = (
        pytest.warns(UserWarning, match=r"zero-center.*densifies")
        if zero_center and count_matrix_format.is_sparse
        else nullcontext()
    )
    with ctx:
        check_rep_mutation(sc.pp.scale, X, zero_center=zero_center)
    with ctx:
        check_rep_results(sc.pp.scale, X, zero_center=zero_center)


def test_scale_array(*, count_matrix_format: _MatrixFormat, zero_center: bool) -> None:
    """Test that running sc.pp.scale on an anndata object and an array returns the same results."""
    X = count_matrix_format(sparse.random(100, 200, density=0.3).toarray())
    adata = AnnData(X=X.copy())

    ctx = (
        pytest.warns(UserWarning, match=r"zero-center.*densifies")
        if zero_center and count_matrix_format.is_sparse
        else nullcontext()
    )
    with ctx:
        sc.pp.scale(adata, zero_center=zero_center)
    with ctx:
        scaled_X = sc.pp.scale(X, zero_center=zero_center, copy=True)
    np.testing.assert_equal(asarray(scaled_X), asarray(adata.X))


# https://github.com/pandas-dev/pandas/issues/61928
@pytest.mark.filterwarnings("ignore:invalid value encountered in cast:RuntimeWarning")
def test_recipe_plotting() -> None:
    sc.settings.autoshow = False
    adata = AnnData(np.random.randint(0, 1000, (1000, 1000)))
    if pkg_version("pandas") < Version("2.2"):
        # https://github.com/pandas-dev/pandas/issues/54661
        warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
    # These shouldn't throw an error
    with pytest.warns(FutureWarning, match=r"sc\.p[pl]\.highly_variable_genes"):
        sc.pp.recipe_seurat(adata.copy(), plot=True)
    with pytest.warns(FutureWarning, match=r"sc\.p[pl]\.highly_variable_genes"):
        sc.pp.recipe_zheng17(adata.copy(), plot=True)


def test_regress_out_ordinal():
    from scipy.sparse import random

    adata = AnnData(random(1000, 100, density=0.6, format="csr"))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    # results using only one processor
    single = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )
    assert adata.X.shape == single.X.shape

    # results using 8 processors
    multi = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=8, copy=True
    )

    np.testing.assert_array_equal(single.X, multi.X)


@pytest.mark.parametrize("dtype", [np.uint32, np.float64, np.uint64])
def test_regress_out_int(dtype):
    adata = pbmc3k()[:200, :200].copy()
    adata.X = adata.X.astype(np.float64 if dtype != np.uint32 else np.float32)
    dtype = adata.X.dtype
    adata.obs["labels"] = pd.Categorical(
        (["A"] * (adata.X.shape[0] - 100)) + (["B"] * 100)
    )
    adata_other = adata.copy()
    adata_other.X = adata_other.X.astype(dtype)
    # results using only one processor
    sc.pp.regress_out(adata, keys=["labels"])
    sc.pp.regress_out(adata_other, keys=["labels"])
    assert_equal(adata_other, adata)
    # This file was generated under scanpy 1.10.3
    ground_truth = np.load(DATA_PATH / "cat_regressor_for_int_input.npy")
    np.testing.assert_allclose(ground_truth, adata_other.X, atol=1e-5, rtol=1e-5)


@pytest.mark.parametrize("dtype", [np.int64, np.float64, np.int32])
def test_regress_out_layer(dtype):
    from scipy.sparse import random

    adata = AnnData(
        random(1000, 100, density=0.6, format="csr", dtype=np.uint16).astype(dtype)
    )
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    if dtype == np.float64:
        dtype_cast = dtype
    if dtype == np.int64:
        dtype_cast = np.float64
    if dtype == np.int32:
        dtype_cast = np.float32
    adata.layers["counts"] = adata.X.copy().astype(dtype_cast)

    single = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )
    assert adata.X.shape == single.X.shape

    layer = sc.pp.regress_out(
        adata, layer="counts", keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )

    np.testing.assert_allclose(single.X, layer.layers["counts"])


def test_regress_out_view():
    from scipy.sparse import random

    adata = AnnData(random(500, 1100, density=0.2, format="csr"))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    subset_adata = adata[:, :1050]
    subset_adata_copy = subset_adata.copy()
    with pytest.warns(UserWarning, match=r"Received a view"):
        sc.pp.regress_out(subset_adata, keys=["n_counts", "percent_mito"])
    sc.pp.regress_out(subset_adata_copy, keys=["n_counts", "percent_mito"])
    assert_equal(subset_adata, subset_adata_copy)
    assert not subset_adata.is_view


def test_regress_out_categorical():
    import pandas as pd
    from scipy.sparse import random

    adata = AnnData(random(1000, 100, density=0.6, format="csr"))
    # create a categorical column
    adata.obs["batch"] = pd.Categorical(np.random.randint(1, 4, size=adata.X.shape[0]))

    multi = sc.pp.regress_out(adata, keys="batch", n_jobs=8, copy=True)
    assert adata.X.shape == multi.X.shape


def test_regress_out_constants():
    adata = AnnData(np.hstack((np.full((10, 1), 0.0), np.full((10, 1), 1.0))))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata_copy = adata.copy()

    sc.pp.regress_out(adata, keys=["n_counts", "percent_mito"])
    assert_equal(adata, adata_copy)


@pytest.mark.parametrize(
    ("keys", "test_file", "atol"),
    [
        (["n_counts", "percent_mito"], "regress_test_small.npy", 0.0),
        (["bulk_labels"], "regress_test_small_cat.npy", 1e-6),
    ],
)
def test_regress_out_reproducible(keys, test_file, atol):
    adata = sc.datasets.pbmc68k_reduced()
    adata = adata.raw.to_adata()[:200, :200].copy()
    sc.pp.regress_out(adata, keys=keys)
    # This file was generated from the original implementation in version 1.10.3
    # Now we compare new implementation with the old one
    tester = np.load(DATA_PATH / test_file)
    np.testing.assert_allclose(adata.X, tester, atol=atol)


def test_regress_out_constants_equivalent():
    # Tests that constant values don't change results
    # (since support for constant values is implemented by us)
    from sklearn.datasets import make_blobs

    X, cat = make_blobs(100, 20)
    a = sc.AnnData(np.hstack([X, np.zeros((100, 5))]), obs={"cat": pd.Categorical(cat)})
    b = sc.AnnData(X, obs={"cat": pd.Categorical(cat)})

    sc.pp.regress_out(a, "cat")
    sc.pp.regress_out(b, "cat")

    np.testing.assert_equal(a[:, b.var_names].X, b.X)


@pytest.mark.parametrize("replace", [True, False], ids=["replace", "no_replace"])
def test_downsample_counts_per_cell(
    *, count_matrix_format: _MatrixFormat, replace: bool, dtype: DTypeLike
) -> None:
    TARGET = 1000
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata = AnnData(X=count_matrix_format(X).astype(dtype))
    with pytest.raises(ValueError, match=r"Must specify exactly one"):
        sc.pp.downsample_counts(
            adata, counts_per_cell=TARGET, total_counts=TARGET, replace=replace
        )
    with pytest.raises(ValueError, match=r"Must specify exactly one"):
        sc.pp.downsample_counts(adata, replace=replace)
    initial_totals = np.ravel(adata.X.sum(axis=1))
    adata = sc.pp.downsample_counts(
        adata, counts_per_cell=TARGET, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if isinstance(adata.X, CSBase):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGET)
    assert all(initial_totals >= new_totals)
    assert all(
        initial_totals[initial_totals <= TARGET] == new_totals[initial_totals <= TARGET]
    )
    if not replace:
        assert np.all(X >= adata.X)
    assert X.dtype == adata.X.dtype


@pytest.mark.parametrize("replace", [True, False], ids=["replace", "no_replace"])
def test_downsample_counts_per_cell_multiple_targets(
    *, count_matrix_format: _MatrixFormat, replace: bool, dtype: DTypeLike
) -> None:
    TARGETS = np.random.randint(500, 1500, 1000)
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata = AnnData(X=count_matrix_format(X).astype(dtype))
    initial_totals = np.ravel(adata.X.sum(axis=1))
    with pytest.raises(ValueError, match=r"counts_per_cell.*length as number of obs"):
        sc.pp.downsample_counts(adata, counts_per_cell=[40, 10], replace=replace)
    adata = sc.pp.downsample_counts(
        adata, counts_per_cell=TARGETS, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if isinstance(adata.X, CSBase):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGETS)
    assert all(initial_totals >= new_totals)
    assert all(
        initial_totals[initial_totals <= TARGETS]
        == new_totals[initial_totals <= TARGETS]
    )
    if not replace:
        assert np.all(X >= adata.X)
    assert X.dtype == adata.X.dtype


@pytest.mark.parametrize("replace", [True, False], ids=["replace", "no_replace"])
def test_downsample_total_counts(
    *, count_matrix_format: _MatrixFormat, replace: bool, dtype: DTypeLike
) -> None:
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata_orig = AnnData(X=count_matrix_format(X))
    total = X.sum()
    target = np.floor_divide(total, 10)
    initial_totals = np.ravel(adata_orig.X.sum(axis=1))
    adata = sc.pp.downsample_counts(
        adata_orig, total_counts=target, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if isinstance(adata.X, CSBase):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert adata.X.sum() == target
    assert all(initial_totals >= new_totals)
    if not replace:
        assert np.all(X >= adata.X)
        adata = sc.pp.downsample_counts(
            adata_orig, total_counts=total + 10, replace=False, copy=True
        )
        assert (adata.X == X).all()
    assert X.dtype == adata.X.dtype


def test_recipe_weinreb():
    # Just tests for failure for now
    adata = pbmc68k_reduced().raw.to_adata()
    adata.X = adata.X.toarray()

    orig = adata.copy()
    sc.pp.recipe_weinreb17(adata, log=False, copy=True)
    assert_equal(orig, adata)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize(
    ("max_cells", "max_counts", "min_cells", "min_counts"),
    [
        (100, None, None, None),
        (None, 100, None, None),
        (None, None, 20, None),
        (None, None, None, 20),
    ],
)
def test_filter_genes(array_type, max_cells, max_counts, min_cells, min_counts):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    sc.pp.filter_genes(
        adata,
        max_cells=max_cells,
        max_counts=max_counts,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    sc.pp.filter_genes(
        adata_casted,
        max_cells=max_cells,
        max_counts=max_counts,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    adata_casted.X = conv.to_dense(adata_casted.X, to_cpu_memory=True)
    adata.X = conv.to_dense(adata.X)
    assert_allclose(adata_casted.X, adata.X, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize(
    ("max_genes", "max_counts", "min_genes", "min_counts"),
    [
        pytest.param(100, None, None, None, id="max_genes"),
        pytest.param(None, 100, None, None, id="max_counts"),
        pytest.param(None, None, 20, None, id="min_genes"),
        pytest.param(None, None, None, 20, id="min_counts"),
    ],
)
def test_filter_cells(array_type, max_genes, max_counts, min_genes, min_counts):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    if pkg_version("pandas") < Version("2.2"):
        # https://github.com/pandas-dev/pandas/issues/54661
        warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
    sc.pp.filter_cells(
        adata,
        max_genes=max_genes,
        max_counts=max_counts,
        min_genes=min_genes,
        min_counts=min_counts,
    )
    sc.pp.filter_cells(
        adata_casted,
        max_genes=max_genes,
        max_counts=max_counts,
        min_genes=min_genes,
        min_counts=min_counts,
    )
    adata_casted.X = conv.to_dense(adata_casted.X, to_cpu_memory=True)
    adata.X = conv.to_dense(adata.X)
    assert_allclose(adata_casted.X, adata.X, rtol=1e-5, atol=1e-5)

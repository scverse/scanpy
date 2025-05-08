from __future__ import annotations

import itertools
import string
from operator import mul, truediv
from types import ModuleType
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata.tests.helpers import asarray
from packaging.version import Version
from scipy import sparse

from scanpy._compat import CSBase, DaskArray, pkg_version
from scanpy._utils import (
    axis_mul_or_truediv,
    axis_sum,
    check_nonnegative_integers,
    descend_classes_and_funcs,
    elem_mul,
    is_constant,
)
from scanpy._utils.random import (
    ith_k_tuple,
    legacy_numpy_gen,
    random_k_tuples,
    random_str,
)
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import (
    ARRAY_TYPES,
    ARRAY_TYPES_DASK,
    ARRAY_TYPES_SPARSE,
    ARRAY_TYPES_SPARSE_DASK_UNSUPPORTED,
)

if TYPE_CHECKING:
    from typing import Any


def test_descend_classes_and_funcs():
    # create module hierarchy
    a = ModuleType("a")
    a.b = ModuleType("a.b")

    # populate with classes
    a.A = type("A", (), {})
    a.A.__module__ = a.__name__
    a.b.B = type("B", (), {})
    a.b.B.__module__ = a.b.__name__

    # create a loop to check if that gets caught
    a.b.a = a

    assert {a.A, a.b.B} == set(descend_classes_and_funcs(a, "a"))


def test_axis_mul_or_truediv_badop():
    dividend = np.array([[0, 1.0, 1.0], [1.0, 0, 1.0]])
    divisor = np.array([0.1, 0.2])
    with pytest.raises(ValueError, match=".*not one of truediv or mul"):
        axis_mul_or_truediv(dividend, divisor, op=np.add, axis=0)


def test_axis_mul_or_truediv_bad_out():
    dividend = sparse.csr_matrix(np.array([[0, 1.0, 1.0], [1.0, 0, 1.0]]))  # noqa: TID251
    divisor = np.array([0.1, 0.2])
    with pytest.raises(ValueError, match="`out` argument provided but not equal to X"):
        axis_mul_or_truediv(dividend, divisor, op=truediv, out=dividend.copy(), axis=0)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("op", [truediv, mul])
def test_scale_row(array_type, op):
    dividend = array_type(asarray([[0, 1.0, 1.0], [1.0, 0, 1.0]]))
    divisor = np.array([0.1, 0.2])
    if op is mul:
        divisor = 1 / divisor
    expd = np.array([[0, 10.0, 10.0], [5.0, 0, 5.0]])
    out = dividend if isinstance(dividend, CSBase | np.ndarray) else None
    res = asarray(axis_mul_or_truediv(dividend, divisor, op=op, axis=0, out=out))
    np.testing.assert_array_equal(res, expd)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("op", [truediv, mul])
def test_scale_column(array_type, op):
    dividend = array_type(asarray([[0, 1.0, 2.0], [3.0, 0, 4.0]]))
    divisor = np.array([0.1, 0.2, 0.5])
    if op is mul:
        divisor = 1 / divisor
    expd = np.array([[0, 5.0, 4.0], [30.0, 0, 8.0]])
    out = dividend if isinstance(dividend, CSBase | np.ndarray) else None
    res = asarray(axis_mul_or_truediv(dividend, divisor, op=op, axis=1, out=out))
    np.testing.assert_array_equal(res, expd)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_divide_by_zero(array_type):
    dividend = array_type(asarray([[0, 1.0, 2.0], [3.0, 0, 4.0]]))
    divisor = np.array([0.1, 0.2, 0.0])
    expd = np.array([[0, 5.0, 2.0], [30.0, 0, 4.0]])
    res = asarray(
        axis_mul_or_truediv(
            dividend, divisor, op=truediv, axis=1, allow_divide_by_zero=False
        )
    )
    np.testing.assert_array_equal(res, expd)
    res = asarray(
        axis_mul_or_truediv(
            dividend, divisor, op=truediv, axis=1, allow_divide_by_zero=True
        )
    )
    expd = np.array([[0, 5.0, np.inf], [30.0, 0, np.inf]])
    np.testing.assert_array_equal(res, expd)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_SPARSE)
def test_scale_out_with_dask_or_sparse_raises(array_type):
    dividend = array_type(asarray([[0, 1.0, 2.0], [3.0, 0, 4.0]]))
    divisor = np.array([0.1, 0.2, 0.5])
    if isinstance(dividend, DaskArray):
        with pytest.raises(
            TypeError if "dask" in array_type.__name__ else ValueError,
            match="`out`*",
        ):
            axis_mul_or_truediv(dividend, divisor, op=truediv, axis=1, out=dividend)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_DASK)
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("op", [truediv, mul])
def test_scale_rechunk(array_type, axis, op):
    import dask.array as da

    dividend = array_type(
        asarray([[0, 1.0, 2.0], [3.0, 0, 4.0], [3.0, 0, 4.0]])
    ).rechunk(((3,), (3,)))
    divisor = da.from_array(np.array([0.1, 0.2, 0.5]), chunks=(1,))
    if op is mul:
        divisor = 1 / divisor
    if axis == 1:
        expd = np.array([[0, 5.0, 4.0], [30.0, 0, 8.0], [30.0, 0, 8.0]])
    else:
        expd = np.array([[0, 10.0, 20.0], [15.0, 0, 20.0], [6.0, 0, 8.0]])
    out = dividend if isinstance(dividend, CSBase | np.ndarray) else None
    with pytest.warns(UserWarning, match="Rechunking scaling_array*"):
        res = asarray(axis_mul_or_truediv(dividend, divisor, op=op, axis=axis, out=out))
    np.testing.assert_array_equal(res, expd)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_elem_mul(array_type):
    m1 = array_type(asarray([[0, 1, 1], [1, 0, 1]]))
    m2 = array_type(asarray([[2, 2, 1], [3, 2, 0]]))
    expd = np.array([[0, 2, 1], [3, 0, 0]])
    res = asarray(elem_mul(m1, m2))
    np.testing.assert_array_equal(res, expd)


@pytest.mark.parametrize(
    "array_type", [*ARRAY_TYPES, pytest.param(sparse.coo_matrix, id="scipy_coo")]
)
def test_axis_sum(array_type):
    m1 = array_type(asarray([[0, 1, 1], [1, 0, 1]]))
    expd_0 = np.array([1, 1, 2])
    expd_1 = np.array([2, 2])
    res_0 = asarray(axis_sum(m1, axis=0))
    res_1 = asarray(axis_sum(m1, axis=1))
    if "matrix" in array_type.__name__:  # for sparse since dimension is kept
        res_0 = res_0.ravel()
        res_1 = res_1.ravel()
    np.testing.assert_array_equal(res_0, expd_0)
    np.testing.assert_array_equal(res_1, expd_1)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize(
    ("array_value", "expected"),
    [
        pytest.param(
            np.random.poisson(size=(100, 100)).astype(np.float64),
            True,
            id="poisson-float64",
        ),
        pytest.param(
            np.random.poisson(size=(100, 100)).astype(np.uint32),
            True,
            id="poisson-uint32",
        ),
        pytest.param(np.random.normal(size=(100, 100)), False, id="normal"),
        pytest.param(np.array([[0, 0, 0], [0, -1, 0], [0, 0, 0]]), False, id="middle"),
    ],
)
def test_check_nonnegative_integers(array_type, array_value, expected):
    X = array_type(array_value)

    received = check_nonnegative_integers(X)
    if isinstance(X, DaskArray):
        assert isinstance(received, DaskArray)
        # compute
        received = received.compute()
        assert not isinstance(received, DaskArray)
    if isinstance(received, np.bool_):
        # convert to python bool
        received = received.item()
    assert received is expected


# TODO: Make it work for sparse-in-dask
@pytest.mark.parametrize("array_type", ARRAY_TYPES_SPARSE_DASK_UNSUPPORTED)
def test_is_constant(array_type):
    constant_inds = [1, 3]
    A = np.arange(20).reshape(5, 4)
    A[constant_inds, :] = 10
    A = array_type(A)
    AT = array_type(A.T)

    assert not is_constant(A)
    assert not np.any(is_constant(A, axis=0))
    np.testing.assert_array_equal(
        [False, True, False, True, False], is_constant(A, axis=1)
    )

    assert not is_constant(AT)
    assert not np.any(is_constant(AT, axis=1))
    np.testing.assert_array_equal(
        [False, True, False, True, False], is_constant(AT, axis=0)
    )


@needs.dask
@pytest.mark.parametrize(
    ("axis", "expected"),
    [
        pytest.param(None, False, id="None"),
        pytest.param(0, [True, True, False, False], id="0"),
        pytest.param(1, [False, False, True, True, False, True], id="1"),
    ],
)
@pytest.mark.parametrize("block_type", [np.array, sparse.csr_matrix])  # noqa: TID251
def test_is_constant_dask(request: pytest.FixtureRequest, axis, expected, block_type):
    import dask.array as da

    if (
        isinstance(block_type, type)
        and issubclass(block_type, CSBase)
        and (axis is None or pkg_version("dask") < Version("2023.2.0"))
    ):
        reason = "Dask has weak support for scipy sparse matrices"
        # This test is flaky for old dask versions, but when `axis=None` it reliably fails
        request.applymarker(pytest.mark.xfail(reason=reason, strict=axis is None))

    x_data = [
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
    ]
    x = da.from_array(np.array(x_data), chunks=2).map_blocks(block_type)
    result = is_constant(x, axis=axis).compute()
    np.testing.assert_array_equal(expected, result)


@pytest.mark.parametrize("seed", [0, 1, 1256712675])
@pytest.mark.parametrize("pass_seed", [True, False], ids=["pass_seed", "set_seed"])
@pytest.mark.parametrize("func", ["choice"])
def test_legacy_numpy_gen(*, seed: int, pass_seed: bool, func: str):
    np.random.seed(seed)
    state_before = np.random.get_state(legacy=False)

    arrs: dict[bool, np.ndarray] = {}
    states_after: dict[bool, dict[str, Any]] = {}
    for direct in [True, False]:
        if not pass_seed:
            np.random.seed(seed)
        arrs[direct] = _mk_random(func, direct=direct, seed=seed if pass_seed else None)
        states_after[direct] = np.random.get_state(legacy=False)

    np.testing.assert_array_equal(arrs[True], arrs[False])
    np.testing.assert_equal(
        *states_after.values(), err_msg="both should affect global state the same"
    )
    # they should affect the global state
    with pytest.raises(AssertionError):
        np.testing.assert_equal(states_after[True], state_before)


def _mk_random(func: str, *, direct: bool, seed: int | None) -> np.ndarray:
    if direct and seed is not None:
        np.random.seed(seed)
    gen = np.random if direct else legacy_numpy_gen(seed)
    match func:
        case "choice":
            arr = np.arange(1000)
            return gen.choice(arr, size=(100, 100))
        case _:
            pytest.fail(f"Unknown {func=}")


def test_ith_k_tuple() -> None:
    """Test that the k-tuples appear in the expected order."""
    np.testing.assert_equal(
        ith_k_tuple(np.arange(2**3), n=2, k=3),
        list(itertools.product(range(2), repeat=3)),
    )


def test_random_k_tuples() -> None:
    """Test that random k-tuples are unique."""
    tups = random_k_tuples(n=26, k=6, size=10_000)
    assert tups.shape == (10_000, 6)
    assert tups.dtype == np.int64
    unique = np.unique(tups, axis=0)
    assert len(unique) == len(tups)


def test_random_str_0d() -> None:
    string = random_str(length=3, alphabet="01")
    assert string.shape == ()
    assert string.dtype == np.dtype("U3")
    assert str(string) in {"000", "001", "010", "011", "100", "101", "110", "111"}


def test_random_str() -> None:
    strings = random_str(size=26**2, length=2, alphabet=string.ascii_lowercase)
    assert strings.shape == (26**2,)
    assert strings.dtype == np.dtype("U2")
    unique = np.unique(strings, axis=0)
    assert len(unique) == len(strings)

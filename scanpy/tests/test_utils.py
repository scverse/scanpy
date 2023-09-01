from types import ModuleType
import pytest
from scipy.sparse import csr_matrix
import numpy as np

from scanpy._utils import (
    descend_classes_and_funcs,
    check_nonnegative_integers,
    is_constant,
    lazy_and,
    lazy_or,
)
from scanpy.testing._pytest.marks import needs
from scanpy._compat import DaskArray


def test_descend_classes_and_funcs():
    # create module hierarchy
    a = ModuleType('a')
    a.b = ModuleType('a.b')

    # populate with classes
    a.A = type('A', (), {})
    a.A.__module__ = a.__name__
    a.b.B = type('B', (), {})
    a.b.B.__module__ = a.b.__name__

    # create a loop to check if that gets caught
    a.b.a = a

    assert {a.A, a.b.B} == set(descend_classes_and_funcs(a, 'a'))


array_kinds = [pytest.param(DaskArray, marks=[needs('dask')], id='dask'), bool]


@pytest.mark.parametrize(
    ('lazy_op', 'left', 'right', 'expected'),
    [
        pytest.param(lazy_and, False, pytest.fail, False, id='false-and-fail'),
        pytest.param(lazy_and, True, lambda: 'hi!', 'hi!', id='true-and-pass'),
        pytest.param(lazy_or, True, pytest.fail, True, id='true-or-fail'),
        pytest.param(lazy_or, False, lambda: 'hi!', 'hi!', id='false-or-pass'),
    ],
)
@pytest.mark.parametrize('kind_left', array_kinds)
@pytest.mark.parametrize('kind_right', array_kinds)
def test_lazy_bool(lazy_op, kind_left, left, kind_right, right, expected):
    if kind_left is DaskArray:
        import dask.array as da

        left_unwrapped = left

        def left():
            return da.array(left_unwrapped)

    if kind_right is DaskArray:
        import dask.array as da

        right_unwrapped = right

        def right():
            return da.array(True).map_blocks(
                lambda _: right_unwrapped(), meta=np.bool_(True)
            )

    assert lazy_op(left, right) == expected


@pytest.mark.parametrize(
    ('array_value', 'expected'),
    [
        pytest.param(
            np.random.poisson(size=(100, 100)).astype(np.float64),
            True,
            id='poisson-float64',
        ),
        pytest.param(
            np.random.poisson(size=(100, 100)).astype(np.uint32),
            True,
            id='poisson-uint32',
        ),
        pytest.param(np.random.normal(size=(100, 100)), False, id='normal'),
        pytest.param(np.array([[0, 0, 0], [0, -1, 0], [0, 0, 0]]), False, id='middle'),
    ],
)
def test_check_nonnegative_integers(array_type, array_value, expected):
    X = array_type(array_value)

    received = check_nonnegative_integers(X)
    if isinstance(X, DaskArray):
        assert isinstance(received, DaskArray)
        # compute
        received = received.compute()
    # convert to python bool
    assert received.item() is expected


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


@needs('dask')
@pytest.mark.parametrize(
    ('axis', 'expected'),
    [
        pytest.param(None, False, id='None'),
        pytest.param(0, [True, True, False, False], id='0'),
        pytest.param(1, [False, False, True, True, False, True], id='1'),
    ],
)
@pytest.mark.parametrize('block_type', [np.array, csr_matrix])
def test_is_constant_dask(axis, expected, block_type):
    import dask.array as da

    if (axis is None) and block_type is csr_matrix:
        pytest.skip('Dask has weak support for scipy sparse matrices')

    x_data = [
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
    ]
    x = da.from_array(np.array(x_data), chunks=2).map_blocks(block_type)

    np.testing.assert_array_equal(expected, is_constant(x, axis=axis))

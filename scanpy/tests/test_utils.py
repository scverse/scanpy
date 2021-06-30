from types import ModuleType
from scipy.sparse import csr_matrix, csc_matrix
import numpy as np

from scanpy._utils import descend_classes_and_funcs, check_nonnegative_integers

from anndata.tests.helpers import assert_equal, asarray
import pytest


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


def test_check_nonnegative_integers():

    X = np.random.poisson(size=(100, 100)).astype(np.float64)
    assert check_nonnegative_integers(X) is True
    assert check_nonnegative_integers(-X) is False

    X_ = X + np.random.normal(size=(100, 100))
    assert check_nonnegative_integers(X_) is False

    X = csr_matrix(X)
    assert check_nonnegative_integers(X) is True
    assert check_nonnegative_integers(-X) is False

    X_ = csr_matrix(X_)
    assert check_nonnegative_integers(X_) is False


@pytest.mark.parametrize(
    'array_type', [asarray, csr_matrix, csc_matrix], ids=lambda x: x.__name__
)
def test_is_constant(array_type):
    from scanpy._utils import is_constant

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

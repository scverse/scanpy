from __future__ import annotations

from functools import partial, singledispatch
from numbers import Integral
from typing import Any

import numpy as np
from numpy.typing import NDArray
from numba import njit
from scipy import sparse

from ..._compat import DaskArray


def _check_axis_supported(axis: Any) -> None:
    if not isinstance(axis, Integral):
        raise TypeError("axis must be integer or None.")
    assert axis in (0, 1)


@singledispatch
def is_constant(a, axis: int | None = None) -> bool | NDArray[np.bool_]:
    """
    Check whether values in array are constant.

    Params
    ------
    a
        Array to check
    axis
        Axis to reduce over.


    Returns
    -------
    Boolean array, True values were constant.

    Example
    -------

    >>> a = np.array([[0, 1], [0, 0]])
    >>> a
    array([[0, 1],
           [0, 0]])
    >>> is_constant(a)
    False
    >>> is_constant(a, axis=0)
    array([ True, False])
    >>> is_constant(a, axis=1)
    array([False,  True])
    """
    raise NotImplementedError()


@is_constant.register(np.ndarray)
def _(a: NDArray, axis: int | None = None) -> bool | NDArray[np.bool_]:
    # Should eventually support nd, not now.
    if axis is None:
        return np.array_equal(a, a.flat[0])
    _check_axis_supported(axis)
    if axis == 0:
        return _is_constant_rows(a.T)
    elif axis == 1:
        return _is_constant_rows(a)


def _is_constant_rows(a: NDArray) -> NDArray[np.bool_]:
    b = np.broadcast_to(a[:, 0][:, np.newaxis], a.shape)
    return (a == b).all(axis=1)


@is_constant.register(sparse.csr_matrix)
def _(a: sparse.csr_matrix, axis: int | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        if len(a.data) == np.multiply(*a.shape):
            return is_constant(a.data)
        else:
            return (a.data == 0).all()
    if not isinstance(axis, Integral):
        raise TypeError("axis must be integer or None.")
    assert axis in (0, 1)
    if axis == 1:
        return _is_constant_csr_rows(a.data, a.indices, a.indptr, a.shape)
    elif axis == 0:
        a = a.T.tocsr()
        return _is_constant_csr_rows(a.data, a.indices, a.indptr, a.shape)


@njit
def _is_constant_csr_rows(
    data: NDArray[np.number],
    indices: NDArray[np.integer],
    indptr: NDArray[np.integer],
    shape: tuple[int, int],
):
    N = len(indptr) - 1
    result = np.ones(N, dtype=np.bool_)
    for i in range(N):
        start = indptr[i]
        stop = indptr[i + 1]
        if stop - start == shape[1]:
            val = data[start]
        else:
            val = 0
        for j in range(start, stop):
            if data[j] != val:
                result[i] = False
                break
    return result


@is_constant.register(sparse.csc_matrix)
def _(a: sparse.csc_matrix, axis: int | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        if len(a.data) == np.multiply(*a.shape):
            return is_constant(a.data)
        else:
            return (a.data == 0).all()
    _check_axis_supported(axis)
    if axis == 0:
        return _is_constant_csr_rows(a.data, a.indices, a.indptr, a.shape[::-1])
    elif axis == 1:
        a = a.T.tocsc()
        return _is_constant_csr_rows(a.data, a.indices, a.indptr, a.shape[::-1])


@is_constant.register(DaskArray)
def _(a: DaskArray, axis: int | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        v = a[tuple(0 for _ in range(a.ndim))]
        return (a == v).all()
    _check_axis_supported(axis)
    # TODO: use overlapping blocks and reduction instead of `drop_axis`
    return a.map_blocks(partial(is_constant, axis=axis), drop_axis=axis)

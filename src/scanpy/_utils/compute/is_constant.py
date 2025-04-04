from __future__ import annotations

from functools import partial, singledispatch, wraps
from numbers import Integral
from typing import TYPE_CHECKING, overload

import numba
import numpy as np

from ..._compat import CSCBase, CSRBase, DaskArray, _register_union, njit

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal, TypeVar

    from numpy.typing import NDArray

    from ..._compat import CSBase

    _Array = NDArray | DaskArray | CSBase

    C = TypeVar("C", bound=Callable)


def _check_axis_supported(wrapped: C) -> C:
    @wraps(wrapped)
    def func(a, axis=None):
        if axis is not None:
            if not isinstance(axis, Integral):
                msg = "axis must be integer or None."
                raise TypeError(msg)
            if axis not in (0, 1):
                msg = "We only support axis 0 and 1 at the moment"
                raise NotImplementedError(msg)
        return wrapped(a, axis)

    return func


@overload
def is_constant(a: _Array, axis: None = None) -> bool: ...
@overload
def is_constant(a: _Array, axis: Literal[0, 1]) -> NDArray[np.bool_]: ...


@_check_axis_supported
@singledispatch
def is_constant(
    a: NDArray, axis: Literal[0, 1] | None = None
) -> bool | NDArray[np.bool_]:
    """Check whether values in array are constant.

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
def _(a: NDArray, axis: Literal[0, 1] | None = None) -> bool | NDArray[np.bool_]:
    # Should eventually support nd, not now.
    if axis is None:
        return bool((a == a.flat[0]).all())
    if axis == 0:
        return _is_constant_rows(a.T)
    elif axis == 1:
        return _is_constant_rows(a)


def _is_constant_rows(a: NDArray) -> NDArray[np.bool_]:
    b = np.broadcast_to(a[:, 0][:, np.newaxis], a.shape)
    return (a == b).all(axis=1)


@_register_union(is_constant, CSRBase)
def _(a: CSRBase, axis: Literal[0, 1] | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        if len(a.data) == np.multiply(*a.shape):
            return is_constant(a.data)
        else:
            return (a.data == 0).all()
    if axis == 1:
        return _is_constant_csr_rows(a.data, a.indptr, a.shape)
    elif axis == 0:
        a = a.T.tocsr()
        return _is_constant_csr_rows(a.data, a.indptr, a.shape)


@njit
def _is_constant_csr_rows(
    data: NDArray[np.number],
    indptr: NDArray[np.integer],
    shape: tuple[int, int],
) -> NDArray[np.bool_]:
    n = len(indptr) - 1
    result = np.ones(n, dtype=np.bool_)
    for i in numba.prange(n):
        start = indptr[i]
        stop = indptr[i + 1]
        val = data[start] if stop - start == shape[1] else 0
        for j in range(start, stop):
            if data[j] != val:
                result[i] = False
                break
    return result


@_register_union(is_constant, CSCBase)
def _(a: CSCBase, axis: Literal[0, 1] | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        if len(a.data) == np.multiply(*a.shape):
            return is_constant(a.data)
        else:
            return (a.data == 0).all()
    if axis == 0:
        return _is_constant_csr_rows(a.data, a.indptr, a.shape[::-1])
    elif axis == 1:
        a = a.T.tocsc()
        return _is_constant_csr_rows(a.data, a.indptr, a.shape[::-1])


@is_constant.register(DaskArray)
def _(a: DaskArray, axis: Literal[0, 1] | None = None) -> bool | NDArray[np.bool_]:
    if axis is None:
        v = a[tuple(0 for _ in range(a.ndim))].compute()
        return (a == v).all()
    # TODO: use overlapping blocks and reduction instead of `drop_axis`
    return a.map_blocks(
        partial(is_constant, axis=axis),
        drop_axis=axis,
        meta=np.array([], dtype=a.dtype),
    )

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .._compat import DaskArray

if TYPE_CHECKING:
    from typing import TypeAlias, Union, Callable, Any

    _EagerBool: TypeAlias = Union[bool, np.bool_]
    _BoolScalar: TypeAlias = Union[_EagerBool, DaskArray]
    _Left: TypeAlias = Union[_BoolScalar, Callable[[], _BoolScalar]]
    _Right: TypeAlias = Union[_BoolScalar, Callable[[], _EagerBool]]


def _call_or_return(maybe_cb: Any, validate: bool = False):
    if not callable(maybe_cb):
        return maybe_cb
    rv = maybe_cb()
    if validate and isinstance(rv, DaskArray):
        msg = (
            "Use lazy_*(_, da.array(...)) instead of lazy_*(_, lambda: da.array(...))."
        )
        raise AssertionError(msg)
    return rv


def _lazy_bool_op(left: _Left, right: _Right, /, *, cmp) -> _BoolScalar:
    left = _call_or_return(left)
    if not isinstance(left, DaskArray):
        return cmp(left, right)

    def chain(l: _EagerBool) -> _EagerBool:
        rv = cmp(l, right)
        return rv.compute() if isinstance(rv, DaskArray) else rv

    return left.map_blocks(chain, meta=np.bool_(True))


def lazy_and(left: _Left, right: _Right, /) -> _BoolScalar:
    def cmp(l, r):
        return l and _call_or_return(r, validate=True)

    return _lazy_bool_op(left, right, cmp=cmp)


def lazy_or(left: _Left, right: _Right, /) -> _BoolScalar:
    def cmp(l, r):
        return l or _call_or_return(r, validate=True)

    return _lazy_bool_op(left, right, cmp=cmp)


def get_ufuncs(data: np.ndarray | DaskArray):
    if isinstance(data, DaskArray):
        import dask.array as da

        return da
    return np

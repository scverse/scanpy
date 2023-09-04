from __future__ import annotations
from functools import partial

from typing import TYPE_CHECKING

import numpy as np

from .._compat import DaskArray

if TYPE_CHECKING:
    from typing import TypeAlias, Union, Callable, Any

    _EagerBool: TypeAlias = Union[bool, np.bool_]
    _BoolScalar: TypeAlias = Union[_EagerBool, DaskArray]


def _lazy_bool_op(
    left: _BoolScalar, right: Callable[[], _BoolScalar], /, *, cmp
) -> _BoolScalar:
    if callable(left):
        raise TypeError('only right operand should be callable')
    if not callable(right):
        raise TypeError('right argument must be callable')

    if not isinstance(left, DaskArray):
        return cmp(left, right)

    def chain(l: _EagerBool) -> _EagerBool:
        rv = cmp(l, right)
        return rv.compute() if isinstance(rv, DaskArray) else rv

    return left.map_blocks(chain, meta=np.bool_(True))


lazy_and = partial(_lazy_bool_op, cmp=lambda l, r: l and r())
lazy_or = partial(_lazy_bool_op, cmp=lambda l, r: l or r())


def get_ufuncs(data: np.ndarray | DaskArray):
    if isinstance(data, DaskArray):
        import dask.array as da

        return da
    return np

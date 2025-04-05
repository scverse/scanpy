from __future__ import annotations

from typing import TYPE_CHECKING, overload

import numpy as np

from .._compat import DaskArray

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from .._compat import ZappyArray


@overload
def materialize_as_ndarray(a: ArrayLike) -> np.ndarray: ...


@overload
def materialize_as_ndarray(a: tuple[ArrayLike]) -> tuple[np.ndarray]: ...


@overload
def materialize_as_ndarray(
    a: tuple[ArrayLike, ArrayLike],
) -> tuple[np.ndarray, np.ndarray]: ...


@overload
def materialize_as_ndarray(
    a: tuple[ArrayLike, ArrayLike, ArrayLike],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...


def materialize_as_ndarray(
    a: DaskArray | ArrayLike | tuple[ArrayLike | ZappyArray | DaskArray, ...],
) -> tuple[np.ndarray] | np.ndarray:
    """Compute distributed arrays and convert them to numpy ndarrays."""
    if isinstance(a, DaskArray):
        return a.compute()
    if not isinstance(a, tuple):
        return np.asarray(a)

    if not any(isinstance(arr, DaskArray) for arr in a):
        return tuple(np.asarray(arr) for arr in a)

    import dask.array as da

    return da.compute(*a, sync=True)

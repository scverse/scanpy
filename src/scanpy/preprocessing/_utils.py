from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
from sklearn.random_projection import sample_without_replacement

from .._compat import njit

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray

    from .._compat import CSBase
    from .._utils.random import _LegacyRandom


def sample_comb(
    dims: tuple[int, ...],
    nsamp: int,
    *,
    random_state: _LegacyRandom = None,
    method: Literal[
        "auto", "tracking_selection", "reservoir_sampling", "pool"
    ] = "auto",
) -> NDArray[np.int64]:
    """Randomly sample indices from a grid, without repeating the same tuple."""
    idx = sample_without_replacement(
        np.prod(dims), nsamp, random_state=random_state, method=method
    )
    return np.vstack(np.unravel_index(idx, dims)).T


def _to_dense(X: CSBase, order: Literal["C", "F"] = "C") -> NDArray:
    """Numba kernel for np.toarray() function."""
    out = np.zeros(X.shape, dtype=X.dtype, order=order)
    if X.format == "csr":
        _to_dense_csr_numba(X.indptr, X.indices, X.data, out, X.shape)
    elif X.format == "csc":
        _to_dense_csc_numba(X.indptr, X.indices, X.data, out, X.shape)
    else:
        out = X.toarray(order=order)
    return out


@njit
def _to_dense_csc_numba(
    indptr: NDArray,
    indices: NDArray,
    data: NDArray,
    X: NDArray,
    shape: tuple[int, int],
) -> None:
    for c in numba.prange(X.shape[1]):
        for i in range(indptr[c], indptr[c + 1]):
            X[indices[i], c] = data[i]


@njit
def _to_dense_csr_numba(
    indptr: NDArray,
    indices: NDArray,
    data: NDArray,
    X: NDArray,
    shape: tuple[int, int],
) -> None:
    for r in numba.prange(shape[0]):
        for i in range(indptr[r], indptr[r + 1]):
            X[r, indices[i]] = data[i]

from __future__ import annotations

import random
from collections.abc import Sequence
from contextlib import contextmanager
from typing import TYPE_CHECKING

import numpy as np
from sklearn.utils import check_random_state

from . import ensure_igraph

if TYPE_CHECKING:
    from collections.abc import Generator

    from numpy.typing import NDArray


__all__ = [
    "RNGLike",
    "SeedLike",
    "ith_k_tuple",
    "random_k_tuples",
    "random_strings",
]

SeedLike = int | np.integer | Sequence[int] | np.random.SeedSequence
RNGLike = np.random.Generator | np.random.BitGenerator


# Compatibility with igraph’s RNG


class RNGIgraph:
    """Random number generator for ipgraph so global seed is not changed.

    See :func:`igraph.set_random_number_generator` for the requirements.
    """

    def __init__(self, random_state: int | np.random.RandomState = 0) -> None:
        self._rng = check_random_state(random_state)

    def __getattr__(self, attr: str):
        return getattr(self._rng, "normal" if attr == "gauss" else attr)


@contextmanager
def set_igraph_random_state(
    random_state: int | np.random.RandomState,
) -> Generator[None, None, None]:
    ensure_igraph()
    import igraph

    rng = RNGIgraph(random_state)
    try:
        igraph.set_random_number_generator(rng)
        yield None
    finally:
        igraph.set_random_number_generator(random)


# legacy subsampling TODO: remove


def subsample(
    X: np.ndarray,
    subsample: int = 1,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Subsample a fraction of 1/subsample samples from the rows of X.

    Parameters
    ----------
    X
        Data array.
    subsample
        1/subsample is the fraction of data sampled, n = X.shape[0]/subsample.
    seed
        Seed for sampling.

    Returns
    -------
    Xsampled
        Subsampled X.
    rows
        Indices of rows that are stored in Xsampled.

    """
    if subsample == 1 and seed == 0:
        return X, np.arange(X.shape[0], dtype=int)
    if seed == 0:
        # this sequence is defined simply by skipping rows
        # is faster than sampling
        rows = np.arange(0, X.shape[0], subsample, dtype=int)
        n = rows.size
        Xsampled = np.array(X[rows])
    else:
        if seed < 0:
            msg = f"Invalid seed value < 0: {seed}"
            raise ValueError(msg)
        n = int(X.shape[0] / subsample)
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X, n=n)
    return Xsampled, rows


def subsample_n(
    X: np.ndarray, n: int = 0, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Subsample n samples from rows of array.

    Parameters
    ----------
    X
        Data array.
    n
        Sample size.
    seed
        Seed for sampling.

    Returns
    -------
    Xsampled
        Subsampled X.
    rows
        Indices of rows that are stored in Xsampled.

    """
    if n < 0:
        msg = "n must be greater 0"
        raise ValueError(msg)
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0], size=n, replace=False)
    Xsampled = X[rows]
    return Xsampled, rows


# random k-tuples


def ith_k_tuple(
    indices: NDArray[np.integer], /, *, n: int, k: int
) -> NDArray[np.int64]:
    """Calculate the k-tuple corresponding to the given :func:`itertools.product` index.

    Given the `n**k` possible k-tuples, this function returns the k-tuples that
    `[list(product(range(n), repeat=k))[i] for i in indices]` would evaluate to.

    Parameters
    ----------
    indices
        The tuple indices (must all be in `range(n**k)`).
    n
        The number of possible choices.
    k
        The length of each tuple.

    Returns
    -------
    A 2D array where each row is the k-tuple corresponding to the index.

    Raises
    ------
    ValueError
        If any of the indices are out of range.
    """
    if np.any((indices < 0) | (indices >= n**k)):  # pragma: no cover
        msg = f"Indices are out of range({n**k})."
        raise ValueError(msg)

    power_of_n = n ** np.arange(k - 1, -1, -1)  # [n^(k-1), n^(k-2), ..., n^0]
    return (indices[:, None] // power_of_n) % n


def random_k_tuples(
    size: int, *, n: int, k: int, rng: SeedLike | RNGLike | None = None
) -> NDArray[np.int64]:
    """Draw `size` distinct `k`-tuples of values in `range(n)` from the `n**k` possible tuples.

    Parameters
    ----------
    size
        The number of distinct k-tuples to draw.
    n
        The number of possible choices.
    k
        The length of the tuple.
    rng
        The random number generator to use.

    Returns
    -------
    A 2D array where each row is a distinct k-tuple.

    Raises
    ------
    ValueError
        If `size` is greater than the total number of possible k-tuples.
    """
    rng = np.random.default_rng(rng)

    total_possible = n**k
    if size > total_possible:  # pragma: no cover
        msg = f"Error: Cannot draw {size} unique k-tuples. Total possible is {total_possible}."
        raise ValueError(msg)

    indices = rng.choice(total_possible, size=size, replace=False)
    return ith_k_tuple(indices, n=n, k=k)


def random_strings(
    size: int, *, length: int, alphabet: str, rng: SeedLike | RNGLike | None = None
) -> NDArray[np.str_]:
    """Draw `size` distinct strings of length `k` from the given alphabet.

    Parameters
    ----------
    size
        The number of distinct strings to draw.
    length
        The length of each string.
    alphabet
        The alphabet to draw from.
    rng
        The random number generator to use.

    Returns
    -------
    A 1D array where each element is a distinct string of length `length`
    """
    letters = np.array(list(alphabet), dtype="U1")
    indices = random_k_tuples(size, n=len(letters), k=length, rng=rng)
    return letters[indices].view(f"U{length}").squeeze()

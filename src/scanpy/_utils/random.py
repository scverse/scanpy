from __future__ import annotations

import random
from collections.abc import Sequence
from contextlib import contextmanager
from functools import WRAPPER_ASSIGNMENTS, wraps
from typing import TYPE_CHECKING

import numpy as np

from . import ensure_igraph

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Self

    from numpy.typing import NDArray


__all__ = [
    "RNGLike",
    "SeedLike",
    "_LegacyRandom",
    "_if_legacy_apply_global",
    "accepts_legacy_random_state",
    "ith_k_tuple",
    "legacy_random_state",
    "random_k_tuples",
    "random_str",
]

type SeedLike = int | np.integer | Sequence[int] | np.random.SeedSequence
type RNGLike = np.random.Generator | np.random.BitGenerator
type _LegacyRandom = int | np.random.RandomState | None


###################################
# Compatibility with igraph’s RNG #
###################################


class _RNGIgraph:
    """Random number generator for igraph so global seed is not changed.

    See :func:`igraph.set_random_number_generator` for the requirements.
    """

    def __init__(self, rng: SeedLike | RNGLike | None) -> None:
        self._rng = np.random.default_rng(rng)

    def getrandbits(self, k: int) -> int:
        lims = np.iinfo(np.uint64)
        i = int(self._rng.integers(0, lims.max, dtype=np.uint64))
        return i & ((1 << k) - 1)

    def randint(self, a: int, b: int) -> np.int64:
        return self._rng.integers(a, b + 1)

    def __getattr__(self, attr: str):
        return getattr(self._rng, "normal" if attr == "gauss" else attr)


@contextmanager
def set_igraph_rng(rng: SeedLike | RNGLike | None) -> Generator[None]:
    ensure_igraph()
    import igraph

    ig_rng = _RNGIgraph(rng)
    try:
        igraph.set_random_number_generator(ig_rng)
        yield None
    finally:
        igraph.set_random_number_generator(random)


###################################
# Compatibility with legacy numpy #
###################################


class _FakeRandomGen(np.random.Generator):
    _arg: _LegacyRandom
    _state: np.random.RandomState

    def __init__(
        self, arg: _LegacyRandom, state: np.random.RandomState | None = None
    ) -> None:
        self._arg = arg
        self._state = np.random.RandomState(arg) if state is None else state
        super().__init__(self._state._bit_generator)

    @classmethod
    def wrap_global(
        cls,
        arg: _LegacyRandom = None,
        state: np.random.RandomState | None = None,
    ) -> Self:
        """Create a generator that wraps the global `RandomState` backing the legacy `np.random` functions."""
        if arg is not None:
            if isinstance(arg, np.random.RandomState):
                np.random.set_state(arg.get_state(legacy=False))
                return _FakeRandomGen(arg, state)
            np.random.seed(arg)
        return _FakeRandomGen(arg, np.random.RandomState(np.random.get_bit_generator()))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, _FakeRandomGen):
            return False
        return self._arg == other._arg

    def __hash__(self) -> int:
        return hash((type(self), self._arg))

    @classmethod
    def _delegate(cls) -> None:
        names = dict(integers="randint")
        for name, meth in np.random.Generator.__dict__.items():
            if name.startswith("_") or not callable(meth):
                continue

            def mk_wrapper(name: str, meth):
                # Old pytest versions try to run the doctests
                @wraps(meth, assigned=set(WRAPPER_ASSIGNMENTS) - {"__doc__"})
                def wrapper(self: _FakeRandomGen, *args, **kwargs):
                    return getattr(self._state, name)(*args, **kwargs)

                return wrapper

            setattr(cls, names.get(name, name), mk_wrapper(name, meth))


_FakeRandomGen._delegate()


def _if_legacy_apply_global(rng: np.random.Generator) -> np.random.Generator:
    """Re-apply legacy `random_state` semantics when `rng` is a `_FakeRandomGen`.

    This resets the global legacy RNG from the original `_arg` and returns a
    generator which continues drawing from the same internal state.
    """
    if not isinstance(rng, _FakeRandomGen):
        return rng

    return _FakeRandomGen.wrap_global(rng._arg, rng._state)


def legacy_random_state(
    rng: SeedLike | RNGLike | None, *, always_state: bool = False
) -> _LegacyRandom:
    """Convert a np.random.Generator into a legacy `random_state` argument.

    If `rng` is already a `_FakeRandomGen`, return its original `_arg` attribute.
    """
    if isinstance(rng, _FakeRandomGen):
        return rng._state if always_state else rng._arg
    rng = np.random.default_rng(rng)
    return np.random.RandomState(rng.bit_generator.spawn(1)[0])


def accepts_legacy_random_state[**P, R](
    random_state_default: _LegacyRandom,
) -> callable[[callable[P, R]], callable[P, R]]:
    """Make a function accept `random_state: _LegacyRandom` and pass it as `rng`.

    If the decorated function is called with a `random_state` argument,
    it’ll be wrapped in a :class:`_FakeRandomGen`.
    Passing both ``rng`` and ``random_state`` at the same time is an error.
    If neither is given, ``random_state_default`` is used.
    """

    def decorator(func: callable[P, R]) -> callable[P, R]:
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            match "random_state" in kwargs, "rng" in kwargs:
                case True, True:
                    msg = "Specify at most one of `rng` and `random_state`."
                    raise TypeError(msg)
                case True, False:
                    kwargs["rng"] = _FakeRandomGen(kwargs.pop("random_state"))
                case False, False:
                    kwargs["rng"] = _FakeRandomGen(random_state_default)
            return func(*args, **kwargs)

        return wrapper

    return decorator


###################
# Random k-tuples #
###################


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


def random_str(
    size: int = 1, *, length: int, alphabet: str, rng: SeedLike | RNGLike | None = None
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
    A 0-1D array where each element is a distinct string of length `length`.
    """
    letters = np.array(list(alphabet), dtype="U1")
    indices = random_k_tuples(size, n=len(letters), k=length, rng=rng)
    return letters[indices].view(f"U{length}").squeeze()

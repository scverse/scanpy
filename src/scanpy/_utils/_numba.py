from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING

import numba

from .._compat import warn

if TYPE_CHECKING:
    from collections.abc import Generator


def _resolve_n_threads(n_threads: int) -> int:
    n_max = numba.config.NUMBA_NUM_THREADS
    if n_threads < -1:
        msg = f"n_threads < -1 is not supported for this operation; treating {n_threads} as -1 ({n_max} threads, as per numba.config.NUMBA_NUM_THREADS)."
        warn(msg, UserWarning)
    if n_threads < 0:
        return n_max
    return max(1, min(n_threads, n_max))


@contextmanager
def _numba_thread_limit(n_threads: int | None) -> Generator[int | None]:
    """Temporarily set Numba's thread count and restore it on exit.

    Resolves ``n_threads`` (``-1`` selects ``numba.config.NUMBA_NUM_THREADS``)
    and yields the resolved count. ``None`` leaves Numba unchanged and yields ``None``.
    """
    if n_threads is None:
        yield None
        return

    resolved = _resolve_n_threads(n_threads)
    previous = numba.get_num_threads()
    numba.set_num_threads(resolved)
    try:
        yield resolved
    finally:
        numba.set_num_threads(previous)

from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING

import numba

if TYPE_CHECKING:
    from collections.abc import Generator


@contextmanager
def _numba_thread_limit(n_threads: int | None) -> Generator[None]:
    """Temporarily set Numba's thread count and restore it on exit."""
    if n_threads is None:
        yield
        return

    previous = numba.get_num_threads()
    n_threads = max(1, min(n_threads, numba.config.NUMBA_NUM_THREADS))
    numba.set_num_threads(n_threads)
    try:
        yield
    finally:
        numba.set_num_threads(previous)

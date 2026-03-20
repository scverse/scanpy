from __future__ import annotations

from functools import cache

import numba as nb
import numpy as np
from numba.extending import overload_method
from numba.np.numpy_support import is_nonelike

__all__ = ["add_np_generator_choice"]


@cache  # make sure this is only called once
def add_np_generator_choice():
    """Add `np.random.Generator.choice` to numba if necessary."""

    @nb.njit  # noqa: TID251
    def test(rng: np.random.Generator) -> int:
        return rng.choice(1)

    try:
        test(np.random.default_rng())
    except nb.TypingError:
        _register_np_generator_choice()


def _register_np_generator_choice():

    @overload_method(nb.types.NumPyRandomGeneratorType, "choice")
    def NumPyRandomGeneratorType_choice(inst, a, size=None, replace=True):  # noqa: FBT002, N802
        if is_nonelike(size):

            def impl(inst, a, size=None, replace=True):  # noqa: FBT002
                return inst.integers(0, a)

        else:

            def impl(inst, a, size=None, replace=True):  # noqa: FBT002
                if replace:
                    return inst.integers(0, a, size=size)
                else:
                    pool = np.arange(a)
                    for i in range(size):
                        j = inst.integers(i, a)
                        pool[i], pool[j] = pool[j], pool[i]
                    return pool[:size].copy()

        return impl

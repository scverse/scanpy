from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .._utils.random import _FakeRandomGen

if TYPE_CHECKING:
    from numpy.typing import NDArray


def sample_comb(
    dims: tuple[int, ...], nsamp: int, *, rng: np.random.Generator
) -> NDArray[np.int64]:
    """Randomly sample indices from a grid, without repeating the same tuple."""
    if isinstance(rng, _FakeRandomGen):
        from sklearn.random_projection import sample_without_replacement

        idx = sample_without_replacement(np.prod(dims), nsamp, random_state=rng._arg)
    else:
        idx = rng.choice(np.prod(dims), size=nsamp, replace=False)
    return np.vstack(np.unravel_index(idx, dims)).T

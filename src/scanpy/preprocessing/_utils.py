from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sklearn.random_projection import sample_without_replacement

from .._utils.random import _legacy_random_state

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray


def sample_comb(
    dims: tuple[int, ...],
    nsamp: int,
    *,
    rng: np.random.Generator,
    method: Literal[
        "auto", "tracking_selection", "reservoir_sampling", "pool"
    ] = "auto",
) -> NDArray[np.int64]:
    """Randomly sample indices from a grid, without repeating the same tuple."""
    random_state = _legacy_random_state(rng)
    idx = sample_without_replacement(
        np.prod(dims), nsamp, random_state=random_state, method=method
    )
    return np.vstack(np.unravel_index(idx, dims)).T

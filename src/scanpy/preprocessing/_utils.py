from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sklearn.random_projection import sample_without_replacement

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray

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

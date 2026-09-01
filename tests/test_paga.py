from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from typing import Literal


pytestmark = [needs.igraph]


@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_paga_positions_reproducible(
    subtests: pytest.Subtests, rng_arg: Literal["rng", "random_state"]
) -> None:
    """Check exact reproducibility and effect of random_state on paga positions."""
    # https://github.com/scverse/scanpy/issues/1859
    pbmc = pbmc68k_reduced()
    sc.tl.paga(pbmc, "bulk_labels")

    a = pbmc.copy()
    b = pbmc.copy()
    c = pbmc.copy()

    sc.pl.paga(a, show=False, **{rng_arg: 42})
    sc.pl.paga(b, show=False, **{rng_arg: 42})
    sc.pl.paga(c, show=False, **{rng_arg: 13})

    with subtests.test("reproducible"):
        np.testing.assert_array_equal(a.uns["paga"]["pos"], b.uns["paga"]["pos"])
    with subtests.test("different positions"):
        assert a.uns["paga"]["pos"].tolist() != c.uns["paga"]["pos"].tolist()

from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc
from scanpy._utils.random import _LegacyRng


@pytest.mark.filterwarnings("ignore:.*Observation names are not unique:UserWarning")
def test_sim_toggleswitch() -> None:
    # the reference data was generated using the legacy global RNG seeded with 0
    adata_sim = sc.tl.sim("toggleswitch", rng=_LegacyRng(0))
    adata_ds = sc.datasets.toggleswitch()
    np.testing.assert_allclose(adata_sim.X, adata_ds.X, atol=np.finfo(np.float32).eps)

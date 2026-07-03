from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc


@pytest.mark.filterwarnings("ignore:.*Observation names are not unique:UserWarning")
def test_sim_toggleswitch() -> None:
    adata_sim = sc.tl.sim("toggleswitch")
    adata_ds = sc.datasets.toggleswitch()
    np.allclose(adata_sim.X, adata_ds.X, np.finfo(np.float32).eps)

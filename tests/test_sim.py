from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc


def test_sim_toggleswitch():
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata_sim = sc.tl.sim("toggleswitch")
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata_ds = sc.datasets.toggleswitch()
    np.allclose(adata_sim.X, adata_ds.X, np.finfo(np.float32).eps)

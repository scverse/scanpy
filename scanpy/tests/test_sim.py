from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc


def test_sim_toggleswitch(tmp_write_dir):
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata = sc.tl.sim("toggleswitch")
        np.allclose(adata.X, sc.datasets.toggleswitch().X, np.finfo(np.float32).eps)

from __future__ import annotations

import numpy as np

import scanpy as sc


def test_sim_toggleswitch(tmp_write_dir):
    adata = sc.tl.sim("toggleswitch")
    np.allclose(adata.X, sc.datasets.toggleswitch().X, np.finfo(np.float32).eps)

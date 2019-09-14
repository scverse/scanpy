import scanpy as sc
import numpy as np


def test_sim_toggleswitch():
    adata = sc.tl.sim('toggleswitch')
    np.allclose(adata.X, sc.datasets.toggleswitch().X, np.finfo(np.float32).eps)

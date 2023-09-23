from _pytest.config import directory_arg
from matplotlib.pyplot import get_current_fig_manager
import pandas as pd
import numpy as np
from scipy import sparse

import scanpy as sc

from scanpy.tests.test_embedding_plots import scanpy_adata


def _normalize_codes(codes: np.ndarray) -> np.ndarray:
    """
    Recode [2, 1, 0, 0] into [0, 1, 2, 2]. Used to normalize categorical codes for comparison.
    """
    new_order = np.argsort(pd.unique(codes))
    return new_order[codes]


def test_visium_connectivity():
    adata = scanpy_adata()
    true_groups = _normalize_codes(adata.obs["label"].cat.codes)
    keys = {}

    for i in range(1, 5):
        key = f"spatial_connectivity_{i}_steps"
        sc.pp.visium_connectivity(adata, n_steps=i, key_added=key)
        keys[i] = key

    # Number of rings should not change connected components
    for k in keys.values():
        _, components = sparse.csgraph.connected_components(
            adata.obsp[k], directed=False, return_labels=True
        )
        components = _normalize_codes(components)
        assert np.array_equal(true_groups, components)

    # Each step should only add new edges
    for i in range(1, len(keys) - 1):
        g_cur = adata.obsp[keys[i]]
        g_next = adata.obsp[keys[i + 1]]
        assert (g_cur > g_next).sum() == 0
        assert g_cur.nnz <= g_next.nnz

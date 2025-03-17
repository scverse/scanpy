from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import scanpy.external as sce


def test_cell_demultiplexing():
    import random

    from scipy import stats

    random.seed(52)
    signal = stats.poisson.rvs(1000, 1, 990)
    doublet_signal = stats.poisson.rvs(1000, 1, 10)
    x = np.reshape(stats.poisson.rvs(500, 1, 10000), (1000, 10))
    for idx, signal_count in enumerate(signal):
        col_pos = idx % 10
        x[idx, col_pos] = signal_count

    for idx, signal_count in enumerate(doublet_signal):
        col_pos = (idx % 10) - 1
        x[idx, col_pos] = signal_count

    test_data = AnnData(np.random.randint(0, 100, size=x.shape), obs=pd.DataFrame(x))
    sce.pp.hashsolo(test_data, test_data.obs.columns)

    doublets = ["Doublet"] * 10
    classes = np.repeat(np.arange(10), 98).reshape(98, 10, order="F").ravel().tolist()
    negatives = ["Negative"] * 10
    expected = pd.array(doublets + classes + negatives, dtype="string")
    classification = test_data.obs["Classification"].array.astype("string")
    # This is a bit flaky, so allow some mismatches:
    if (expected != classification).sum() > 3:
        # Compare lists for better error message
        assert classification.tolist() == expected.tolist()

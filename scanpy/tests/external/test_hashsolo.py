from anndata import AnnData
import numpy as np
import scanpy.external as sce


def test_cell_demultiplexing():
    from scipy import stats
    import random

    random.seed(52)
    signal = stats.poisson.rvs(1000, 1, 990)
    doublet_signal = stats.poisson.rvs(1000, 1, 10)
    x = np.reshape(stats.poisson.rvs(5, 1, 10000), (1000, 10))
    for idx, signal_count in enumerate(signal):
        col_pos = idx % 10
        x[idx, col_pos] = signal_count

    for idx, signal_count in enumerate(doublet_signal):
        col_pos = (idx % 10) - 1
        x[idx, col_pos] = signal_count

    test_data = AnnData(np.random.randint(0, 100, size=x.shape), obs=x)
    sce.pp.hashsolo(test_data, test_data.obs.columns)

    doublets = ["Doublet"] * 10
    classes = list(
        np.repeat(np.arange(10), 98).reshape(98, 10, order="F").ravel().astype(str)
    )
    negatives = ["Negative"] * 10
    classification = doublets + classes + negatives
    assert all(test_data.obs["Classification"].astype(str) == classification)

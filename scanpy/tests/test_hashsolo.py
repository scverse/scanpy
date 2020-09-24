import scanpy as sc
from anndata import AnnData
import numpy as np


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
    sc.pp.hashsolo(test_data, test_data.obs.columns)

    doublets = ['Doublet'] * 10
    classes = list(np.repeat(np.arange(10), 98).reshape(98, 10,
                                                        order='F').ravel())
    negatives = ['Negative'] * 10
    classification = doublets + classes + negatives

    assert all(test_data.obs['Classification'] == classification)

    doublets = [2] * 10
    classes = [1] * 980
    negatives = [0] * 10
    classification = doublets + classes + negatives
    ll_results = np.argmax(sc.pp._calculate_log_likelihoods(x, 8)[0],
                           axis=1)
    assert all(ll_results == classification)

    bayes_results = sc.pp._calculate_bayes_rule(x, [.1, .8, .1], 8)
    assert all(bayes_results['most_likely_hypothesis'] == classification)

    singlet_prior = .99999999999999999
    other_prior = (1 - singlet_prior)/2
    bayes_results = sc.pp._calculate_bayes_rule(x,
                                                [other_prior,
                                                 singlet_prior,
                                                 other_prior], 8)
    assert all(bayes_results['most_likely_hypothesis'] == 1)
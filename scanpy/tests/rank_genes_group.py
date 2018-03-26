from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from scipy import sparse as sp
from numpy.random import negative_binomial, binomial, seed

from anndata import AnnData
from scanpy.api.tl import rank_genes_groups


HERE = Path(__file__).parent


# We test results for a simple generic example
# Tests are conducted for sparse and non-sparse AnnData objects.
# Due to minor changes in multiplication implementation for sparse and non-sparse objects,
# results differ (very) slightly


def get_example_data(*, sparse=False):
    # create test object
    adata = AnnData(np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5)))

    # The following construction is inefficient, but makes sure that the same data is used in the sparse case
    if sparse:
        adata.X = sp.csr_matrix(adata.X)

    # Create cluster according to groups
    adata.obs['true_groups'] = pd.Categorical(np.concatenate((
        np.zeros((10,), dtype=int),
        np.ones( (90,), dtype=int),
    )))

    return adata


def get_true_scores():
    # Here, we have saved the true results
    # Note: Default value is on copying = true.
    with Path(HERE, 'objs_t_test.pkl').open('rb') as f:
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with Path(HERE, 'objs_wilcoxon.pkl').open('rb') as f:
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)

    return true_names_t_test,  true_names_wilcoxon, \
           true_scores_t_test, true_scores_wilcoxon


def test_results_dense():
    seed(1234)

    adata = get_example_data()

    true_names_t_test,  true_names_wilcoxon, \
    true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    # Now run the rank_genes_groups, test functioning.
    rank_genes_groups(adata, 'true_groups', n_genes=20, method='t-test')
    print(true_scores_t_test)
    print(adata.uns['rank_genes_groups']['scores'])
    print(true_names_t_test)
    print(adata.uns['rank_genes_groups']['names'])
    assert np.array_equal(true_scores_t_test, adata.uns['rank_genes_groups']['scores'])
    assert np.array_equal(true_names_t_test, adata.uns['rank_genes_groups']['names'])

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='wilcoxon')
    assert np.array_equal(true_scores_wilcoxon, adata.uns['rank_genes_groups']['scores'])
    assert np.array_equal(true_names_wilcoxon, adata.uns['rank_genes_groups']['names'])


def test_results_sparse():
    seed(1234)

    adata_sparse = get_example_data(sparse=True)

    true_names_t_test,  true_names_wilcoxon, \
    true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    rank_genes_groups(adata_sparse, 'true_groups', n_genes=20, method='t-test')

    # Here, we allow a minor error tolerance due to different multiplication for sparse/non-spars objects
    ERROR_TOLERANCE = 5e-7
    max_error = 0
    for i, k in enumerate(adata_sparse.uns['rank_genes_groups']['scores']):
        max_error = max(max_error, abs(
            adata_sparse.uns['rank_genes_groups']['scores'][i][0] - true_scores_t_test[i][0]))
        max_error = max(max_error, abs(
            adata_sparse.uns['rank_genes_groups']['scores'][i][1] - true_scores_t_test[i][1]))

    # assert np.array_equal(true_scores_t_test,adata_sparse.uns['rank_genes_groups']['scores'])
    assert max_error < ERROR_TOLERANCE
    rank_genes_groups(adata_sparse, 'true_groups', n_genes=20, method='wilcoxon')
    assert np.array_equal(true_scores_wilcoxon, adata_sparse.uns['rank_genes_groups']['scores'])
    assert np.array_equal(true_names_wilcoxon, adata_sparse.uns['rank_genes_groups']['names'])

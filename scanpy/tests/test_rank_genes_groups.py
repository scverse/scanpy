from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from scipy import sparse as sp
from numpy.random import negative_binomial, binomial, seed

from anndata import AnnData
from scanpy.tools import rank_genes_groups


HERE = Path(__file__).parent / Path('_data/')


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
        np.ones((90,), dtype=int),
    )))

    return adata


def get_true_scores():
    with Path(HERE, 'objs_t_test.pkl').open('rb') as f:
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with Path(HERE, 'objs_wilcoxon.pkl').open('rb') as f:
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)

    return true_names_t_test, true_names_wilcoxon,\
           true_scores_t_test, true_scores_wilcoxon


def test_results_dense():
    seed(1234)

    adata = get_example_data()

    true_names_t_test, true_names_wilcoxon,\
    true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='t-test')
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_t_test[name], adata.uns['rank_genes_groups']['scores'][name])
    assert np.array_equal(true_names_t_test, adata.uns['rank_genes_groups']['names'])

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='wilcoxon')
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_wilcoxon[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])
    assert np.array_equal(true_names_wilcoxon[:7], adata.uns['rank_genes_groups']['names'][:7])


def test_results_sparse():
    seed(1234)

    adata = get_example_data(sparse=True)

    true_names_t_test, true_names_wilcoxon,\
    true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='t-test')

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='t-test')
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_t_test[name], adata.uns['rank_genes_groups']['scores'][name])
    assert np.array_equal(true_names_t_test, adata.uns['rank_genes_groups']['names'])    

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='wilcoxon')
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_wilcoxon[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])
    assert np.array_equal(true_names_wilcoxon[:7], adata.uns['rank_genes_groups']['names'][:7])


def test_results_layers():
    seed(1234)

    adata = get_example_data(sparse=False)
    adata.layers["to_test"] = adata.X.copy()
    adata.X = adata.X * np.random.randint(0, 2, adata.shape, dtype=bool)

    true_names_t_test, true_names_wilcoxon,\
    true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    # Wilcoxon
    rank_genes_groups(
        adata, 'true_groups', method='wilcoxon',
        layer="to_test", use_raw=False, n_genes=20,
    )
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_wilcoxon[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])

    rank_genes_groups(adata, 'true_groups', method='wilcoxon', n_genes=20)
    for name in true_scores_t_test.dtype.names:
        assert not np.allclose(true_scores_wilcoxon[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])

    # t-test
    rank_genes_groups(
        adata, 'true_groups', method='t-test',
        layer="to_test", use_raw=False, n_genes=20,
    )
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(true_scores_t_test[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])

    rank_genes_groups(adata, 'true_groups', method='t-test', n_genes=20)
    for name in true_scores_t_test.dtype.names:
        assert not np.allclose(true_scores_t_test[name][:7], adata.uns['rank_genes_groups']['scores'][name][:7])

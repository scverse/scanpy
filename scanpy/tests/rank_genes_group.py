import numpy as np
from numpy import ma
from scipy import sparse as sp
# we don’t need this in requirements.txt, as it’s only needed for testing
from pytest import mark
from numpy.random import negative_binomial, binomial, seed
from scanpy.api.tl import rank_genes_groups
import pickle

from scanpy.data_structs.ann_data import AnnData, BoundStructArray, SMP_INDEX

# We test results for a simple generic example
# Tests are conducted for sparse and non-sparse AnnData objects.
# Due to minor changes in multiplication implementation for sparse and non-sparse objects, results differ (very) slightly


def test_results_dense():
    # set seed
    seed(1234)
    # create test object
    adata = AnnData(np.multiply(binomial(1, 0.15, (100, 20)),negative_binomial(2, 0.25, (100, 20))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(binomial(1, 0.9, (10, 5)),negative_binomial(1, 0.5, (10, 5)))

    # Create cluster according to groups

    smp = 'true_groups'
    true_groups = np.zeros((2, 100), dtype=bool)
    true_groups[0, 0:10] = 1
    true_groups[1, 10:100] = 1
    adata.add[smp + '_masks'] = true_groups
    adata.add[smp + '_order'] = np.asarray(['0', '1'])
    # Now run the rank_genes_groups, test functioning.
    # Note: Default value is on copying = true.
    with open('objs_t_test.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with open('objs_wilcoxon.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)
    rank_genes_groups(adata, 'true_groups', n_genes=20, test_type='t_test')
    assert np.array_equal(true_scores_t_test,adata.add['rank_genes_groups_gene_scores'])
    assert np.array_equal(true_names_t_test, adata.add['rank_genes_groups_gene_names'])

    rank_genes_groups(adata, 'true_groups', n_genes=20, test_type='wilcoxon')
    assert np.array_equal(true_scores_wilcoxon,adata.add['rank_genes_groups_gene_scores'])
    assert np.array_equal(true_names_wilcoxon,adata.add['rank_genes_groups_gene_names'])


def test_results_sparse():
    # set seed
    seed(1234)
    # The following construction is inefficient, but makes sure that the same data is used in the sparse case
    adata = AnnData(np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5)))

    adata_sparse = AnnData(sp.csr_matrix(adata.X))

    # Create cluster according to groups

    smp = 'true_groups'
    true_groups = np.zeros((2, 100), dtype=bool)
    true_groups[0, 0:10] = 1
    true_groups[1, 10:100] = 1
    adata_sparse.add[smp + '_masks'] = true_groups
    adata_sparse.add[smp + '_order'] = np.asarray(['0', '1'])

    # Here, we have saved the true results

    # Now run the rank_genes_groups, test functioning.
    # Note: Default value is on copying = true.
    with open('objs_t_test.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with open('objs_wilcoxon.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)
    rank_genes_groups(adata_sparse, 'true_groups', n_genes=20, test_type='t_test')
    # Here, we allow a minor error tolerance due to different multiplication for sparse/non-spars objects
    ERROR_TOLERANCE=5e-7
    max_error=0
    for i, k in enumerate(adata_sparse.add['rank_genes_groups_gene_scores']):
        max_error = max(max_error, abs(
            adata_sparse.add['rank_genes_groups_gene_scores'][i][0] - true_scores_t_test[i][0]))
        max_error = max(max_error, abs(
            adata_sparse.add['rank_genes_groups_gene_scores'][i][1] - true_scores_t_test[i][1]))
    # assert np.array_equal(true_scores_t_test,adata_sparse.add['rank_genes_groups_gene_scores'])
    assert max_error<ERROR_TOLERANCE
    rank_genes_groups(adata_sparse, 'true_groups', n_genes=20, test_type='wilcoxon')
    assert np.array_equal(true_scores_wilcoxon,adata_sparse.add['rank_genes_groups_gene_scores'])
    assert np.array_equal(true_names_wilcoxon, adata_sparse.add['rank_genes_groups_gene_names'])

def test_compute_distribution():
    # set seed
    seed(1234)
    # create test object
    adata = AnnData(np.multiply(binomial(1, 0.15, (100, 20)),negative_binomial(2, 0.25, (100, 20))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(binomial(1, 0.9, (10, 5)),negative_binomial(1, 0.5, (10, 5)))

    # Create cluster according to groups

    smp = 'true_groups'
    true_groups = np.zeros((2, 100), dtype=bool)
    true_groups[0, 0:10] = 1
    true_groups[1, 10:100] = 1
    adata.add[smp + '_masks'] = true_groups
    adata.add[smp + '_order'] = np.asarray(['0', '1'])
    # Now run the rank_genes_groups, test functioning.
    # Note: Default value is on copying = true.
    with open('objs_t_test.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with open('objs_wilcoxon.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)
    rank_genes_groups(adata, 'true_groups', n_genes=20,compute_distribution=True, test_type='t_test')
    assert np.array_equal(true_scores_t_test,adata.add['rank_genes_groups_gene_scores'])
    assert np.array_equal(true_names_t_test, adata.add['rank_genes_groups_gene_names'])

    rank_genes_groups(adata, 'true_groups', n_genes=20,compute_distribution=True, test_type='wilcoxon')
    assert np.array_equal(true_scores_wilcoxon,adata.add['rank_genes_groups_gene_scores'])
    assert np.array_equal(true_names_wilcoxon,adata.add['rank_genes_groups_gene_names'])

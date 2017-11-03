import numpy as np
from numpy import ma
from scipy import sparse as sp
# we don’t need this in requirements.txt, as it’s only needed for testing
from pytest import mark
from numpy.random import negative_binomial, binomial, seed
from scanpy.api.tl import rank_genes_groups

from scanpy.data_structs.ann_data import AnnData, BoundStructArray, SMP_INDEX


# TODO: Useful tests:
# 1. Matching dimenstions: Output agrees with wanted input (especially: rank_gene_groups and rank_gene_groups_order DONE
# 2. Test: Copying: Logging works as required
# 3. Test: Make sure adata is not modified (i.e. when not copying, result after test is equal to the one before     DONE
# 4. Test: Make sure test works for both sparse and non-sparse input (and yield same results)
# 5. Tets: Make


def test_dimensions_dense():
    # The following initialization will be (probably) the same for all  tests
    # set seed
    seed(1234)
    # create test object
    adata = AnnData(np.multiply(binomial(1, 0.15, (100, 20)),negative_binomial(2, 0.25, (100, 20))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(binomial(1, 0.9, (10, 5)),negative_binomial(1, 0.5, (10, 5)))
    # # Create cluster according to groups
    smp = 'true_groups'
    true_groups = np.zeros((2, 100), dtype=bool)
    true_groups[0, 0:10] = 1
    true_groups[1, 10:100] = 1
    adata.add[smp + '_masks'] = true_groups
    adata.add[smp + '_order'] = np.asarray(['0', '1'])


    # Now run the rank_genes_groups, test functioning.
    implemented_methods = ['t_test', 'wilcoxon', 'wilcoxon2', 'wilcoxon3']
    for i, k in enumerate(implemented_methods):
        rank_genes_groups(adata, 'true_groups', n_genes=10, test_type=k)
        assert adata.add['rank_genes_groups_gene_scores'].shape == (10,)
        adata.add['rank_genes_groups_gene_names'].shape == (10,)


    pass


def test_copying_dense():
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
    adata_save=adata
    # Note: Default value is on copying = true.
    # TODO: Loop over methods
    implemented_methods=['t_test', 'wilcoxon', 'wilcoxon2', 'wilcoxon3' ]
    for i, k in enumerate(implemented_methods):
        rank_genes_groups(adata, 'true_groups', n_genes=20, test_type=k)
        assert adata_save == adata
        adata=adata_save


def test_copying_sparse():
    # set seed
    seed(1234)
    # create test object
    adata = AnnData(sp.csr_matrix(np.multiply(binomial(1, 0.15, (100, 20)),negative_binomial(2, 0.25, (100, 20)))))
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = sp.csr_matrix(np.multiply(binomial(1, 0.9, (10, 5)),negative_binomial(1, 0.5, (10, 5))))

    # Create cluster according to groups

    smp = 'true_groups'
    true_groups = np.zeros((2, 100), dtype=bool)
    true_groups[0, 0:10] = 1
    true_groups[1, 10:100] = 1
    adata.add[smp + '_masks'] = true_groups
    adata.add[smp + '_order'] = np.asarray(['0', '1'])

    # Now run the rank_genes_groups, test functioning.
    adata_save=adata
    # Note: Default value is on copying = true.
    # TODO: Loop over methods
    implemented_methods=['t_test', 'wilcoxon', 'wilcoxon2', 'wilcoxon3' ]
    for i, k in enumerate(implemented_methods):
        rank_genes_groups(adata, 'true_groups', n_genes=20, test_type=k)
        assert adata_save == adata
        adata=adata_save

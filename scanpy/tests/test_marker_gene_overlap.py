import scanpy as sc
import numpy as np

def generate_test_data():
    # Create an artificial data set
    test_data = sc.AnnData(X = np.ones((9,10)))
    test_data.uns['rank_genes_groups'] = dict()
    test_data.uns['rank_genes_groups']['names'] = np.rec.fromarrays(
        [['a', 'b','c','d','e'], ['a','f','g','h','i']])
    test_data.uns['rank_genes_groups']['pvals_adj'] = np.rec.fromarrays(
        [[0.001, 0.01, 0.02, 0.05, 0.6], [0.001, 0.01, 0.02, 0.05, 0.6]])

    marker_genes = {'type 1':{'a','b','c'}, 'type 2':{'a','f','g'}}

    return test_data,marker_genes
    

def test_marker_overlap_base():
    # Test all overlap calculations on artificial data
    test_data,marker_genes = generate_test_data()

    t1 = sc.tl.marker_gene_overlap(test_data, marker_genes)
    
    print(t1)

    assert t1.iloc[1,1] == 3.0
    assert t1.iloc[0,0] == 3.0


def test_marker_overlap_normalization():
    test_data,marker_genes = generate_test_data()

    t2 = sc.tl.marker_gene_overlap(test_data, marker_genes, normalize='reference')
    t3 = sc.tl.marker_gene_overlap(test_data, marker_genes, normalize='data')
    print(t2)
    print(t3)

    assert t2.iloc[0,0] == 1.0
    assert t3.iloc[1,1] == 0.6


def test_marker_overlap_methods():
    test_data,marker_genes = generate_test_data()

    t4 = sc.tl.marker_gene_overlap(test_data, marker_genes, method='overlap_coef')
    t5 = sc.tl.marker_gene_overlap(test_data, marker_genes, method='jaccard')
    print(t4)
    print(t5)
    assert t4.iloc[0,0] == 1.0
    assert t5.iloc[0,0] == 0.6


def test_marker_overlap_subsetting():
    test_data,marker_genes = generate_test_data()

    t6 = sc.tl.marker_gene_overlap(test_data, marker_genes, top_n_markers=2)
    t7 = sc.tl.marker_gene_overlap(test_data, marker_genes, adj_pval_threshold=0.01)
    print(t6)
    print(t7)

    assert t6.iloc[0,0] == 2
    assert t7.iloc[0,0] == 1.0


import scanpy as sc
import numpy as np

def test_marker_overlap():
    # Test all overlap calculations on artificial data
    test_data = sc.AnnData(X = np.ones((9,10)))
    test_data.uns['rank_genes_groups'] = dict()
    test_data.uns['rank_genes_groups']['names'] = np.rec.fromarrays(
        [['a', 'b','c','d','e'], ['a','f','g','h','i']])
    test_data.uns['rank_genes_groups']['pvals_adj'] = np.rec.fromarrays(
        [[0.001, 0.01, 0.02, 0.05, 0.6], [0.001, 0.01, 0.02, 0.05, 0.6]])

    marker_genes = {'type 1':{'a','b','c'}, 'type 2':{'a','f','g'}}

    sc.tl.marker_gene_overlap(test_data, marker_genes, key_added='test1')
    sc.tl.marker_gene_overlap(test_data, marker_genes, normalize='reference', key_added='test2')
    sc.tl.marker_gene_overlap(test_data, marker_genes, method='overlap_coef', key_added='test3')
    sc.tl.marker_gene_overlap(test_data, marker_genes, method='jaccard', key_added='test4')
    sc.tl.marker_gene_overlap(test_data, marker_genes, top_n_markers=2, key_added='test5')
    sc.tl.marker_gene_overlap(test_data, marker_genes, adj_pval_threshold=0.01, key_added='test6')

    assert test_data.uns['test1'].iloc[1,1] == 3.0
    assert test_data.uns['test1'].iloc[0,0] == 3.0
    assert test_data.uns['test2'].iloc[0,0] == 0.75
    assert test_data.uns['test3'].iloc[0,0] == 1.0
    assert test_data.uns['test4'].iloc[0,0] == 0.6
    assert test_data.uns['test5'].iloc[0,0] == 2
    assert test_data.uns['test6'].iloc[0,0] == 1.0

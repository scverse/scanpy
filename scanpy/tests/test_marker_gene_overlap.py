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

    t1 = sc.tl.marker_gene_overlap(test_data, marker_genes)
    t2 = sc.tl.marker_gene_overlap(test_data, marker_genes, normalize='reference')
    t3 = sc.tl.marker_gene_overlap(test_data, marker_genes, method='overlap_coef')
    t4 = sc.tl.marker_gene_overlap(test_data, marker_genes, method='jaccard')
    t5 = sc.tl.marker_gene_overlap(test_data, marker_genes, top_n_markers=2)
    t6 = sc.tl.marker_gene_overlap(test_data, marker_genes, adj_pval_threshold=0.01)

    assert t1.iloc[1,1] == 3.0
    assert t1.iloc[0,0] == 3.0
    assert t2.iloc[0,0] == 0.75
    assert t3.iloc[0,0] == 1.0
    assert t4.iloc[0,0] == 0.6
    assert t5.iloc[0,0] == 2
    assert t6.iloc[0,0] == 1.0

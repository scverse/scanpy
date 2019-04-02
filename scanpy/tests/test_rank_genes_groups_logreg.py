import numpy as np
import scanpy as sc


def test_rank_genes_groups_with_renamed_categories():
    adata = sc.datasets.blobs(n_variables=4, n_centers=3, n_observations=200)
    assert np.allclose(adata.X[1], [9.214668, -2.6487126, 4.2020774, 0.51076424])
    
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)

    for method in ['logreg', 't-test']:
        sc.tl.louvain(adata)

        sc.tl.rank_genes_groups(adata, 'louvain', method=method)
        assert adata.uns['rank_genes_groups']['names'].dtype.names == ('0', '1', '2')
        assert adata.uns['rank_genes_groups']['names'][0].tolist() == ('3', '1', '0')

        adata.rename_categories('louvain', ['Zero', 'One', 'Two'])
        assert adata.uns['rank_genes_groups']['names'][0].tolist() == ('3', '1', '0')

        sc.tl.rank_genes_groups(adata, 'louvain', method=method)
        assert adata.uns['rank_genes_groups']['names'][0].tolist() == ('3', '1', '0')

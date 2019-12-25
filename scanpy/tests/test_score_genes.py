import numpy as np
import scanpy as sc
from anndata import AnnData


def test_add_score():
    # TODO: write a test that costs less resources and is more meaningful
    adata = AnnData(np.random.randint(0, 1000, 100000).reshape((100, 1000)))
    gene_names = np.array(
        [''.join(map(chr, np.random.randint(65, 90, 6))) for _ in range(2000)]
    )
    adata.var_names = gene_names[:1000]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    some_genes = np.r_[
        np.unique(gene_names[np.random.randint(0, 1000, 10)]),
        np.unique(gene_names[np.random.randint(1000, 2000, 3)]),
    ]
    sc.tl.score_genes(adata, some_genes, score_name='Test')
    assert adata.obs['Test'].dtype == 'float32'

    adata.X *= 2.0
    adata.raw = adata
    sc.tl.score_genes(adata, some_genes, score_name='Test_raw')
    assert adata.obs['Test_raw'].dtype == 'float32'
    assert adata.obs['Test_raw'].mean() != adata.obs['Test'].mean()

    adata.layers['layer'] = adata.raw.X
    sc.tl.score_genes(adata, some_genes, score_name='Test_layer', layer='layer')
    assert adata.obs['Test_layer'].dtype == 'float32'
    assert np.allclose(adata.obs['Test_raw'], adata.obs['Test_layer'])

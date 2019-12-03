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
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    some_genes = np.r_[
        np.unique(gene_names[np.random.randint(0, 1000, 10)]),
        np.unique(gene_names[np.random.randint(1000, 2000, 3)]),
    ]
    sc.tl.score_genes(adata, some_genes, score_name='Test')
    assert adata.obs['Test'].dtype == 'float32'

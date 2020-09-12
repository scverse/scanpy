import scanpy as sc
import pytest


@pytest.mark.parametrize('metric', ['cosine', 'euclidean'])
def test_neighbors_metric(metric):
    adata = sc.datasets.pbmc68k_reduced()
    sc.pp.neighbors(adata, random_state=0, use_rep="X", metric=metric)

    assert 'neighbors' in adata.uns
    assert adata.obsp["connectivities"].shape == (700, 700)
    assert adata.obsp["distances"].shape == (700, 700)

import scanpy as sc
import pytest
from scipy.sparse import csr_matrix


@pytest.mark.parametrize('metric', ['cosine', 'euclidean', 'hellinger', 'll_dirichlet'])
@pytest.mark.parametrize('adata_type', ['dense', 'sparse'])
def test_neighbors_metric(metric, adata_type):
    adata = sc.datasets.pbmc68k_reduced()
    if adata_type == 'sparse':
        adata.X = csr_matrix(adata.X)

    sc.pp.neighbors(adata, random_state=0, use_rep="X", metric=metric)

    assert 'neighbors' in adata.uns
    assert adata.obsp["connectivities"].shape == (700, 700)
    assert adata.obsp["distances"].shape == (700, 700)

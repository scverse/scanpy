import numpy as np
import scanpy as sc
from scanpy.preprocessing._simple import N_PCS

from sklearn.neighbors import KDTree

@pytest.fixture
def adatas():
    pbmc = sc.datasets.pbmc68k_reduced()
    n_split = 500
    adata_ref = sc.AnnData(pbmc.X[:n_split, :], obs=pbmc.obs.iloc[:n_split])
    adata_new = sc.AnnData(pbmc.X[n_split:, :])

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)

    return adata_ref, adata_new

def test_representation(adatas):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    ing = sc.tl.Ingest(adata_ref)
    ing.transform(adata_new)

    assert ing._use_rep == 'X_pca'
    assert ing._obsm['rep'].shape == (adata_new.n_obs, N_PCS)
    assert ing._pca_centered

    sc.pp.pca(adata_ref, n_comps=30, zero_center=False)
    sc.pp.neighbors(adata_ref)

    ing = sc.tl.Ingest(adata_ref)
    ing.transform(adata_new)

    assert ing._use_rep  = 'X_pca'
    assert ing._obsm['rep'].shape == (adata_new.n_obs, 30)
    assert not ing._pca_centered

    sc.pp.neighbors(adata_ref, use_rep='X')

    ing = sc.tl.Ingest(adata_ref)
    ing.transform(adata_new)

    assert ing._use_rep == 'X'
    assert ing._obsm['rep'] is adata_new.X

def test_neighbors(adatas):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    ing = sc.tl.Ingest(adata_ref)
    ing.transform(adata_new)
    ing.neighbors(k=10)
    indices = ing._indices

    tree = KDTree(adata_ref.obsm['X_pca'])
    true_indices = tree.query(ing._obsm['rep'], 10, return_distance=False)

    num_correct = 0.0
    for i in range(adata_new.n_obs):
        num_correct += np.sum(np.in1d(true_indices[i], indices[i]))
    percent_correct = num_correct / (adata_new.n_obs * 10)

    assert percent_correct > 0.99

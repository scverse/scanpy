import scanpy as sc
import numpy as np

X = np.array([[1, 0], [3, 0], [5, 6], [0, 4]])
n_neighbors = 3
key = 'test'


def test_neighbors_key_added():
    adata = sc.AnnData(X)

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0, key_added=key)

    conns_key = adata.uns[key]['connectivities_key']
    dists_key = adata.uns[key]['distances_key']

    assert adata.uns['neighbors']['params'] == adata.uns[key]['params']
    assert np.allclose(adata.obsp['connectivities'].toarray(), adata.obsp[conns_key].toarray())
    assert np.allclose(adata.obsp['distances'].toarray(), adata.obsp[dists_key].toarray())

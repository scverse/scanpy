import numpy as np
import scanpy as sc


def test_pca_sparse():
    pbmc = sc.datasets.pbmc3k()
    pbmc.X = pbmc.X.astype(np.float64)
    sc.pp.log1p(pbmc)

    implicit = sc.pp.pca(pbmc, pca_sparse=True, dtype=np.float64, copy=True)
    explicit = sc.pp.pca(pbmc, pca_sparse=False, dtype=np.float64, copy=True)

    assert np.allclose(implicit.uns["pca"]["variance"], explicit.uns["pca"]["variance"])
    assert np.allclose(
        implicit.uns["pca"]["variance_ratio"], explicit.uns["pca"]["variance_ratio"]
    )
    assert np.allclose(implicit.obsm['X_pca'], explicit.obsm['X_pca'])

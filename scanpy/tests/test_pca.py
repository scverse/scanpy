import pytest
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy import sparse

import scanpy as sc
from scanpy.tests.fixtures import array_type, float_dtype
from anndata.tests.helpers import assert_equal

A_list = [
    [0, 0, 7, 0, 0],
    [8, 5, 0, 2, 0],
    [6, 0, 0, 2, 5],
    [0, 0, 0, 1, 0],
    [8, 8, 2, 1, 0],
    [0, 0, 0, 4, 5],
]

A_pca = np.array(
    [
        [-4.4783009, 5.55508466, 1.73111572, -0.06029139, 0.17292555],
        [5.4855141, -0.42651191, -0.74776055, -0.74532146, 0.74633582],
        [0.01161428, -4.0156662, 2.37252748, -1.33122372, -0.29044446],
        [-3.61934397, 0.48525412, -2.96861931, -1.16312545, -0.33230607],
        [7.14050048, 1.86330409, -0.05786325, 1.25045782, -0.50213107],
        [-4.53998399, -3.46146476, -0.32940009, 2.04950419, 0.20562023],
    ]
)

A_svd = np.array(
    [
        [-0.77034038, -2.00750922, 6.64603489, -0.39669256, -0.22212097],
        [-9.47135856, -0.6326006, -1.33787112, -0.24894361, -1.02044665],
        [-5.90007339, 4.99658727, 0.70712592, -2.15188849, 0.30430008],
        [-0.19132409, 0.42172251, 0.11169531, 0.50977966, -0.71637566],
        [-11.1286238, -2.73045559, 0.08040596, 1.06850585, 0.74173764],
        [-1.50180389, 5.56886849, 1.64034442, 2.24476032, -0.05109001],
    ]
)


def test_pca_transform(array_type):
    A = array_type(A_list).astype('float32')
    A_pca_abs = np.abs(A_pca)
    A_svd_abs = np.abs(A_svd)

    adata = AnnData(A)

    sc.pp.pca(adata, n_comps=4, zero_center=True, svd_solver='arpack', dtype='float64')

    assert np.linalg.norm(A_pca_abs[:, :4] - np.abs(adata.obsm['X_pca'])) < 2e-05

    sc.pp.pca(
        adata,
        n_comps=5,
        zero_center=True,
        svd_solver='randomized',
        dtype='float64',
        random_state=14,
    )
    assert np.linalg.norm(A_pca_abs - np.abs(adata.obsm['X_pca'])) < 2e-05

    sc.pp.pca(adata, n_comps=4, zero_center=False, dtype='float64', random_state=14)
    assert np.linalg.norm(A_svd_abs[:, :4] - np.abs(adata.obsm['X_pca'])) < 2e-05


def test_pca_shapes():
    """Tests that n_comps behaves correctly"""
    # https://github.com/theislab/scanpy/issues/1051
    adata = AnnData(np.random.randn(30, 20))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (30, 19)

    adata = AnnData(np.random.randn(20, 30))
    sc.pp.pca(adata)
    assert adata.obsm["X_pca"].shape == (20, 19)

    with pytest.raises(ValueError):
        sc.pp.pca(adata, n_comps=100)


def test_pca_sparse(pbmc3k_normalized):
    """
    Tests that implicitly centered pca on sparse arrays returns equivalent results to
    explicit centering on dense arrays.
    """
    pbmc = pbmc3k_normalized

    pbmc_dense = pbmc.copy()
    pbmc_dense.X = pbmc_dense.X.toarray()

    implicit = sc.pp.pca(pbmc, dtype=np.float64, copy=True)
    explicit = sc.pp.pca(pbmc_dense, dtype=np.float64, copy=True)

    assert np.allclose(implicit.uns["pca"]["variance"], explicit.uns["pca"]["variance"])
    assert np.allclose(
        implicit.uns["pca"]["variance_ratio"], explicit.uns["pca"]["variance_ratio"]
    )
    assert np.allclose(implicit.obsm['X_pca'], explicit.obsm['X_pca'])
    assert np.allclose(implicit.varm['PCs'], explicit.varm['PCs'])


# This will take a while to run, but irreproducibility may
# not show up for float32 unless the matrix is large enough
def test_pca_reproducible(pbmc3k_normalized, array_type, float_dtype):
    pbmc = pbmc3k_normalized
    pbmc.X = array_type(pbmc.X)

    a = sc.pp.pca(pbmc, copy=True, dtype=float_dtype, random_state=42)
    b = sc.pp.pca(pbmc, copy=True, dtype=float_dtype, random_state=42)
    c = sc.pp.pca(pbmc, copy=True, dtype=float_dtype, random_state=0)

    assert_equal(a, b)
    # Test that changing random seed changes result
    assert not np.array_equal(a.obsm["X_pca"], c.obsm["X_pca"])


def test_pca_chunked(pbmc3k_normalized):
    # https://github.com/theislab/scanpy/issues/1590
    # But also a more general test

    # Subsetting for speed of test
    pbmc = pbmc3k_normalized[::6].copy()
    pbmc.X = pbmc.X.astype(np.float64)
    chunked = sc.pp.pca(pbmc3k_normalized, chunked=True, copy=True)
    default = sc.pp.pca(pbmc3k_normalized, copy=True)

    # Taking absolute value since sometimes dimensions are flipped
    np.testing.assert_allclose(
        np.abs(chunked.obsm["X_pca"]), np.abs(default.obsm["X_pca"])
    )
    np.testing.assert_allclose(np.abs(chunked.varm["PCs"]), np.abs(default.varm["PCs"]))
    np.testing.assert_allclose(
        np.abs(chunked.uns["pca"]["variance"]), np.abs(default.uns["pca"]["variance"])
    )
    np.testing.assert_allclose(
        np.abs(chunked.uns["pca"]["variance_ratio"]),
        np.abs(default.uns["pca"]["variance_ratio"]),
    )

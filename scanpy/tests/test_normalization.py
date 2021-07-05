import pytest
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy import sparse
import warnings

import scanpy as sc
from scanpy.tests.helpers import (
    check_rep_mutation,
    check_rep_results,
    _prepare_pbmc_testdata,
)
from anndata.tests.helpers import assert_equal, asarray


X_total = [[1, 0], [3, 0], [5, 6]]
X_frac = [[1, 0, 1], [3, 0, 1], [5, 6, 1]]


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    sc.pp.normalize_total(adata, key_added='n_counts')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [3.0, 3.0, 3.0])
    sc.pp.normalize_total(adata, target_sum=1, key_added='n_counts2')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [1.0, 1.0, 1.0])

    adata = AnnData(typ(X_frac, dtype=dtype))
    sc.pp.normalize_total(adata, exclude_highly_expressed=True, max_fraction=0.7)
    assert np.allclose(np.ravel(adata.X[:, 1:3].sum(axis=1)), [1.0, 1.0, 1.0])


@pytest.mark.parametrize('typ', [asarray, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_rep(typ, dtype):
    # Test that layer kwarg works
    X = typ(sparse.random(100, 50, format="csr", density=0.2, dtype=dtype))
    check_rep_mutation(sc.pp.normalize_total, X, fields=["layer"])
    check_rep_results(sc.pp.normalize_total, X, fields=["layer"])


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_layers(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    adata.layers["layer"] = adata.X.copy()
    with pytest.warns(FutureWarning, match=r".*layers.*deprecated"):
        sc.pp.normalize_total(adata, layers=["layer"])
    assert np.allclose(adata.layers["layer"].sum(axis=1), [3.0, 3.0, 3.0])


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_view(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    v = adata[:, :]

    sc.pp.normalize_total(v)
    sc.pp.normalize_total(adata)

    assert not v.is_view
    assert_equal(adata, v)


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_pearson_residuals_inputchecks(sparsity_func, dtype):

    adata = _prepare_pbmc_testdata(sparsity_func, dtype)

    # depending on check_values, warnings should be raised for non-integer data
    if dtype == 'float32':

        adata_noninteger = adata.copy()
        x, y = np.nonzero(adata_noninteger.X)
        adata_noninteger.X[x[0], y[0]] = 0.5

        with pytest.warns(UserWarning) as record:
            sc.pp.normalize_pearson_residuals(
                adata_noninteger.copy(), check_values=True
            )
        assert len(record) == 1
        assert "expects raw count data" in record[0].message.args[0]

        with pytest.warns(None) as record:
            sc.pp.normalize_pearson_residuals(
                adata_noninteger.copy(), check_values=False
            )
        assert len(record) == 0

    # errors should be raised for invalid theta values
    with pytest.raises(ValueError) as record:
        sc.pp.normalize_pearson_residuals(adata.copy(), theta=0)
    with pytest.raises(ValueError) as record:
        sc.pp.normalize_pearson_residuals(adata.copy(), theta=-1)

    # error should be raised for invalid clipping values
    with pytest.raises(ValueError) as record:
        sc.pp.normalize_pearson_residuals(adata.copy(), clip=-1)


@pytest.mark.parametrize(
    'sparsity_func', [np.array, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
@pytest.mark.parametrize('theta', [0.01, 1, 100, np.Inf])
@pytest.mark.parametrize('clip', [None, 1, np.Inf])
@pytest.mark.parametrize('inplace', [True, False])
def test_normalize_pearson_residuals_values(sparsity_func, dtype, theta, clip, inplace):

    # toy data
    X = np.array([[3, 6], [2, 4], [1, 0]])
    ns = np.sum(X, axis=1)
    ps = np.sum(X, axis=0) / np.sum(X)
    mu = np.outer(ns, ps)

    # compute reference residuals
    if np.isinf(theta):
        # Poisson case
        residuals_reference = (X - mu) / np.sqrt(mu)
    else:
        # NB case
        residuals_reference = (X - mu) / np.sqrt(mu + mu ** 2 / theta)

    # compute output to test
    adata = AnnData(sparsity_func(X), dtype=dtype)
    output = sc.pp.normalize_pearson_residuals(
        adata, theta=theta, clip=clip, inplace=inplace
    )

    # handle and test inplace argument
    if inplace:
        output_X = adata.X
        assert output is None
        # check for correct new `adata.uns` keys
        assert np.all(
            np.isin(['pearson_residuals_normalization'], list(adata.uns.keys()))
        )
        assert np.all(
            np.isin(
                ['theta', 'clip', 'computed_on'],
                list(adata.uns['pearson_residuals_normalization'].keys()),
            )
        )

    else:
        output_X = output['X']

    if clip is None:
        # default clipping: compare to sqrt(n) threshold
        clipping_threshold = np.sqrt(adata.shape[0]).astype(np.float32)
        assert np.max(output_X) <= clipping_threshold
        assert np.min(output_X) >= -clipping_threshold
    elif np.isinf(clip):
        # no clipping: compare to raw residuals
        assert np.allclose(output_X, residuals_reference)
    else:
        # custom clipping: compare to custom threshold
        assert np.max(output_X) <= clip
        assert np.min(output_X) >= -clip


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_pearson_residuals_pca(sparsity_func, dtype):

    adata = _prepare_pbmc_testdata(sparsity_func, dtype, small=True)
    n_cells = adata.shape[0]
    n_genes = adata.shape[1]
    n_hvgs = 100
    n_comps_pca = 50
    adata_with_hvgs = adata.copy()
    sc.pp.highly_variable_genes(
        adata_with_hvgs, flavor='pearson_residuals', n_top_genes=n_hvgs
    )
    adata_not_using_hvgs = adata_with_hvgs.copy()

    ### inplace = False ###
    # outputs the (potentially hvg-restricted) adata_pca object
    # PCA on all genes
    adata_pca = sc.pp.normalize_pearson_residuals_pca(
        adata.copy(), inplace=False, n_comps_pca=n_comps_pca
    )
    # PCA on hvgs only
    adata_pca_with_hvgs = sc.pp.normalize_pearson_residuals_pca(
        adata_with_hvgs.copy(), inplace=False, n_comps_pca=n_comps_pca
    )
    # PCA again on all genes (hvg use supressed)
    adata_pca_not_using_hvgs = sc.pp.normalize_pearson_residuals_pca(
        adata_not_using_hvgs.copy(),
        inplace=False,
        n_comps_pca=n_comps_pca,
        use_highly_variable=False,
    )

    # for both cases, check adata_pca keys are complete
    for ad in [adata_pca, adata_pca_with_hvgs, adata_pca_not_using_hvgs]:
        assert np.all(
            np.isin(
                ['pearson_residuals_normalization', 'pearson_residuals_pca'],
                list(ad.uns.keys()),
            )
        )
        assert np.all(np.isin(['X_pearson_residuals_pca'], list(ad.obsm.keys())))
        assert np.all(np.isin(['PCs'], list(ad.varm.keys())))
        assert ad.obsm['X_pearson_residuals_pca'].shape == (n_cells, n_comps_pca)

    # check adata shape to see if all genes or only HVGs are in the returned adata
    assert adata_pca.shape == (n_cells, n_genes)
    assert adata_pca_with_hvgs.shape == (n_cells, n_hvgs)  # only HVGs retained
    assert adata_pca_not_using_hvgs.shape == (n_cells, n_genes)

    # check PC shapes to see whether or not HVGs were used for PCA
    assert adata_pca.varm['PCs'].shape == (n_genes, n_comps_pca)
    assert adata_pca_with_hvgs.varm['PCs'].shape == (
        n_hvgs,
        n_comps_pca,
    )  # only HVGs used
    assert adata_pca_not_using_hvgs.varm['PCs'].shape == (n_genes, n_comps_pca)

    ### inplace = True ###
    # modifies the input adata object
    # PCA on all genes
    sc.pp.normalize_pearson_residuals_pca(adata, inplace=True, n_comps_pca=n_comps_pca)
    # PCA on hvgs only
    sc.pp.normalize_pearson_residuals_pca(
        adata_with_hvgs, inplace=True, n_comps_pca=n_comps_pca
    )
    # PCA again on all genes (hvg use supressed)
    sc.pp.normalize_pearson_residuals_pca(
        adata_not_using_hvgs,
        inplace=True,
        n_comps_pca=n_comps_pca,
        use_highly_variable=False,
    )

    for ad in [adata, adata_with_hvgs, adata_not_using_hvgs]:
        # check adata_pca keys are complete
        assert np.all(
            np.isin(
                ['pearson_residuals_normalization', 'pearson_residuals_pca'],
                list(ad.uns.keys()),
            )
        )
        assert np.all(np.isin(['X_pearson_residuals_pca'], list(ad.obsm.keys())))
        # check shapes: adata should always retains original shape
        assert ad.shape == (n_cells, n_genes)
        assert ad.obsm['X_pearson_residuals_pca'].shape == (n_cells, n_comps_pca)

    # check PC shapes to see whether or not HVGs were used for PCA
    assert adata.uns['pearson_residuals_pca']['PCs'].shape == (n_genes, n_comps_pca)
    assert adata_with_hvgs.uns['pearson_residuals_pca']['PCs'].shape == (
        n_hvgs,
        n_comps_pca,
    )
    assert adata_not_using_hvgs.uns['pearson_residuals_pca']['PCs'].shape == (
        n_genes,
        n_comps_pca,
    )


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_pearson_residuals_recipe(sparsity_func, dtype):
    adata = _prepare_pbmc_testdata(sparsity_func, dtype, small=True)
    n_cells = adata.shape[0]
    n_genes = adata.shape[1]
    n_hvgs = 100
    n_comps_pca = 50
    adata_with_hvgs = adata.copy()
    sc.pp.highly_variable_genes(
        adata_with_hvgs, flavor='pearson_residuals', n_top_genes=n_hvgs
    )

    ### inplace = False ###
    # outputs the (potentially hvg-restricted) adata_pca object
    # PCA on all genes
    adata_pca, hvg = sc.pp.recipe_pearson_residuals(
        adata.copy(), inplace=False, n_comps_pca=n_comps_pca, n_top_genes=n_hvgs
    )

    # for both cases, check adata_pca keys are complete
    assert np.all(
        np.isin(
            ['pearson_residuals_normalization', 'pearson_residuals_pca'],
            list(adata_pca.uns.keys()),
        )
    )
    assert np.all(np.isin(['X_pearson_residuals_pca'], list(adata_pca.obsm.keys())))
    assert np.all(np.isin(['PCs'], list(adata_pca.varm.keys())))
    assert adata_pca.obsm['X_pearson_residuals_pca'].shape == (n_cells, n_comps_pca)

    # check adata shape
    assert adata_pca.shape == (n_cells, n_hvgs)
    # check PC shapes to check that HVGs were used for PCA
    assert adata_pca.varm['PCs'].shape == (n_hvgs, n_comps_pca)

    # check hvg df
    assert np.all(
        np.isin(
            [
                'means',
                'variances',
                'residual_variances',
                'highly_variable_rank',
                'highly_variable',
            ],
            list(hvg.columns),
        )
    )
    assert np.sum(hvg['highly_variable']) == n_hvgs
    assert hvg.shape[0] == n_genes

    ### inplace = True ###
    # modifies the input adata object
    # PCA on all genes
    sc.pp.recipe_pearson_residuals(
        adata, inplace=True, n_comps_pca=n_comps_pca, n_top_genes=n_hvgs
    )

    assert np.all(
        np.isin(
            ['pearson_residuals_normalization', 'pearson_residuals_pca'],
            list(adata.uns.keys()),
        )
    )
    assert np.all(np.isin(['X_pearson_residuals_pca'], list(adata.obsm.keys())))
    assert adata.shape == (n_cells, n_genes)
    assert adata.obsm['X_pearson_residuals_pca'].shape == (n_cells, n_comps_pca)

    # check PC shapes to see whether or not HVGs were used for PCA
    assert adata.uns['pearson_residuals_pca']['PCs'].shape == (n_hvgs, n_comps_pca)

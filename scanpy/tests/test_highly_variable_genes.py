import pytest
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from scipy.sparse import csr_matrix
from scanpy.tests.helpers import _prepare_pbmc_testdata
import warnings

FILE = Path(__file__).parent / Path('_scripts/seurat_hvg.csv')
FILE_V3 = Path(__file__).parent / Path('_scripts/seurat_hvg_v3.csv.gz')
FILE_V3_BATCH = Path(__file__).parent / Path('_scripts/seurat_hvg_v3_batch.csv')


def test_highly_variable_genes_basic():
    adata = sc.datasets.blobs()
    sc.pp.highly_variable_genes(adata)

    adata = sc.datasets.blobs()
    np.random.seed(0)
    adata.obs['batch'] = np.random.binomial(3, 0.5, size=(adata.n_obs))
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    sc.pp.highly_variable_genes(adata, batch_key='batch')
    assert 'highly_variable_nbatches' in adata.var.columns
    assert 'highly_variable_intersection' in adata.var.columns

    adata = sc.datasets.blobs()
    adata.obs['batch'] = np.random.binomial(4, 0.5, size=(adata.n_obs))
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    sc.pp.highly_variable_genes(adata, batch_key='batch', n_top_genes=3)
    assert 'highly_variable_nbatches' in adata.var.columns
    assert adata.var['highly_variable'].sum() == 3

    sc.pp.highly_variable_genes(adata)
    no_batch_hvg = adata.var.highly_variable.copy()
    assert no_batch_hvg.any()
    adata.obs['batch'] = 'batch'
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    sc.pp.highly_variable_genes(adata, batch_key='batch')
    assert np.all(no_batch_hvg == adata.var.highly_variable)
    assert np.all(adata.var.highly_variable_intersection == adata.var.highly_variable)

    adata.obs["batch"] = "a"
    adata.obs.batch.loc[::2] = "b"
    sc.pp.highly_variable_genes(adata, batch_key="batch")
    assert adata.var["highly_variable"].any()

    colnames = [
        'means',
        'dispersions',
        'dispersions_norm',
        'highly_variable_nbatches',
        'highly_variable_intersection',
        'highly_variable',
    ]
    hvg_df = sc.pp.highly_variable_genes(adata, batch_key="batch", inplace=False)
    assert np.all(np.isin(colnames, hvg_df.columns))


def _residual_var_reference(adata, clip=None, theta=100):
    sc.experimental.pp.normalize_pearson_residuals(adata, clip=clip, theta=theta)
    residuals = adata.X
    return np.var(residuals, axis=0)


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_highly_variable_genes_pearson_residuals_inputchecks(sparsity_func, dtype):

    adata = _prepare_pbmc_testdata(sparsity_func, dtype, small=True)

    # depending on check_values, warnings should be raised for non-integer data
    if dtype == 'float32':

        adata_noninteger = adata.copy()
        x, y = np.nonzero(adata_noninteger.X)
        adata_noninteger.X[x[0], y[0]] = 0.5

        # expecting 0 no-int warnings
        with warnings.catch_warnings(record=True) as record:
            sc.experimental.pp.highly_variable_genes(
                adata_noninteger.copy(),
                flavor='pearson_residuals',
                n_top_genes=100,
                check_values=False,
            )

        warning_msgs = [w.message.args[0] for w in record]
        assert (
            "`flavor='pearson_residuals'` expects raw count data, but non-integers were found."
            not in warning_msgs
        )

        # expecting 1 no-int warning
        with pytest.warns(
            UserWarning,
            match="`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
        ) as record:
            sc.experimental.pp.highly_variable_genes(
                adata_noninteger.copy(),
                flavor='pearson_residuals',
                n_top_genes=100,
                check_values=True,
            )

    # errors should be raised for invalid theta values
    with pytest.raises(
        ValueError, match='Pearson residuals require theta > 0'
    ) as record:
        sc.experimental.pp.highly_variable_genes(
            adata.copy(), flavor='pearson_residuals', n_top_genes=100, theta=0
        )

    with pytest.raises(
        ValueError, match='Pearson residuals require theta > 0'
    ) as record:
        sc.experimental.pp.highly_variable_genes(
            adata.copy(), flavor='pearson_residuals', n_top_genes=100, theta=-1
        )

    # error should be raised for invalid clipping values
    with pytest.raises(
        ValueError, match='Pearson residuals require `clip>=0` or `clip=None`.'
    ) as record:
        sc.experimental.pp.highly_variable_genes(
            adata.copy(), flavor='pearson_residuals', n_top_genes=100, clip=-1
        )


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
@pytest.mark.parametrize('subset', [True, False])
@pytest.mark.parametrize('clip', [None, np.Inf, 30])
@pytest.mark.parametrize('theta', [100, np.Inf])
@pytest.mark.parametrize('n_top_genes', [100, 200])
def test_highly_variable_genes_pearson_residuals_values(
    subset, sparsity_func, dtype, clip, theta, n_top_genes
):
    adata = _prepare_pbmc_testdata(sparsity_func, dtype, small=True)
    # cleanup var
    adata.var.drop(columns=adata.var.columns, inplace=True)
    # compute reference output
    residual_variances_reference = _residual_var_reference(
        adata.copy(), clip=clip, theta=theta
    )
    if subset:
        # lazyly sort by residual variance and take top N
        top_n_idx = np.argsort(-residual_variances_reference)[:n_top_genes]
        # (results in sorted "gene order" in reference)
        residual_variances_reference = residual_variances_reference[top_n_idx]

    # compute output to be tested
    output_df = sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=False,
        clip=clip,
        theta=theta,
    )

    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=True,
        clip=clip,
        theta=theta,
    )

    pd.testing.assert_frame_equal(output_df, adata.var)

    # consistency with normalization method
    if subset:
        # sort values before comparing as reference is sorted as well for subset case
        sort_output_idx = np.argsort(-output_df['residual_variances'].values)
        assert np.allclose(
            output_df['residual_variances'].values[sort_output_idx],
            residual_variances_reference,
        )
    else:
        assert np.allclose(
            output_df['residual_variances'].values, residual_variances_reference
        )


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
@pytest.mark.parametrize('subset', [True, False])
@pytest.mark.parametrize('n_top_genes', [1000, 500])
def test_highly_variable_genes_pearson_residuals_general(
    subset, sparsity_func, dtype, n_top_genes
):

    adata = _prepare_pbmc_testdata(sparsity_func, dtype)
    # cleanup var
    adata.var.drop(columns=adata.var.columns, inplace=True)
    # compute reference output
    residual_variances_reference = _residual_var_reference(adata.copy())
    if subset:
        # lazily sort by residual variance and take top N
        top_n_idx = np.argsort(-residual_variances_reference)[:n_top_genes]
        # (results in sorted "gene order" in reference)
        residual_variances_reference = residual_variances_reference[top_n_idx]
    # compute output to be tested
    output_df = sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=False,
    )

    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=True,
    )

    pd.testing.assert_frame_equal(output_df, adata.var)

    # check output is complete
    for key in [
        'highly_variable',
        'means',
        'variances',
        'residual_variances',
        'highly_variable_rank',
    ]:
        assert key in output_df.keys()

    # check residual variances
    assert output_df['residual_variances'].values.dtype is np.dtype('float32')
    # consistency with normalization method
    if subset:
        # sort values before comparing as reference is sorted as well for subset case
        sort_output_idx = np.argsort(-output_df['residual_variances'].values)
        assert np.allclose(
            output_df['residual_variances'].values[sort_output_idx],
            residual_variances_reference,
        )
    else:
        assert np.allclose(
            output_df['residual_variances'].values, residual_variances_reference
        )

    # check hvg flag
    assert output_df['highly_variable'].values.dtype is np.dtype('bool')
    assert np.sum(output_df['highly_variable']) == n_top_genes
    hvg_idx = np.where(output_df['highly_variable'])[0]
    topn_idx = np.sort(
        np.argsort(-output_df['residual_variances'].values)[:n_top_genes]
    )
    assert np.all(hvg_idx == topn_idx)

    # check ranks
    assert np.nanmin(output_df['highly_variable_rank'].values) == 0
    assert np.nanmax(output_df['highly_variable_rank'].values) <= n_top_genes - 1


@pytest.mark.parametrize(
    'sparsity_func', [csr_matrix.toarray, csr_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
@pytest.mark.parametrize('subset', [True, False])
@pytest.mark.parametrize('n_top_genes', [1000, 500])
def test_highly_variable_genes_pearson_residuals_batch(
    subset, n_top_genes, sparsity_func, dtype
):
    adata = _prepare_pbmc_testdata(sparsity_func, dtype)
    # cleanup var
    adata.var.drop(columns=adata.var.columns, inplace=True)
    n_genes = adata.shape[1]

    output_df = sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        batch_key='batch',
        subset=subset,
        inplace=False,
    )

    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor='pearson_residuals',
        n_top_genes=n_top_genes,
        batch_key='batch',
        subset=subset,
        inplace=True,
    )

    pd.testing.assert_frame_equal(output_df, adata.var)

    # check output is complete
    for key in [
        'highly_variable',
        'means',
        'variances',
        'residual_variances',
        'highly_variable_rank',
        'highly_variable_nbatches',
        'highly_variable_intersection',
    ]:
        assert key in output_df.keys()

    # check hvg flag
    assert output_df['highly_variable'].values.dtype is np.dtype('bool')
    assert np.sum(output_df['highly_variable']) == n_top_genes

    # check intersection flag
    nbatches = len(np.unique(adata.obs['batch']))
    assert output_df['highly_variable_intersection'].values.dtype is np.dtype('bool')
    assert np.sum(output_df['highly_variable_intersection']) <= n_top_genes * nbatches
    assert np.all(output_df['highly_variable'][output_df.highly_variable_intersection])

    # check ranks (with batch_key these are the median of within-batch ranks)
    assert output_df['highly_variable_rank'].values.dtype is np.dtype('float32')
    assert np.nanmax(output_df['highly_variable_rank'].values) <= n_top_genes - 1

    # check nbatches
    assert output_df['highly_variable_nbatches'].values.dtype is np.dtype('int')
    assert np.min(output_df['highly_variable_nbatches'].values) >= 0
    assert np.max(output_df['highly_variable_nbatches'].values) <= nbatches

    # check subsetting
    if subset:
        assert len(output_df) == n_top_genes
    else:
        assert len(output_df) == n_genes


def test_higly_variable_genes_compare_to_seurat():
    seurat_hvg_info = pd.read_csv(FILE, sep=' ')

    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    pbmc.var_names_make_unique()

    sc.pp.normalize_per_cell(pbmc, counts_per_cell_after=1e4)
    sc.pp.log1p(pbmc)
    sc.pp.highly_variable_genes(
        pbmc, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True
    )

    np.testing.assert_array_equal(
        seurat_hvg_info['highly_variable'], pbmc.var['highly_variable']
    )

    # (still) Not equal to tolerance rtol=2e-05, atol=2e-05
    # np.testing.assert_allclose(4, 3.9999, rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(
        seurat_hvg_info['means'],
        pbmc.var['means'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions'],
        pbmc.var['dispersions'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions_norm'],
        pbmc.var['dispersions_norm'],
        rtol=2e-05,
        atol=2e-05,
    )


def test_higly_variable_genes_compare_to_seurat_v3():
    seurat_hvg_info = pd.read_csv(
        FILE_V3, sep=' ', dtype={"variances_norm": np.float64}
    )

    pbmc = sc.datasets.pbmc3k()
    pbmc.var_names_make_unique()

    pbmc_dense = pbmc.copy()
    pbmc_dense.X = pbmc_dense.X.toarray()

    sc.pp.highly_variable_genes(pbmc, n_top_genes=1000, flavor='seurat_v3')
    sc.pp.highly_variable_genes(pbmc_dense, n_top_genes=1000, flavor='seurat_v3')

    np.testing.assert_array_equal(
        seurat_hvg_info['highly_variable'], pbmc.var['highly_variable']
    )
    np.testing.assert_allclose(
        seurat_hvg_info['variances'],
        pbmc.var['variances'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['variances_norm'],
        pbmc.var['variances_norm'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        pbmc_dense.var['variances_norm'],
        pbmc.var['variances_norm'],
        rtol=2e-05,
        atol=2e-05,
    )

    batch = np.zeros((len(pbmc)), dtype=int)
    batch[1500:] = 1
    pbmc.obs["batch"] = batch
    df = sc.pp.highly_variable_genes(
        pbmc, n_top_genes=4000, flavor='seurat_v3', batch_key="batch", inplace=False
    )
    df.sort_values(
        ["highly_variable_nbatches", "highly_variable_rank"],
        ascending=[False, True],
        na_position="last",
        inplace=True,
    )
    df = df.iloc[:4000]
    seurat_hvg_info_batch = pd.read_csv(
        FILE_V3_BATCH, sep=' ', dtype={"variances_norm": np.float64}
    )

    # ranks might be slightly different due to many genes having same normalized var
    seu = pd.Index(seurat_hvg_info_batch['x'].values)
    assert len(seu.intersection(df.index)) / 4000 > 0.95

    sc.pp.log1p(pbmc)
    with pytest.warns(
        UserWarning,
        match="`flavor='seurat_v3'` expects raw count data, but non-integers were found.",
    ):
        sc.pp.highly_variable_genes(pbmc, n_top_genes=1000, flavor='seurat_v3')


def test_filter_genes_dispersion_compare_to_seurat():
    seurat_hvg_info = pd.read_csv(FILE, sep=' ')

    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    pbmc.var_names_make_unique()

    sc.pp.normalize_per_cell(pbmc, counts_per_cell_after=1e4)
    sc.pp.filter_genes_dispersion(
        pbmc,
        flavor='seurat',
        log=True,
        subset=False,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
    )

    np.testing.assert_array_equal(
        seurat_hvg_info['highly_variable'], pbmc.var['highly_variable']
    )

    # (still) Not equal to tolerance rtol=2e-05, atol=2e-05:
    # np.testing.assert_allclose(4, 3.9999, rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(
        seurat_hvg_info['means'],
        pbmc.var['means'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions'],
        pbmc.var['dispersions'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions_norm'],
        pbmc.var['dispersions_norm'],
        rtol=2e-05,
        atol=2e-05,
    )


def test_highly_variable_genes_batches():
    adata = sc.datasets.pbmc68k_reduced()
    adata[:100, :100].X = np.zeros((100, 100))

    adata.obs['batch'] = ['0' if i < 100 else '1' for i in range(adata.n_obs)]
    adata_1 = adata[adata.obs.batch.isin(['0']), :]
    adata_2 = adata[adata.obs.batch.isin(['1']), :]

    sc.pp.highly_variable_genes(
        adata,
        batch_key='batch',
        flavor='cell_ranger',
        n_top_genes=200,
    )

    sc.pp.filter_genes(adata_1, min_cells=1)
    sc.pp.filter_genes(adata_2, min_cells=1)
    hvg1 = sc.pp.highly_variable_genes(
        adata_1, flavor='cell_ranger', n_top_genes=200, inplace=False
    )
    hvg2 = sc.pp.highly_variable_genes(
        adata_2, flavor='cell_ranger', n_top_genes=200, inplace=False
    )

    assert np.isclose(
        adata.var['dispersions_norm'][100],
        0.5 * hvg1['dispersions_norm'][0] + 0.5 * hvg2['dispersions_norm'][100],
    )
    assert np.isclose(
        adata.var['dispersions_norm'][101],
        0.5 * hvg1['dispersions_norm'][1] + 0.5 * hvg2['dispersions_norm'][101],
    )
    assert np.isclose(
        adata.var['dispersions_norm'][0], 0.5 * hvg2['dispersions_norm'][0]
    )

    colnames = [
        'means',
        'dispersions',
        'dispersions_norm',
        'highly_variable',
    ]

    assert np.all(np.isin(colnames, hvg1.columns))


from scanpy.preprocessing._utils import _get_mean_var


def test_seurat_v3_mean_var_output_with_batchkey():
    pbmc = sc.datasets.pbmc3k()
    pbmc.var_names_make_unique()
    n_cells = pbmc.shape[0]
    batch = np.zeros((n_cells), dtype=int)
    batch[1500:] = 1
    pbmc.obs["batch"] = batch

    # true_mean, true_var = _get_mean_var(pbmc.X)
    true_mean = np.mean(pbmc.X.toarray(), axis=0)
    true_var = np.var(pbmc.X.toarray(), axis=0, dtype=np.float64, ddof=1)

    result_df = sc.pp.highly_variable_genes(
        pbmc, batch_key='batch', flavor='seurat_v3', n_top_genes=4000, inplace=False
    )
    np.testing.assert_allclose(true_mean, result_df['means'], rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(true_var, result_df['variances'], rtol=2e-05, atol=2e-05)

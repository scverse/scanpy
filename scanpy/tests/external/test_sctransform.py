import scanpy as sc
import scanpy.external as sce

def test_sct():
    ad = sc.datasets.paul15()[:200, :2001].copy()

    sce.pp.sctransform(ad, verbose=False)
    assert ad.var.highly_variable.sum() > 1000
    assert 'highly_variable_sct_residual_var' in ad.var_keys()
    assert 'sct_corrected' in ad.layers

    ad.obs['batch'] = 'A'
    ad.obs.batch.iloc[:100] = 'B'
    ad.var.drop(
        ['highly_variable', 'highly_variable_sct_residual_var'],
        axis=1,
        inplace=True,
    )
    ad = sce.pp.sctransform(
        ad,
        batch_key='batch',
        store_residuals=True,
        verbose=False,
        inplace=False,
    )

    assert ad.var.highly_variable.sum() > 1000
    assert 'highly_variable_sct_residual_var' in ad.var_keys()
    assert 'sct_corrected' in ad.layers
    assert 'sct_residuals' in ad.layers

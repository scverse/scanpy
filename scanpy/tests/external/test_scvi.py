import pytest
import numpy as numpy
import scanpy.external as sce
import numpy as np

from anndata import AnnData


def test_scvi():
    n_samples = 4
    batch1 = np.random.randint(1, 5, size=(n_samples, 7))
    batch2 = np.random.randint(1, 5, size=(n_samples, 7))
    ad1 = AnnData(batch1)
    ad2 = AnnData(batch2)
    adata = ad1.concatenate(ad2, batch_categories=['test1', 'test2'])
    n_latent = 30
    sce.pp.scvi(
        adata,
        use_cuda=False,
        n_epochs=1,
        n_latent=n_latent,
        return_posterior=True,
        batch_key='batch',
    )
    assert adata.obsm['X_scvi'].shape == (n_samples * 2, n_latent)
    assert adata.obsm['X_scvi_denoised'].shape == adata.shape
    assert adata.obsm['X_scvi_sample_rate'].shape == adata.shape

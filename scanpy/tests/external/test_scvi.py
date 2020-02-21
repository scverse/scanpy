import pytest
import numpy as numpy
import scanpy.external as sce
import numpy as np

from anndata import AnnData


def test_scvi():
    n_samples = 4
    raw_data = np.random.randint(1, 5, size=(n_samples, 7))
    ad = AnnData(raw_data)
    ad.raw = ad.copy()
    n_latent = 30
    ad = sce.pp.scvi(ad, use_cuda=False, n_epochs=1, n_latent=n_latent)
    assert ad.obsm['X_scvi'].shape == (n_samples, n_latent)
    assert ad.obsm['X_scvi_denoised'].shape == raw_data.shape
    assert ad.obsm['X_scvi_sample_rate'].shape == raw_data.shape

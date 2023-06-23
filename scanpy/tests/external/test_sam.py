import pytest
import scanpy as sc
import scanpy.external as sce
import numpy as np

from scanpy.testing._helpers.data import pbmc3k


pytest.importorskip("samalg")


def test_sam():
    adata_ref = pbmc3k()
    ix = np.random.choice(adata_ref.shape[0], size=200, replace=False)
    adata = adata_ref[ix, :].copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sce.tl.sam(adata, inplace=True)
    uns_keys = list(adata.uns.keys())
    obsm_keys = list(adata.obsm.keys())
    assert all(['sam' in uns_keys, 'X_umap' in obsm_keys, 'neighbors' in uns_keys])

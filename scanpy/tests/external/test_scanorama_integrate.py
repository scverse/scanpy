from importlib.util import find_spec

import pytest

import scanpy as sc
import scanpy.external as sce


@pytest.mark.skipif(not find_spec("scanorama"), reason="needs module `scanorama`")
def test_scanorama_integrate():
    """
    Test that Scanorama integration works.

    This is a very simple test that just checks to see if the Scanorama
    integrate wrapper succesfully added a new field to ``adata.obsm``
    and makes sure it has the same dimensions as the original PCA table.
    """
    adata = sc.datasets.pbmc68k_reduced()
    sc.tl.pca(adata)
    adata.obs['batch'] = 350 * ['a'] + 350 * ['b']
    sce.pp.scanorama_integrate(adata, 'batch', approx=False)
    assert adata.obsm['X_scanorama'].shape == adata.obsm['X_pca'].shape

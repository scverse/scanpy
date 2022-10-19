import pytest

import scanpy as sc
import scanpy.external as sce

pytest.importorskip("scalex")


def test_scalex_integrate():
    """
    Test that SCALEX integrate works.

    This is a very simple test that just checks to see if the SCALEX
    integrate wrapper succesfully added a new field to ``adata.obsm``
    """
    adata = sc.datasets.pbmc3k()

    adata.obs['batch'] = 1350 * ['a'] + 1350 * ['b']
    sce.pp.scalex_integrate(adata)
    assert adata.obsm['X_scalex_umap'].shape[0] == adata.shape[0]

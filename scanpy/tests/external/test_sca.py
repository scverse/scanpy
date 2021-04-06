import pytest
import scanpy as sc
import scanpy.external as sce

pytest.importorskip("shannonca")


def test_sca():
    """
    Test that SCA outputs the correct dimensionality
    """
    adata = sc.datasets.pbmc3k()
    sce.tl.sca(adata, iters=1, key_added='sca', n_comps=10)
    assert adata.obsm['X_sca'].shape[1] == 10, 'wrong dimension output!'

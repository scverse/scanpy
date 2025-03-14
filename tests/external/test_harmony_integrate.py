from __future__ import annotations

import scanpy as sc
import scanpy.external as sce
from testing.scanpy._helpers.data import pbmc3k
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.harmonypy]


def test_harmony_integrate():
    """Test that Harmony integrate works.

    This is a very simple test that just checks to see if the Harmony
    integrate wrapper succesfully added a new field to ``adata.obsm``
    and makes sure it has the same dimensions as the original PCA table.
    """
    adata = pbmc3k()
    sc.pp.recipe_zheng17(adata)
    sc.pp.pca(adata)
    adata.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]
    sce.pp.harmony_integrate(adata, "batch")
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape

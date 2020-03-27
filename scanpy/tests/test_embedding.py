from importlib.util import find_spec

import numpy as np
import pytest
from sklearn.utils.testing import assert_array_almost_equal

import scanpy as sc


def test_umap_init_dtype():
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc = pbmc[:100, :].copy()
    sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(np.float32))
    embed1 = pbmc.obsm["X_umap"].copy()
    sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(np.float64))
    embed2 = pbmc.obsm["X_umap"].copy()
    assert_array_almost_equal(embed1, embed2)
    assert_array_almost_equal(embed1, embed2)


needs_fa2 = pytest.mark.skipif(not find_spec("fa2"), reason="needs module `fa2`")


@pytest.mark.parametrize(
    "layout", [pytest.param("fa", marks=needs_fa2), "fr"],
)
def test_umap_init_paga(layout):
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc = pbmc[:100, :].copy()
    sc.tl.paga(pbmc)
    sc.pl.paga(pbmc, layout=layout, show=False)
    sc.tl.umap(pbmc, init_pos="paga")


def test_umap_args():
    pbmc = sc.datasets.pbmc68k_reduced()
    if "connectivities" in pbmc.uns["neighbors"]:  # Backwards compat for older anndata
        pbmc.obsp["connectivities"] = pbmc.uns["neighbors"].pop("connectivities")
        pbmc.obsp["distances"] = pbmc.uns["neighbors"].pop("distances")

    sc.tl.umap(pbmc, key_added="umap_orig")

    del pbmc.uns["neighbors"]

    with pytest.warns(RuntimeWarning):
        sc.tl.umap(pbmc, obsp="connectivities", key_added="umap_no-neighbors")
    pbmc.obsp["new_connect"] = pbmc.obsp.pop("connectivities")

    with pytest.warns(RuntimeWarning):
        sc.tl.umap(pbmc, obsp="new_connect", key_added="umap_obsp-key")

    assert np.allclose(pbmc.obsm["umap_orig"], pbmc.obsm["umap_no-neighbors"])
    assert np.allclose(pbmc.obsm["umap_orig"], pbmc.obsm["umap_obsp-key"])

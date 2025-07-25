from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

n_neighbors = 5
key = "test"


@pytest.fixture(scope="session")
def adata_session() -> sc.AnnData:
    adata = sc.AnnData(pbmc68k_reduced().X)
    sc.pp.pca(adata)
    return adata


@pytest.fixture
def adata(adata_session: sc.AnnData) -> sc.AnnData:
    return adata_session.copy()


def test_neighbors_key_added(adata: sc.AnnData) -> None:
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0, key_added=key)

    conns_key = adata.uns[key]["connectivities_key"]
    dists_key = adata.uns[key]["distances_key"]

    assert adata.uns["neighbors"]["params"] == adata.uns[key]["params"]
    assert np.allclose(
        adata.obsp["connectivities"].toarray(), adata.obsp[conns_key].toarray()
    )
    assert np.allclose(
        adata.obsp["distances"].toarray(), adata.obsp[dists_key].toarray()
    )


def test_neighbors_pca_keys_added_without_previous_pca_run(adata: sc.AnnData) -> None:
    del adata.uns["pca"]
    del adata.obsm["X_pca"]
    with pytest.warns(
        UserWarning,
        match=r".*Falling back to preprocessing with `sc.pp.pca` and default params",
    ):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0)
    assert "pca" in adata.uns


# test functions with neighbors_key and obsp
@needs.igraph
@pytest.mark.parametrize("field", ["neighbors_key", "obsp"])
def test_neighbors_key_obsp(adata, field):
    adata1 = adata.copy()

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=0)
    sc.pp.neighbors(adata1, n_neighbors=n_neighbors, random_state=0, key_added=key)

    if field == "neighbors_key":
        arg = {field: key}
    else:
        arg = {field: adata1.uns[key]["connectivities_key"]}

    sc.tl.draw_graph(adata, layout="fr", random_state=1)
    sc.tl.draw_graph(adata1, layout="fr", random_state=1, **arg)

    assert adata.uns["draw_graph"]["params"] == adata1.uns["draw_graph"]["params"]
    assert np.allclose(adata.obsm["X_draw_graph_fr"], adata1.obsm["X_draw_graph_fr"])

    sc.tl.leiden(adata, flavor="igraph", random_state=0)
    sc.tl.leiden(adata1, flavor="igraph", random_state=0, **arg)

    assert adata.uns["leiden"]["params"] == adata1.uns["leiden"]["params"]
    assert np.all(adata.obs["leiden"] == adata1.obs["leiden"])

    # no obsp in umap, paga
    if field == "neighbors_key":
        sc.tl.umap(adata, random_state=0)
        sc.tl.umap(adata1, random_state=0, neighbors_key=key)

        assert adata.uns["umap"]["params"] == adata1.uns["umap"]["params"]
        assert np.allclose(adata.obsm["X_umap"], adata1.obsm["X_umap"])

        sc.tl.paga(adata, groups="leiden")
        sc.tl.paga(adata1, groups="leiden", neighbors_key=key)

        assert np.allclose(
            adata.uns["paga"]["connectivities"].toarray(),
            adata1.uns["paga"]["connectivities"].toarray(),
        )
        assert np.allclose(
            adata.uns["paga"]["connectivities_tree"].toarray(),
            adata1.uns["paga"]["connectivities_tree"].toarray(),
        )

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


n_neighbors = 5
key = "test"
uns_key = f"{key}_neighbors"


@pytest.fixture(scope="session")
def adata_session() -> sc.AnnData:
    adata = sc.AnnData(pbmc68k_reduced().X)
    sc.pp.pca(adata)
    return adata


@pytest.fixture
def adata(adata_session: sc.AnnData) -> sc.AnnData:
    return adata_session.copy()


@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_neighbors_key_added(
    adata: sc.AnnData, rng_arg: Literal["rng", "random_state"]
) -> None:
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, **{rng_arg: 0})
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, **{rng_arg: 0}, out=sc.Prefix(key))

    conns_key = adata.uns[uns_key]["connectivities_key"]
    dists_key = adata.uns[uns_key]["distances_key"]

    assert adata.uns["neighbors"]["params"] == adata.uns[uns_key]["params"]
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


LEGACY = pytest.mark.filterwarnings(
    "ignore:.*argument (neighbors_key|obsp) is deprecated:FutureWarning"
)


@needs.igraph
@pytest.mark.parametrize(
    "field",
    [
        pytest.param("prefix", id="use-prefix"),
        pytest.param("obsp_acc", marks=needs.anndata_acc, id="use-obsp"),
        pytest.param("neighbors_key", marks=LEGACY, id="legacy-neighbors_key"),
        pytest.param("obsp", marks=LEGACY, id="legacy-obsp"),
    ],
)
@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_neighbors_key_obsp(
    subtests: pytest.Subtests,
    adata: AnnData,
    field: Literal["prefix", "obsp_acc", "neighbors_key", "obsp"],
    rng_arg: Literal["rng", "random_state"],
) -> None:
    """Functions using the neighbors graph find a non-default one via `use` (or the legacy args)."""
    adata1 = adata.copy()

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, **{rng_arg: 0})
    sc.pp.neighbors(adata1, n_neighbors=n_neighbors, **{rng_arg: 0}, out=sc.Prefix(key))
    conns_key = adata1.uns[uns_key]["connectivities_key"]

    match field:
        case "prefix":
            arg = dict(use=sc.Prefix(key))
        case "obsp_acc":
            from anndata.acc import A

            arg = dict(use=A.obsp[conns_key])
        case "neighbors_key":
            arg = dict(neighbors_key=uns_key)
        case "obsp":
            arg = dict(obsp=conns_key)
    # `umap` and `paga` need the `.uns` entry, not just a graph
    uns_arg = None if field in {"obsp_acc", "obsp"} else arg

    sc.tl.draw_graph(adata, layout="fr", **{rng_arg: 1})
    sc.tl.draw_graph(adata1, layout="fr", **{rng_arg: 1}, **arg)

    with subtests.test("draw_graph"):
        assert adata.uns["draw_graph"]["params"] == adata1.uns["draw_graph"]["params"]
        assert np.allclose(
            adata.obsm["X_draw_graph_fr"], adata1.obsm["X_draw_graph_fr"]
        )

    sc.tl.leiden(adata, flavor="igraph", **{rng_arg: 0})
    sc.tl.leiden(adata1, flavor="igraph", **{rng_arg: 0}, **arg)

    with subtests.test("leiden"):
        assert adata.uns["leiden"]["params"] == adata1.uns["leiden"]["params"]
        assert np.all(adata.obs["leiden"] == adata1.obs["leiden"])

    if uns_arg is not None:
        sc.tl.umap(adata, **{rng_arg: 0})
        sc.tl.umap(adata1, **{rng_arg: 0}, **uns_arg)

        with subtests.test("umap"):
            assert adata.uns["umap"]["params"] == adata1.uns["umap"]["params"]
            assert np.allclose(adata.obsm["X_umap"], adata1.obsm["X_umap"])

        sc.tl.paga(adata, groups="leiden")
        sc.tl.paga(adata1, groups="leiden", **uns_arg)

        with subtests.test("paga"):
            assert np.allclose(
                adata.uns["paga"]["connectivities"].toarray(),
                adata1.uns["paga"]["connectivities"].toarray(),
            )
            assert np.allclose(
                adata.uns["paga"]["connectivities_tree"].toarray(),
                adata1.uns["paga"]["connectivities_tree"].toarray(),
            )

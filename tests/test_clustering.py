from __future__ import annotations

from contextlib import nullcontext
from functools import partial
from typing import TYPE_CHECKING

import pandas as pd
import pytest
from sklearn.metrics.cluster import normalized_mutual_info_score

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


@pytest.fixture
def adata_neighbors() -> AnnData:
    return pbmc68k_reduced()


def _leiden_copy(adata: AnnData, **kwargs) -> AnnData:
    adata = adata.copy()
    sc.tl.leiden(adata, **kwargs)
    return adata


@pytest.fixture(
    params=[
        pytest.param("igraph", marks=needs.igraph),
        pytest.param("leidenalg", marks=needs.leidenalg),
    ]
)
def flavor(request: pytest.FixtureRequest) -> Literal["igraph", "leidenalg"]:
    return request.param


@needs.leidenalg
@needs.igraph
@pytest.mark.parametrize("resolution", [1, 2])
@pytest.mark.parametrize("n_iterations", [-1, 3])
def test_leiden_basic(
    adata_neighbors: AnnData,
    flavor: Literal["igraph", "leidenalg"],
    resolution: float,
    n_iterations: int,
) -> None:
    with (
        nullcontext()
        if flavor == "igraph"
        else pytest.warns(
            UserWarning, match=r"The `igraph` implementation of leiden clustering"
        )
    ):
        sc.tl.leiden(
            adata_neighbors,
            flavor=flavor,
            resolution=resolution,
            n_iterations=n_iterations,
            directed=(flavor == "leidenalg"),
            out=sc.Prefix("custom"),
        )
    assert "custom_leiden" in adata_neighbors.obs
    assert adata_neighbors.uns["custom_leiden"]["params"]["resolution"] == resolution
    assert (
        adata_neighbors.uns["custom_leiden"]["params"]["n_iterations"] == n_iterations
    )


@needs.leidenalg
@needs.igraph
@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_leiden_random_state(
    subtests: pytest.Subtests,
    adata_neighbors: AnnData,
    flavor: Literal["igraph", "leidenalg"],
    rng_arg: Literal["rng", "random_state"],
) -> None:
    is_leiden_alg = flavor == "leidenalg"
    n_iterations = 2 if is_leiden_alg else -1
    adata_1, adata_1_again, adata_2 = (
        _leiden_copy(
            adata_neighbors,
            flavor=flavor,
            directed=is_leiden_alg,
            n_iterations=n_iterations,
            **{rng_arg: seed},
        )
        for seed in (1, 1, 42)
    )
    with subtests.test("reproducible"):
        pd.testing.assert_series_equal(
            adata_1.obs["leiden"], adata_1_again.obs["leiden"]
        )
        assert (
            pytest.approx(adata_1.uns["leiden"]["modularity"])
            == adata_1_again.uns["leiden"]["modularity"]
        )
    with subtests.test("different clustering"):
        assert not adata_2.obs["leiden"].equals(adata_1_again.obs["leiden"])
        assert (
            pytest.approx(adata_2.uns["leiden"]["modularity"])
            != adata_1_again.uns["leiden"]["modularity"]
        )


@needs.igraph
def test_leiden_igraph_directed(adata_neighbors):
    with pytest.raises(ValueError, match=r"Cannot use igraph’s leiden.*directed"):
        sc.tl.leiden(adata_neighbors, flavor="igraph", directed=True)


@needs.igraph
def test_leiden_wrong_flavor(adata_neighbors):
    with pytest.raises(ValueError, match=r"flavor must be.*'igraph'.*'leidenalg'.*but"):
        sc.tl.leiden(adata_neighbors, flavor="foo")


@needs.igraph
@needs.leidenalg
def test_leiden_igraph_partition_type(adata_neighbors):
    import leidenalg

    with pytest.raises(ValueError, match=r"Do not pass in partition_type"):
        sc.tl.leiden(
            adata_neighbors,
            flavor="igraph",
            partition_type=leidenalg.RBConfigurationVertexPartition,
        )


@needs.leidenalg
@needs.igraph
def test_leiden_equal_defaults_same_args(adata_neighbors):
    """Ensure the two implementations are the same for the same args."""
    leiden_alg_clustered = _leiden_copy(
        adata_neighbors, flavor="leidenalg", n_iterations=2
    )
    igraph_clustered = _leiden_copy(
        adata_neighbors, flavor="igraph", directed=False, n_iterations=2
    )
    assert (
        normalized_mutual_info_score(
            leiden_alg_clustered.obs["leiden"], igraph_clustered.obs["leiden"]
        )
        > 0.9
    )


@needs.leidenalg
@needs.igraph
def test_leiden_equal_defaults(adata_neighbors):
    """Ensure that the old leidenalg defaults are close enough to the current default outputs."""
    leiden_alg_clustered = _leiden_copy(
        adata_neighbors, flavor="leidenalg", directed=True
    )
    igraph_clustered = _leiden_copy(
        adata_neighbors, flavor="igraph", n_iterations=2, directed=False
    )
    assert (
        normalized_mutual_info_score(
            leiden_alg_clustered.obs["leiden"], igraph_clustered.obs["leiden"]
        )
        > 0.9
    )


@needs.igraph
def test_leiden_objective_function(adata_neighbors):
    """Ensure that popping this as a `clustering_kwargs` and using it does not error out."""
    sc.tl.leiden(
        adata_neighbors,
        objective_function="modularity",
        flavor="igraph",
        directed=False,
    )


@needs.anndata_acc
@needs.igraph
@pytest.mark.parametrize(
    ("clustering", "key"),
    [
        pytest.param(
            partial(sc.tl.leiden, flavor="leidenalg"),
            "leiden",
            marks=needs.leidenalg,
            id="leiden",
        ),
    ],
)
def test_clustering_subset(adata_neighbors, clustering, key):
    from anndata.acc import A

    clustering(adata_neighbors, out=A.obs[key])

    for c in adata_neighbors.obs[key].unique():
        print("Analyzing cluster ", c)
        cells_in_c = adata_neighbors.obs[key] == c
        ncells_in_c = adata_neighbors.obs[key].value_counts().loc[c]
        key_sub = f"{key}_sub"
        clustering(
            adata_neighbors,
            restrict_to=(key, [c]),
            out=A.obs[key_sub],
        )
        # Get new clustering labels
        new_partition = adata_neighbors.obs[key_sub]

        cat_counts = new_partition[cells_in_c].value_counts()

        # Only original cluster's cells assigned to new categories
        assert cat_counts.sum() == ncells_in_c

        # Original category's cells assigned only to new categories
        nonzero_cat = cat_counts[cat_counts > 0].index
        common_cat = nonzero_cat.intersection(adata_neighbors.obs[key].cat.categories)
        assert len(common_cat) == 0


@needs.anndata_acc
@pytest.mark.parametrize(
    ("clustering", "default_key", "default_res", "custom_resolutions"),
    [
        pytest.param(
            partial(sc.tl.leiden, flavor="leidenalg"),
            "leiden",
            0.8,
            [0.9, 1.1],
            marks=needs.leidenalg,
            id="leiden",
        ),
    ],
)
def test_clustering_custom_key(
    adata_neighbors, clustering, default_key, default_res, custom_resolutions
):
    from anndata.acc import A

    custom_keys = [f"{default_key}_{res}" for res in custom_resolutions]

    # Run clustering with default key, then custom keys
    clustering(adata_neighbors, resolution=default_res)
    for key, res in zip(custom_keys, custom_resolutions, strict=True):
        clustering(adata_neighbors, resolution=res, out=A.obs[key])

    # ensure that all clustering parameters are added to user provided keys and not overwritten
    assert adata_neighbors.uns[default_key]["params"]["resolution"] == default_res
    for key, res in zip(custom_keys, custom_resolutions, strict=True):
        assert adata_neighbors.uns[key]["params"]["resolution"] == res


@needs.anndata_acc
@needs.igraph
def test_leiden_out_acc(adata_neighbors: AnnData) -> None:
    """An accessor `out` names the column, and `.uns` follows it."""
    from anndata.acc import A

    sc.tl.leiden(adata_neighbors, flavor="igraph", out=A.obs["mine"])
    assert "mine" in adata_neighbors.obs
    assert "modularity" in adata_neighbors.uns["mine"]
    assert "leiden" not in adata_neighbors.uns
    # an explicit accessor isn’t renamed when restricting
    sc.tl.leiden(
        adata_neighbors, flavor="igraph", restrict_to=("mine", ["0"]), out="obs.leiden"
    )
    assert "leiden" in adata_neighbors.obs
    assert "leiden_R" not in adata_neighbors.obs
    # neither is a prefixed one
    sc.tl.leiden(
        adata_neighbors,
        flavor="igraph",
        restrict_to=("mine", ["0"]),
        out=sc.Prefix("p"),
    )
    assert "p_leiden" in adata_neighbors.obs
    assert "p_leiden" in adata_neighbors.uns


@needs.igraph
def test_leiden_out_true(adata_neighbors: AnnData) -> None:
    """`out=True` returns the labels and writes nothing."""
    obs_before, uns_before = adata_neighbors.obs.copy(), set(adata_neighbors.uns)
    labels = sc.tl.leiden(adata_neighbors, flavor="igraph", out=True)
    assert isinstance(labels, pd.Categorical)
    assert len(labels) == adata_neighbors.n_obs
    pd.testing.assert_frame_equal(adata_neighbors.obs, obs_before)
    assert set(adata_neighbors.uns) == uns_before


@needs.anndata_acc
@needs.igraph
@pytest.mark.parametrize("restrict", [False, True], ids=["plain", "restrict"])
def test_leiden_key_added_deprecated(
    adata_neighbors: AnnData, *, restrict: bool
) -> None:
    """`key_added='k'` is `out=A.obs['k']`, except for the `_R` suffix for the default."""
    from anndata.acc import A

    kw = dict(flavor="igraph", restrict_to=("louvain", ["0"]) if restrict else None)
    key = "leiden_R" if restrict else "leiden"
    expected = adata_neighbors.copy()
    sc.tl.leiden(expected, out=A.obs[key], **kw)
    with pytest.warns(FutureWarning, match=r"argument key_added is deprecated"):
        sc.tl.leiden(adata_neighbors, key_added="leiden", **kw)
    assert key in adata_neighbors.obs
    assert adata_neighbors.obs.keys().equals(expected.obs.keys())
    assert adata_neighbors.uns.keys() == expected.uns.keys()

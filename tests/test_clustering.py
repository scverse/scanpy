from __future__ import annotations

import pandas as pd
import pytest
from sklearn.metrics.cluster import normalized_mutual_info_score

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


@pytest.fixture
def adata_neighbors():
    return pbmc68k_reduced()


FLAVORS = [
    pytest.param("igraph", marks=needs.igraph),
    pytest.param("leidenalg", marks=needs.leidenalg),
]


@needs.leidenalg
@needs.igraph
@pytest.mark.parametrize("flavor", FLAVORS)
@pytest.mark.parametrize("resolution", [1, 2])
@pytest.mark.parametrize("n_iterations", [-1, 3])
def test_leiden_basic(adata_neighbors, flavor, resolution, n_iterations):
    sc.tl.leiden(
        adata_neighbors,
        flavor=flavor,
        resolution=resolution,
        n_iterations=n_iterations,
        directed=(flavor == "leidenalg"),
        key_added="leiden_custom",
    )
    assert adata_neighbors.uns["leiden_custom"]["params"]["resolution"] == resolution
    assert (
        adata_neighbors.uns["leiden_custom"]["params"]["n_iterations"] == n_iterations
    )


@needs.leidenalg
@needs.igraph
@pytest.mark.parametrize("flavor", FLAVORS)
def test_leiden_random_state(adata_neighbors, flavor):
    is_leiden_alg = flavor == "leidenalg"
    n_iterations = 2 if is_leiden_alg else -1
    adata_1 = sc.tl.leiden(
        adata_neighbors,
        flavor=flavor,
        random_state=1,
        copy=True,
        directed=is_leiden_alg,
        n_iterations=n_iterations,
    )
    adata_1_again = sc.tl.leiden(
        adata_neighbors,
        flavor=flavor,
        random_state=1,
        copy=True,
        directed=is_leiden_alg,
        n_iterations=n_iterations,
    )
    adata_2 = sc.tl.leiden(
        adata_neighbors,
        flavor=flavor,
        random_state=2,
        copy=True,
        directed=is_leiden_alg,
        n_iterations=n_iterations,
    )
    pd.testing.assert_series_equal(adata_1.obs["leiden"], adata_1_again.obs["leiden"])
    assert not adata_2.obs["leiden"].equals(adata_1_again.obs["leiden"])


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
    leiden_alg_clustered = sc.tl.leiden(
        adata_neighbors, flavor="leidenalg", copy=True, n_iterations=2
    )
    igraph_clustered = sc.tl.leiden(
        adata_neighbors, flavor="igraph", copy=True, directed=False, n_iterations=2
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
    leiden_alg_clustered = sc.tl.leiden(
        adata_neighbors, flavor="leidenalg", directed=True, copy=True
    )
    igraph_clustered = sc.tl.leiden(
        adata_neighbors, copy=True, n_iterations=2, directed=False
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


@needs.igraph
@pytest.mark.parametrize(
    ("clustering", "key"),
    [
        pytest.param(sc.tl.louvain, "louvain", marks=needs.louvain),
        pytest.param(sc.tl.leiden, "leiden", marks=needs.leidenalg),
    ],
)
def test_clustering_subset(adata_neighbors, clustering, key):
    clustering(adata_neighbors, key_added=key)

    for c in adata_neighbors.obs[key].unique():
        print("Analyzing cluster ", c)
        cells_in_c = adata_neighbors.obs[key] == c
        ncells_in_c = adata_neighbors.obs[key].value_counts().loc[c]
        key_sub = str(key) + "_sub"
        clustering(
            adata_neighbors,
            restrict_to=(key, [c]),
            key_added=key_sub,
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


@needs.louvain
@needs.igraph
def test_louvain_basic(adata_neighbors):
    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")


@needs.louvain
@pytest.mark.parametrize("random_state", [10, 999])
@pytest.mark.parametrize("resolution", [0.9, 1.1])
def test_louvain_custom_key(adata_neighbors, resolution, random_state):
    sc.tl.louvain(
        adata_neighbors,
        key_added="louvain_custom",
        random_state=random_state,
        resolution=resolution,
    )
    assert (
        adata_neighbors.uns["louvain_custom"]["params"]["random_state"] == random_state
    )
    assert adata_neighbors.uns["louvain_custom"]["params"]["resolution"] == resolution


@needs.louvain
@needs.igraph
def test_partition_type(adata_neighbors):
    import louvain

    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)


@pytest.mark.parametrize(
    ("clustering", "default_key", "default_res", "custom_resolutions"),
    [
        pytest.param(sc.tl.leiden, "leiden", 0.8, [0.9, 1.1], marks=needs.leidenalg),
        pytest.param(sc.tl.louvain, "louvain", 0.8, [0.9, 1.1], marks=needs.louvain),
    ],
)
def test_clustering_custom_key(
    adata_neighbors, clustering, default_key, default_res, custom_resolutions
):
    custom_keys = [f"{default_key}_{res}" for res in custom_resolutions]

    # Run clustering with default key, then custom keys
    clustering(adata_neighbors, resolution=default_res)
    for key, res in zip(custom_keys, custom_resolutions, strict=True):
        clustering(adata_neighbors, resolution=res, key_added=key)

    # ensure that all clustering parameters are added to user provided keys and not overwritten
    assert adata_neighbors.uns[default_key]["params"]["resolution"] == default_res
    for key, res in zip(custom_keys, custom_resolutions, strict=True):
        assert adata_neighbors.uns[key]["params"]["resolution"] == res

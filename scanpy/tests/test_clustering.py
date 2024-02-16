from __future__ import annotations

import pytest
from sklearn.metrics.cluster import normalized_mutual_info_score

import scanpy as sc
from scanpy.testing._helpers.data import pbmc68k_reduced
from scanpy.testing._pytest.marks import needs


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
    )
    assert adata_neighbors.uns["leiden"]["params"]["resolution"] == resolution
    assert adata_neighbors.uns["leiden"]["params"]["n_iterations"] == n_iterations


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
    assert (adata_1.obs["leiden"] == adata_1_again.obs["leiden"]).all()
    # This random state produces different categories so can't check the arrays against each other.
    assert (adata_2.obs["leiden"] != adata_1_again.obs["leiden"]).any()


@needs.igraph
def test_leiden_igraph_directed(adata_neighbors):
    with pytest.raises(ValueError):
        sc.tl.leiden(adata_neighbors, flavor="igraph")


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
    "clustering,key",
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
@needs.igraph
def test_partition_type(adata_neighbors):
    import louvain

    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)

from __future__ import annotations

import pytest

import scanpy as sc
from scanpy.testing._helpers.data import pbmc68k_reduced
from scanpy.testing._pytest.marks import needs
from sklearn.metrics.cluster import normalized_mutual_info_score

@pytest.fixture
def adata_neighbors():
    return pbmc68k_reduced()


@needs.leidenalg
@pytest.mark.parametrize("use_igraph,directed", [pytest.param(True, False), pytest.param(False, True), pytest.param(False, False)])
@pytest.mark.parametrize("use_weights", [True, False])
@pytest.mark.parametrize("resolution", [1, 2])
@pytest.mark.parametrize("n_iterations", [-1, 3])
def test_leiden_basic(adata_neighbors, use_igraph, directed, use_weights, resolution, n_iterations):
    sc.tl.leiden(adata_neighbors, use_igraph=use_igraph, use_weights=use_weights, resolution=resolution, n_iterations=n_iterations, directed=directed)
    assert adata_neighbors.uns['leiden']['params']['resolution'] == resolution
    assert adata_neighbors.uns['leiden']['params']['n_iterations'] == n_iterations

def test_leiden_igraph_directed(adata_neighbors):
    with pytest.raises(ValueError):
        sc.tl.leiden(adata_neighbors, directed=True)

@needs.leidenalg
def test_leiden_equal_defaults(adata_neighbors):
    leiden_alg_clustered = sc.tl.leiden(adata_neighbors, use_igraph=False, copy=True)
    igraph_clustered = sc.tl.leiden(adata_neighbors, use_igraph=True, copy=True)
    assert normalized_mutual_info_score(leiden_alg_clustered.obs['leiden'],igraph_clustered.obs['leiden']) > .9

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
def test_louvain_basic(adata_neighbors):
    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")


@needs.louvain
def test_partition_type(adata_neighbors):
    import louvain

    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)

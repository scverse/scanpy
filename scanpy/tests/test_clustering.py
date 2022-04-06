from importlib.util import find_spec

import pytest
import scanpy as sc
from scanpy.tests._data._cached_datasets import pbmc68k_reduced


@pytest.fixture
def adata_neighbors():
    return pbmc68k_reduced()


needs_louvain = pytest.mark.skipif(
    not find_spec("louvain"), reason="needs module `louvain`"
)
needs_leiden = pytest.mark.skipif(
    not find_spec("leidenalg"), reason="needs module `leidenalg`"
)


@needs_leiden
def test_leiden_basic(adata_neighbors):
    sc.tl.leiden(adata_neighbors)


@pytest.mark.parametrize(
    'clustering,key',
    [
        pytest.param(sc.tl.louvain, 'louvain', marks=needs_louvain),
        pytest.param(sc.tl.leiden, 'leiden', marks=needs_leiden),
    ],
)
def test_clustering_subset(adata_neighbors, clustering, key):
    clustering(adata_neighbors, key_added=key)

    for c in adata_neighbors.obs[key].unique():
        print('Analyzing cluster ', c)
        cells_in_c = adata_neighbors.obs[key] == c
        ncells_in_c = adata_neighbors.obs[key].value_counts().loc[c]
        key_sub = str(key) + '_sub'
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


def test_louvain_basic(adata_neighbors):
    pytest.importorskip("louvain")

    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")


def test_partition_type(adata_neighbors):
    louvain = pytest.importorskip("louvain")

    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)

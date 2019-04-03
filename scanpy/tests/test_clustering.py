import pytest
import scanpy as sc


@pytest.fixture
def adata_neighbors():
    return sc.datasets.pbmc68k_reduced()


def test_leiden_basic(adata_neighbors):
    sc.tl.leiden(adata_neighbors)


def test_louvain_basic(adata_neighbors):
    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")


def test_partition_type(adata_neighbors):
    import louvain
    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)

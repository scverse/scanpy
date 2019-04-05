import numpy as np
import pytest
import scanpy as sc


@pytest.fixture
def adata_neighbors():
    adata = sc.read(
        './data/pbmc3k_raw.h5ad',
        backup_url='http://falexwolf.de/data/pbmc3k_raw.h5ad',
    )
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    return adata


def test_leiden_basic(adata_neighbors):
    sc.tl.leiden(adata_neighbors)

@pytest.mark.parametrize('clustering,key', [
    (sc.tl.louvain, 'louvain'),
    (sc.tl.leiden, 'leiden'),
])
def test_clustering_subset(adata_neighbors, clustering, key):
    clustering(adata_neighbors, key_added=key)
    print(adata_neighbors)
    print(adata_neighbors.obs[key].value_counts())

    for c in np.unique(adata_neighbors.obs[key]):
        print('Analyzing cluster ', c)
        cells_in_c = adata_neighbors.obs[key] == c
        ncells_in_c = adata_neighbors.obs[key].value_counts().loc[c]
        key_sub = str(key) + '_sub'
        clustering(
            adata_neighbors,
            restrict_to=(key, [c]),
            key_added=key_sub
        ) 
        new_partition = adata_neighbors.obs[key_sub]

        # Old category should not be present
        categories = new_partition.value_counts().index
        assert c not in categories

        # Cells in the original category are assigned only to new categories
        print('Categories')
        cat_counts = new_partition[cells_in_c].value_counts()
        print(cat_counts)

        ## Only cells of the original cluster have been assigned to
        ## new categories
        assert cat_counts.sum() == ncells_in_c

        print('Non-zero categories')
        nonzero_cat = cat_counts[cat_counts > 0].index
        print(nonzero_cat)

        common_cat = nonzero_cat & adata_neighbors.obs[key].cat.categories
        print('Common categories')
        print(common_cat)
        ## Only new categories
        assert len(common_cat) == 0


def test_louvain_basic(adata_neighbors):
    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")


def test_partition_type(adata_neighbors):
    import louvain
    sc.tl.louvain(adata_neighbors, partition_type=louvain.RBERVertexPartition)
    sc.tl.louvain(adata_neighbors, partition_type=louvain.SurpriseVertexPartition)

import pytest
import scanpy.api as sc

@pytest.fixture
def adata_neighbors():
    adata = sc.read('./data/pbmc3k_raw.h5ad',
                    backup_url='http://falexwolf.de/data/pbmc3k_raw.h5ad')
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    return adata

def test_nofail(adata_neighbors): 
    sc.tl.louvain(adata_neighbors)
    sc.tl.louvain(adata_neighbors, use_weights=True)
    sc.tl.louvain(adata_neighbors, use_weights=True, flavor="igraph")
    sc.tl.louvain(adata_neighbors, flavor="igraph")
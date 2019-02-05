import scanpy as sc

myColors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
            '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe']

def pre_preprocessing(path):
    adata = sc.read(path)
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_per_cell(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    return adata

def clustering(adata, resolution=0.8):
    sc.tl.leiden(adata, resolution=resolution)

def test_tSNE(adata, palette, random_state=10, n_components=3):
    sc.tl.tsne(adata, random_state=random_state, n_components=n_components)
    sc.pl.tsne(adata, color='leiden', components=['1,2', '2,3'], palette=palette)
    sc.pl.tsne(adata, color='leiden', projection='3d', palette=palette)

if __name__ == "__main__":

    path = '../datasets/10x_pbmc68k_reduced.h5ad'
    test_data = pre_preprocessing(path)
    clustering(test_data)
    test_tSNE(test_data, myColors)

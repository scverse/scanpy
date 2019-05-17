# -*- coding: utf-8 -*-
import pandas as pd
import scanpy as sc

# Scanpy
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

#####
# We have already performed weighted sampling and PCA on weighted sampling.
# We will red the following files

# 1. Sampled matrix to use all genes when we find marker genes, Rows are data points and columns are features
# 2. Weighted matrix , where each entry is a weight of rows/data point in Sampled matrix
# 3. Sampled PCA matrix, we also performed PCA in separately because according to my knowledge scanpy do not
#    support sparse PCA for weighted samples. We use this for clustering.

# 1.
fileName =  '_data/C1000.txt'
raw_data = sc.read_text(fileName, delimiter = ",")

# 2.
weightName = '_data/W1000.txt'
weights = pd.read_csv(weightName, header=None, sep=',')

# 3.
adata = sc.read_text('_data/C1000_WPCA10.txt', delimiter = ",")
adata.raw = raw_data

sc.pp.neighbors(adata, n_neighbors=10, use_rep='X') # use_rep='X' , n_pcs=10
sc.tl.umap(adata, min_dist=0.5)

# Louvain Clustering
#sc.tl.louvain(adata, resolution=1.0)
#sc.pl.umap(adata, color=['louvain'], s = 30)
sc.tl.leiden(adata, resolution=1.0)
sc.pl.umap(adata, color=['leiden'], s = 30)

marker_genes = ['RPS12',  'LTB', 'IL32', 'LDHB', 'CD74', 'LYZ', 'FTL', 'CCL5', 'NKG7', 'GNLY',
                'CST3', 'FCER1G', 'AIF1', 'CD79B']

def test_genes_ranking():

    # 'weights' is None by default
    # Find Marker Genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', weights=weights)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', weights=weights)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

def test_heatmap():
    sc.pl.heatmap(adata, marker_genes, groupby='leiden', figsize=(5, 8), use_raw=True, vmin=-3, vmax=3, cmap='bwr',
                  var_group_rotation=0, dendrogram=True, weights=weights)


def test_dotplot():
    sc.pl.dotplot(adata, marker_genes, groupby='leiden', use_raw=True, dendrogram=True, weights=weights)


def test_stacked_violin():
    sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', dendrogram=True, weights=weights)

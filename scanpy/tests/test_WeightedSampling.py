# -*- coding: utf-8 -*-
import pandas as pd
import scanpy as sc

from pathlib import Path
import os.path

from matplotlib.testing import setup
import matplotlib.pyplot as pl
from matplotlib.testing.compare import compare_images
setup()

HERE: Path = Path(__file__).parent
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/_images/'
DATA = HERE / '_data/weighted_sampled'

sc.pl.set_rcParams_defaults()
sc.set_figure_params(dpi=40, color_map='viridis')

# Scanpy
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

#####
# Test images are saved under the folder ./figures
# if test images need to be updated, simply copy them from
# the ./figures folder to ./_images/


#####
# We have already performed weighted sampling and PCA on weighted sampling.
# We will red the following files

# 1. Sampled matrix to use all genes when we find marker genes, Rows are data points and columns are features
# 2. Weighted matrix , where each entry is a weight of rows/data point in Sampled matrix
# 3. Sampled PCA matrix, we also performed PCA in separately because according to my knowledge scanpy do not
#    support sparse PCA for weighted samples. We use this for clustering.

# 1.
fileName = DATA / 'C1000.txt'
raw_data = sc.read_text(fileName, delimiter = ",")

# 2.
weightName = DATA / 'W1000.txt'
weights = pd.read_csv(weightName, header=None, sep=',')

# 3.
pcaName = DATA / 'C1000_WPCA10.txt'
adata = sc.read_text(pcaName, delimiter = ",")
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

def save_and_compare_images(basename, tolerance=20):
    if not os.path.exists('./figures/'): os.makedirs('./figures/')
    outname = './figures/' + basename + '.png'
    pl.savefig(outname, dpi=40)
    pl.close()
    res = compare_images(ROOT + '/' + basename + '.png', outname, tolerance)
    assert res is None, res

def test_weighted_genes_ranking_t_test():

    # 'weights' is None by default
    # Find Marker Genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', weights=weights)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    save_and_compare_images('weighted_genes_ranking_t_test')
    #assert np.all(adata.uns['t-test']['names'][:5] == ['RPS12', 'LTB', 'IL32', 'LDHB', 'CD74'])

def test_weighted_genes_ranking_wilcoxon():
    # 'weights' is None by default
    # Find Marker Genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', weights=weights)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    save_and_compare_images('weighted_genes_ranking_wilcoxon')
    # assert np.all(adata.uns['wilcoxon']['names'][:5] == ['RPS12', 'LTB', 'IL32', 'LDHB', 'CD74'])

def test_weighted_dotplot():
    sc.pl.dotplot(adata, marker_genes, groupby='leiden', use_raw=True, dendrogram=True, weights=weights)
    save_and_compare_images('weighted_dotplot')

def test_weighted_heatmap():
    sc.pl.heatmap(adata, marker_genes, groupby='leiden', figsize=(5, 8), use_raw=True, vmin=-3, vmax=3, cmap='bwr',
                  var_group_rotation=0, dendrogram=True, weights=weights)
    save_and_compare_images('weighted_heatmap')

def test_weighted_stacked_violin():
    sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', dendrogram=True, weights=weights)
    save_and_compare_images('weighted_stacked_violin')

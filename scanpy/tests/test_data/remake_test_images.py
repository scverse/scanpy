# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl

import scanpy.api as sc

#sc.set_figure_params()

adata = sc.datasets.krumsiek11()

# make heatmap
sc.pl.heatmap(adata, adata.var_names, 'cell_type', use_raw=False)
pl.savefig('master_heatmap.png', dpi=80)
pl.close()

# make heatmap with continues data
adata.obs['Gata2'] = adata.X[:, 0]
sc.pl.heatmap(adata, adata.var_names, 'Gata2', use_raw=False,
              num_categories=4, figsize=(4.5, 5))
pl.savefig('master_heatmap2.png', dpi=80)
pl.close()

# make dotplot
sc.pl.dotplot(adata, adata.var_names, 'cell_type', use_raw=False)
pl.savefig('master_dotplot.png', dpi=80)
pl.close()

# make dotplot with continuous data
sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
              num_categories=7, figsize=(7, 2.5))
pl.savefig('master_dotplot2.png', dpi=80)
pl.close()

sc.pl.matrixplot(adata, adata.var_names, 'cell_type', use_raw=False)
pl.savefig('master_matrixplot.png', dpi=80)
pl.close()

sc.pl.matrixplot(adata, adata.var_names, 'Gata2', use_raw=False,
                 num_categories=4, figsize=(8, 2.5), cmap='viridis')
pl.savefig('master_matrixplot2.png', dpi=80)
pl.close()

# make stacked violing plot
sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, color='blue')
pl.title("image may have cut labels.\nThis is ok for test")
pl.savefig('master_stacked_violin.png', dpi=80)
pl.close()

# make stacked violing plot with swapped axes
sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, swap_axes=True, figsize=(3, 5))
pl.savefig('master_stacked_violin_swapped_axes.png', dpi=80)
pl.close()


####
# tests based on pbmc68k_reduced dataset

pbmc = sc.datasets.pbmc68k_reduced()
sc.pl.violin(pbmc, ['n_genes', 'percent_mito', 'n_counts'], stripplot=False, multi_panel=True, jitter=False)
pl.savefig("master_violin_multi_panel.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3)
pl.savefig("master_ranked_genes_sharey.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3, sharey=False)
pl.savefig("master_ranked_genes.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=5)
pl.savefig("master_ranked_genes_heatmap.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups_stacked_violin(pbmc, n_genes=3)
pl.savefig("master_ranked_genes_stacked_violin.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4)
pl.savefig("master_ranked_genes_dotplot.png", dpi=80)
pl.close()

sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=5)
pl.savefig('master_ranked_genes_matrixplot.png', dpi=80)
pl.close()

sc.pl.rank_genes_groups_violin(pbmc, groups=pbmc.obs.bulk_labels.cat.categories[0], n_genes=5,
                               jitter=False, strip=False)
pl.savefig("master_ranked_genes_violin.png", dpi=80)
pl.close()

sc.pl.umap(pbmc, color='louvain')
pl.savefig("master_umap.png", dpi=80)
pl.close()

sc.pl.umap(pbmc, color=['LYZ', 'CD79A'], size=20, alpha=0.5)
pl.savefig("master_umap_gene_expr.png", dpi=80)
pl.close()

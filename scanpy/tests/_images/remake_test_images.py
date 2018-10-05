# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl

import scanpy.api as sc

sc.set_figure_params(dpi=80, color_map='viridis')


def make_heatmaps():
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


def make_dotplots():
    # make dotplot
    adata = sc.datasets.krumsiek11()
    sc.pl.dotplot(adata, adata.var_names, 'cell_type', use_raw=False)
    pl.savefig('master_dotplot.png', dpi=80)
    pl.close()

    # make dotplot with continuous data
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=7, figsize=(7, 2.5))
    pl.savefig('master_dotplot2.png', dpi=80)
    pl.close()


def make_matrix_plots():
    adata = sc.datasets.krumsiek11()
    sc.pl.matrixplot(adata, adata.var_names, 'cell_type', use_raw=False)
    pl.savefig('master_matrixplot.png', dpi=80)
    pl.close()

    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.matrixplot(adata, adata.var_names, 'Gata2', use_raw=False,
                     num_categories=4, figsize=(8, 2.5), cmap='viridis')
    pl.savefig('master_matrixplot2.png', dpi=80)
    pl.close()


def make_violin_plots():
    # make stacked violin plot
    adata = sc.datasets.krumsiek11()
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, color='blue')
    pl.title("image may have cut labels.\nThis is ok for test")
    pl.savefig('master_stacked_violin.png', dpi=80)
    pl.close()

    # make stacked violing plot with swapped axes
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, swap_axes=True, figsize=(3, 5))
    pl.savefig('master_stacked_violin_swapped_axes.png', dpi=80)
    pl.close()


def make_violin_images_multi():
    ####
    # tests based on pbmc68k_reduced dataset
    sc.pl.set_rcParams_defaults()
    sc.set_figure_params(dpi=80, color_map='viridis')

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pl.violin(pbmc, ['n_genes', 'percent_mito', 'n_counts'],
                 stripplot=False, multi_panel=True, jitter=False)
    pl.savefig("master_violin_multi_panel.png", dpi=80)
    pl.close()


def make_rank_genes_groups_plots():
    pbmc = sc.datasets.pbmc68k_reduced()
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


def make_scatterplots():
    pbmc = sc.datasets.pbmc68k_reduced()

    # pca
    sc.pl.pca(pbmc, color='bulk_labels')
    pl.savefig("master_pca.png", dpi=80)
    pl.close()

    # 3d
    sc.pl.pca(pbmc, color='bulk_labels', projection='3d')
    pl.savefig("master_3dprojection.png", dpi=80)
    pl.close()

    # multipanel
    sc.pl.pca(pbmc, color=['CD3D', 'CD79A'], components=['1,2', '1,3'], vmax=5, use_raw=False, vmin=-5, cmap='seismic')
    pl.savefig("master_multipanel.png", dpi=80)
    pl.close()

    # tsne
    sc.tl.tsne(pbmc, random_state=0, n_pcs=30)
    sc.pl.tsne(pbmc, color=['CD3D', 'louvain'])
    pl.savefig("master_tsne.png", dpi=80)
    pl.close()

    sc.pl.umap(pbmc, color=['louvain'],
               palette=['b', 'g', 'r', 'yellow', 'black', 'gray', 'lightblue'], frameon=False)
    pl.savefig("master_umap.png", dpi=80)
    pl.close()

    sc.pl.umap(pbmc, color=['LYZ', 'CD79A'], s=20, alpha=0.5, frameon=False)
    pl.savefig("master_umap_gene_expr.png", dpi=80)
    pl.close()

    #  edges = True
    sc.pp.neighbors(pbmc)
    sc.pl.umap(pbmc, color='louvain', edges=True, edges_width=0.1, s=50)
    pl.savefig("master_umap_with_edges.png", dpi=80)
    pl.close()

    # diffmap
    sc.tl.diffmap(pbmc)
    sc.pl.diffmap(pbmc, components='all', color=['CD3D'])
    pl.savefig("master_diffmap.png", dpi=80)
    pl.close()


make_heatmaps()
make_dotplots()
make_matrix_plots()
make_violin_plots()
make_violin_images_multi()
make_rank_genes_groups_plots()
make_scatterplots()


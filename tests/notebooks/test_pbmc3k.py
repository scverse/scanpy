# *First compiled on May 5, 2017. Updated August 14, 2018.*
# # Clustering 3k PBMCs following a Seurat Tutorial
#
# This started out with a demonstration that Scanpy would allow to reproduce most of Seurat's
# ([Satija *et al.*, 2015](https://doi.org/10.1038/nbt.3192)) clustering tutorial as described on
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html (July 26, 2017), which we gratefully acknowledge.
# In the meanwhile, we have added and removed several pieces.
#
# The data consists in *3k PBMCs from a Healthy Donor* and is freely available from 10x Genomics
# ([here](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
# from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)).
from __future__ import annotations

from functools import partial
from pathlib import Path

import numpy as np
from matplotlib.testing import setup

setup()

import scanpy as sc
from testing.scanpy._pytest.marks import needs

HERE: Path = Path(__file__).parent
ROOT = HERE / "_images_pbmc3k"


@needs.leidenalg
def test_pbmc3k(image_comparer):  # noqa: PLR0915
    # ensure violin plots and other non-determinstic plots have deterministic behavior
    np.random.seed(0)
    save_and_compare_images = partial(image_comparer, ROOT, tol=20)
    adata = sc.datasets.pbmc3k()

    # Preprocessing

    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    save_and_compare_images("highest_expr_genes")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    mito_genes = [name for name in adata.var_names if name.startswith("MT-")]
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse to transform to a dense array after summing
    adata.obs["percent_mito"] = (
        np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    )
    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1

    sc.pl.violin(
        adata,
        ["n_genes", "n_counts", "percent_mito"],
        jitter=False,
        multi_panel=True,
        show=False,
    )
    save_and_compare_images("violin")

    sc.pl.scatter(adata, x="n_counts", y="percent_mito", show=False)
    save_and_compare_images("scatter_1")
    sc.pl.scatter(adata, x="n_counts", y="n_genes", show=False)
    save_and_compare_images("scatter_2")

    adata = adata[adata.obs["n_genes"] < 2500, :]
    adata = adata[adata.obs["percent_mito"] < 0.05, :]

    adata.raw = sc.pp.log1p(adata, copy=True)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    filter_result = sc.pp.filter_genes_dispersion(
        adata.X,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
    )
    sc.pl.filter_genes_dispersion(filter_result, show=False)
    save_and_compare_images("filter_genes_dispersion")

    adata = adata[:, filter_result.gene_subset]
    sc.pp.log1p(adata)
    sc.pp.regress_out(adata, ["n_counts", "percent_mito"])
    sc.pp.scale(adata, max_value=10)

    # PCA

    sc.pp.pca(adata, svd_solver="arpack")
    sc.pl.pca(adata, color="CST3", show=False)
    save_and_compare_images("pca")

    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    save_and_compare_images("pca_variance_ratio")

    # UMAP

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    # sc.tl.umap(adata)  # umaps lead to slight variations

    # sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False, show=False)
    # save_and_compare_images('umap_1')

    # Clustering the graph

    sc.tl.leiden(
        adata,
        resolution=0.9,
        random_state=1,
        directed=False,
        n_iterations=2,
        flavor="igraph",
    )

    # sc.pl.umap(adata, color=["leiden", "CST3", "NKG7"], show=False)
    # save_and_compare_images("umap_2")
    sc.pl.scatter(adata, "CST3", "NKG7", color="leiden", show=False)
    save_and_compare_images("scatter_3")

    # Finding marker genes
    # Due to incosistency with our test runner vs local, these clusters need to
    # be pre-annotated as the numbers for each cluster are not consistent.
    marker_genes = [
        *["RP11-18H21.1", "GZMK", "CD79A", "FCGR3A"],
        *["GNLY", "S100A8", "FCER1A", "PPBP"],
    ]
    data_df = adata[:, marker_genes].to_df()
    data_df["leiden"] = adata.obs["leiden"]
    max_idxs = data_df.groupby("leiden", observed=True).mean().idxmax()
    assert not max_idxs[marker_genes][
        max_idxs[marker_genes].duplicated(keep=False)
    ].tolist(), "Not all marker genes are unique per cluster"
    leiden_relabel = {
        max_idxs[marker_gene]: str(i) for i, marker_gene in enumerate(marker_genes)
    }
    adata.obs["leiden_old"] = adata.obs["leiden"].copy()
    adata.rename_categories(
        "leiden", [leiden_relabel[key] for key in sorted(leiden_relabel.keys())]
    )
    # ensure that the column can be sorted for consistent plotting since it is by default unordered
    adata.obs["leiden"] = adata.obs["leiden"].cat.reorder_categories(
        list(map(str, range(len(adata.obs["leiden"].cat.categories)))), ordered=True
    )

    sc.tl.rank_genes_groups(adata, "leiden")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    save_and_compare_images("rank_genes_groups_1")

    sc.tl.rank_genes_groups(adata, "leiden", method="logreg")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    save_and_compare_images("rank_genes_groups_2")

    sc.tl.rank_genes_groups(adata, "leiden", groups=["0"], reference="1")
    sc.pl.rank_genes_groups(adata, groups="0", n_genes=20, show=False)
    save_and_compare_images("rank_genes_groups_3")

    # gives a strange error, probably due to jitter or something
    # sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
    # save_and_compare_images('rank_genes_groups_4')

    new_cluster_names = [
        *["CD4 T cells", "CD8 T cells", "B cells", "NK cells"],
        *["FCGR3A+ Monocytes", "CD14+ Monocytes", "Dendritic cells", "Megakaryocytes"],
    ]
    adata.rename_categories("leiden", new_cluster_names)

    # sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, show=False)
    # save_and_compare_images('umap_3')
    sc.pl.violin(
        adata, ["CST3", "NKG7", "PPBP"], groupby="leiden", rotation=90, show=False
    )
    save_and_compare_images("violin_2")

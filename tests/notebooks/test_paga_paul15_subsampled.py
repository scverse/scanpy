# PAGA for hematopoiesis in mouse [(Paul *et al.*, 2015)](https://doi.org/10.1016/j.cell.2015.11.013)
# Hematopoiesis: trace myeloid and erythroid differentiation for data of [Paul *et al.* (2015)](https://doi.org/10.1016/j.cell.2015.11.013).
#
# This is the subsampled notebook for testing.
from __future__ import annotations

from functools import partial
from pathlib import Path

import numpy as np
import pytest
from matplotlib.testing import setup

import scanpy as sc
from testing.scanpy._helpers.data import paul15
from testing.scanpy._pytest.marks import needs

HERE: Path = Path(__file__).parent
ROOT = HERE / "_images_paga_paul15_subsampled"


@pytest.mark.skip(reason="Broken, needs fixing")
@needs.igraph
@needs.louvain
def test_paga_paul15_subsampled(image_comparer, plt):
    setup()
    save_and_compare_images = partial(image_comparer, ROOT, tol=25)

    adata = paul15()
    sc.pp.subsample(adata, n_obs=200)
    del adata.uns["iroot"]
    adata.X = adata.X.astype("float64")

    # Preprocessing and Visualization
    sc.pp.recipe_zheng17(adata)
    sc.pp.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")

    sc.tl.diffmap(adata)
    sc.tl.diffmap(adata)  # See #1262
    sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_diffmap")
    sc.tl.draw_graph(adata)

    sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")

    # TODO: currently needs skip if louvain isn't installed, do major rework

    # Clustering and PAGA
    sc.tl.louvain(adata, resolution=1.0)
    sc.tl.paga(adata, groups="louvain")
    # sc.pl.paga(adata, color=['louvain', 'Hba-a2', 'Elane', 'Irf8'])
    # sc.pl.paga(adata, color=['louvain', 'Itga2b', 'Prss34'])

    adata.obs["louvain_anno"] = adata.obs["louvain"]
    sc.tl.paga(adata, groups="louvain_anno")

    PAGA_CONNECTIVITIES = np.array(
        [
            [0.0, 0.128553, 0.0, 0.07825, 0.0, 0.0, 0.238741, 0.0, 0.0, 0.657049],
            [
                *[0.128553, 0.0, 0.480676, 0.257505, 0.533036],
                *[0.043871, 0.0, 0.032903, 0.0, 0.087743],
            ],
        ]
    )

    assert np.allclose(
        adata.uns["paga"]["connectivities"].toarray()[:2],
        PAGA_CONNECTIVITIES,
        atol=1e-4,
    )

    sc.pl.paga(adata, threshold=0.03)

    # !!!! no clue why it doesn't produce images with the same shape
    # save_and_compare_images('paga')

    sc.tl.draw_graph(adata, init_pos="paga")
    sc.pl.paga_compare(
        adata,
        threshold=0.03,
        title="",
        right_margin=0.2,
        size=10,
        edge_width_scale=0.5,
        legend_fontsize=12,
        fontsize=12,
        frameon=False,
        edges=True,
    )

    # slight deviations because of graph drawing
    # save_and_compare_images('paga_compare')

    adata.uns["iroot"] = np.flatnonzero(adata.obs["louvain_anno"] == "3")[0]
    sc.tl.dpt(adata)
    gene_names = [
        "Gata2",
        "Gata1",
        "Klf1",
        "Hba-a2",  # erythroid
        "Elane",
        "Cebpe",  # neutrophil
        "Irf8",
    ]  # monocyte

    paths = [
        ("erythrocytes", [3, 9, 0, 6]),
        ("neutrophils", [3, 1, 2]),
        ("monocytes", [3, 1, 4, 5]),
    ]

    adata.obs["distance"] = adata.obs["dpt_pseudotime"]

    _, axs = plt.subplots(
        ncols=3, figsize=(6, 2.5), gridspec_kw={"wspace": 0.05, "left": 0.12}
    )
    plt.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
    for ipath, (descr, path) in enumerate(paths):
        _, data = sc.pl.paga_path(
            adata,
            path,
            gene_names,
            show_node_names=False,
            ax=axs[ipath],
            ytick_fontsize=12,
            left_margin=0.15,
            n_avg=50,
            annotations=["distance"],
            show_yticks=ipath == 0,
            show_colorbar=False,
            color_map="Greys",
            color_maps_annotations={"distance": "viridis"},
            title=f"{descr} path",
            return_data=True,
            show=False,
        )
        # add a test for this at some point
        # data.to_csv(f"./write/paga_path_{descr}.csv")

    save_and_compare_images("paga_path")

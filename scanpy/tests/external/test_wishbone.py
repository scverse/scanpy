import pytest

import scanpy as sc
import scanpy.external as sce

pytest.importorskip("wishbone")


def test_run_wishbone():
    adata = sc.datasets.pbmc3k()
    sc.preprocessing.normalize_per_cell(adata)
    sc.preprocessing.pca(adata)
    sc.tools.tsne(adata=adata, n_pcs=5, perplexity=30)
    sc.pp.neighbors(adata, n_pcs=15, n_neighbors=10)
    sc.tl.diffmap(adata, n_comps=10)

    sce.tl.wishbone(
        adata=adata,
        start_cell='ACAAGAGACTTATC-1',
        components_list=[2, 3],
        num_waypoints=150,
    )
    assert all(
        [k in adata.obs for k in ['trajectory_wishbone', 'branch_wishbone']]
    ), "Run Wishbone Error!"

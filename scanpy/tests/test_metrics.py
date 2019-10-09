from operator import eq
import numpy as np
import scanpy as sc
from scipy import sparse


def test_gearys_c():
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.layers["raw"] = pbmc.raw.X.copy()
    g = pbmc.uns["neighbors"]["connectivities"]

    assert eq(
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"]),
        sc.metrics.gearys_c(pbmc, vals=pbmc.obs["percent_mito"])
    )

    assert eq(  # Test that series and vectors return same value
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"]),
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"].values),
    )

    assert np.array_equal(
        sc.metrics.gearys_c(pbmc, obsm="X_pca"),
        sc.metrics.gearys_c(g, pbmc.obsm["X_pca"].T)
    )

    # Test case with perfectly seperated groups
    connected = np.zeros(100)
    connected[np.random.choice(100, size=30, replace=False)] = 1
    graph = np.zeros((100, 100))
    graph[np.ix_(connected.astype(bool), connected.astype(bool))] = 1
    graph[np.ix_(~connected.astype(bool), ~connected.astype(bool))] = 1
    graph = sparse.csr_matrix(graph)

    assert sc.metrics.gearys_c(graph, connected) == 0.
    adata = sc.AnnData(
        sparse.csr_matrix((100, 100)), obsp={"connectivities": graph}
    )
    assert sc.metrics.gearys_c(adata, vals=connected) == 0.


def test_confusion_matrix():
    mtx = sc.metrics.confusion_matrix(["a", "b"], ["c", "d"], normalize=False)
    assert mtx.loc["a", "c"] == 1
    assert mtx.loc["a", "d"] == 0
    assert mtx.loc["b", "d"] == 1
    assert mtx.loc["b", "c"] == 0

    mtx = sc.metrics.confusion_matrix(["a", "b"], ["c", "d"], normalize=True)
    assert mtx.loc["a", "c"] == 1.
    assert mtx.loc["a", "d"] == 0
    assert mtx.loc["b", "d"] == 1.
    assert mtx.loc["b", "c"] == 0

    mtx = sc.metrics.confusion_matrix(["a", "a", "b", "b"], ["c", "d", "c", "d"], normalize=True)
    assert np.all(mtx == .5)
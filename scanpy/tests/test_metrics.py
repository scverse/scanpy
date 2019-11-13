from operator import eq
from string import ascii_letters

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def test_gearys_c():
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.layers["raw"] = pbmc.raw.X.copy()
    g = pbmc.uns["neighbors"]["connectivities"]

    assert eq(
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"]),
        sc.metrics.gearys_c(pbmc, vals=pbmc.obs["percent_mito"]),
    )

    assert eq(  # Test that series and vectors return same value
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"]),
        sc.metrics.gearys_c(g, pbmc.obs["percent_mito"].values),
    )

    assert np.array_equal(
        sc.metrics.gearys_c(pbmc, obsm="X_pca"),
        sc.metrics.gearys_c(g, pbmc.obsm["X_pca"].T),
    )

    all_genes = sc.metrics.gearys_c(pbmc, layer="raw")
    first_gene = sc.metrics.gearys_c(
        pbmc,
        vals=pbmc.obs_vector(pbmc.var_names[0], layer="raw")
    )

    assert np.allclose(all_genes[0], first_gene)

    assert np.allclose(
        sc.metrics.gearys_c(pbmc, layer="raw"),
        sc.metrics.gearys_c(pbmc, vals=pbmc.layers["raw"].T.toarray())
    )

    # Test case with perfectly seperated groups
    connected = np.zeros(100)
    connected[np.random.choice(100, size=30, replace=False)] = 1
    graph = np.zeros((100, 100))
    graph[np.ix_(connected.astype(bool), connected.astype(bool))] = 1
    graph[np.ix_(~connected.astype(bool), ~connected.astype(bool))] = 1
    graph = sparse.csr_matrix(graph)

    assert sc.metrics.gearys_c(graph, connected) == 0.0
    assert sc.metrics.gearys_c(graph, connected) \
        == sc.metrics.gearys_c(graph, sparse.csr_matrix(connected))
    adata = sc.AnnData(
        sparse.csr_matrix((100, 100)), obsp={"connectivities": graph}
    )
    assert sc.metrics.gearys_c(adata, vals=connected) == 0.0


def test_confusion_matrix():
    mtx = sc.metrics.confusion_matrix(["a", "b"], ["c", "d"], normalize=False)
    assert mtx.loc["a", "c"] == 1
    assert mtx.loc["a", "d"] == 0
    assert mtx.loc["b", "d"] == 1
    assert mtx.loc["b", "c"] == 0

    mtx = sc.metrics.confusion_matrix(["a", "b"], ["c", "d"], normalize=True)
    assert mtx.loc["a", "c"] == 1.0
    assert mtx.loc["a", "d"] == 0.0
    assert mtx.loc["b", "d"] == 1.0
    assert mtx.loc["b", "c"] == 0.0

    mtx = sc.metrics.confusion_matrix(
        ["a", "a", "b", "b"], ["c", "d", "c", "d"], normalize=True
    )
    assert np.all(mtx == 0.5)


def test_confusion_matrix_randomized():
    chars = np.array(list(ascii_letters))
    pos = np.random.choice(len(chars), size=np.random.randint(50, 150))
    a = chars[pos]
    b = np.random.permutation(chars)[pos]
    df = pd.DataFrame({"a": a, "b": b})

    pd.testing.assert_frame_equal(
        sc.metrics.confusion_matrix("a", "b", df),
        sc.metrics.confusion_matrix(df["a"], df["b"]),
    )
    pd.testing.assert_frame_equal(
        sc.metrics.confusion_matrix(df["a"].values, df["b"].values),
        sc.metrics.confusion_matrix(a, b),
    )

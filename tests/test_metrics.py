from __future__ import annotations

import warnings
from functools import partial
from string import ascii_letters
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import threadpoolctl
from scipy import sparse
from scipy.sparse import csr_matrix  # noqa: TID251

import scanpy as sc
from scanpy.metrics import modularity, modularity_adata
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES

if TYPE_CHECKING:
    from collections.abc import Generator


@pytest.fixture(scope="session", params=[sc.metrics.gearys_c, sc.metrics.morans_i])
def metric(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture(params=["single-threaded", "multi-threaded"])
def _threading(request: pytest.FixtureRequest) -> Generator[None, None, None]:
    if request.param == "single-threaded":
        with threadpoolctl.threadpool_limits(limits=1):
            yield
    elif request.param == "multi-threaded":
        yield


@pytest.mark.usefixtures("_threading")
def test_consistency(metric) -> None:
    pbmc = pbmc68k_reduced()
    pbmc.layers["raw"] = pbmc.raw.X.copy()
    g = pbmc.obsp["connectivities"]
    equality_check = partial(np.testing.assert_allclose, atol=1e-11)

    # This can fail
    equality_check(
        metric(g, pbmc.obs["percent_mito"]),
        metric(g, pbmc.obs["percent_mito"]),
    )
    equality_check(
        metric(g, pbmc.obs["percent_mito"]),
        metric(pbmc, vals=pbmc.obs["percent_mito"]),
    )

    equality_check(  # Test that series and vectors return same value
        metric(g, pbmc.obs["percent_mito"]),
        metric(g, pbmc.obs["percent_mito"].values),
    )

    equality_check(
        metric(pbmc, obsm="X_pca"),
        metric(g, pbmc.obsm["X_pca"].T),
    )

    all_genes = metric(pbmc, layer="raw")
    first_gene = metric(pbmc, vals=pbmc.obs_vector(pbmc.var_names[0], layer="raw"))

    np.testing.assert_allclose(all_genes[0], first_gene, rtol=1e-9)

    # Test that results are similar for sparse and dense reps of same data
    equality_check(
        metric(pbmc, layer="raw"),
        metric(pbmc, vals=pbmc.layers["raw"].T.toarray()),
    )


@pytest.mark.parametrize(
    ("metric", "size", "expected"),
    [
        pytest.param(sc.metrics.gearys_c, 30, 0.0, id="gearys_c"),
        pytest.param(sc.metrics.morans_i, 50, 1.0, id="morans_i"),
    ],
)
def test_correctness(metric, size, expected):
    # Test case with perfectly seperated groups
    connected = np.zeros(100)
    connected[np.random.choice(100, size=size, replace=False)] = 1
    graph = np.zeros((100, 100))
    graph[np.ix_(connected.astype(bool), connected.astype(bool))] = 1
    graph[np.ix_(~connected.astype(bool), ~connected.astype(bool))] = 1
    graph = sparse.csr_matrix(graph)  # noqa: TID251

    np.testing.assert_equal(metric(graph, connected), expected)
    np.testing.assert_equal(
        metric(graph, connected),
        metric(graph, sparse.csr_matrix(connected)),  # noqa: TID251
    )
    # Checking that obsp works
    adata = sc.AnnData(sparse.csr_matrix((100, 100)), obsp={"connectivities": graph})  # noqa: TID251
    np.testing.assert_equal(metric(adata, vals=connected), expected)


@pytest.mark.usefixtures("_threading")
@pytest.mark.parametrize(
    "array_type", [*ARRAY_TYPES, pytest.param(sparse.coo_matrix, id="scipy_coo")]
)
def test_graph_metrics_w_constant_values(
    request: pytest.FixtureRequest, metric, array_type
):
    if "dask" in array_type.__name__:
        reason = "DaskArray not yet supported"
        request.applymarker(pytest.mark.xfail(reason=reason))

    # https://github.com/scverse/scanpy/issues/1806
    pbmc = pbmc68k_reduced()
    XT = pbmc.raw.X.T.copy()
    g = pbmc.obsp["connectivities"].copy()
    equality_check = partial(np.testing.assert_allclose, atol=1e-11)

    const_inds = np.random.choice(XT.shape[0], 10, replace=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
        XT_zero_vals = XT.copy()
        XT_zero_vals[const_inds, :] = 0
        XT_zero_vals = array_type(XT_zero_vals)
        XT_const_vals = XT.copy()
        XT_const_vals[const_inds, :] = 42
        XT_const_vals = array_type(XT_const_vals)

    results_full = metric(g, array_type(XT))
    # TODO: Check for warnings
    with pytest.warns(
        UserWarning, match=r"10 variables were constant, will return nan for these"
    ):
        results_const_zeros = metric(g, XT_zero_vals)
    with pytest.warns(
        UserWarning, match=r"10 variables were constant, will return nan for these"
    ):
        results_const_vals = metric(g, XT_const_vals)

    assert not np.isnan(results_full).any()
    equality_check(results_const_zeros, results_const_vals)
    np.testing.assert_array_equal(np.nan, results_const_zeros[const_inds])
    np.testing.assert_array_equal(np.nan, results_const_vals[const_inds])

    non_const_mask = ~np.isin(np.arange(XT.shape[0]), const_inds)
    equality_check(results_full[non_const_mask], results_const_zeros[non_const_mask])


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


def test_confusion_matrix_api():
    data = pd.DataFrame(
        {"a": np.random.randint(5, size=100), "b": np.random.randint(5, size=100)}
    )
    expected = sc.metrics.confusion_matrix(data["a"], data["b"])

    pd.testing.assert_frame_equal(expected, sc.metrics.confusion_matrix("a", "b", data))

    pd.testing.assert_frame_equal(
        expected, sc.metrics.confusion_matrix("a", data["b"], data)
    )

    pd.testing.assert_frame_equal(
        expected, sc.metrics.confusion_matrix(data["a"], "b", data)
    )


# Test 1: Sample graph with clear community structure (dense & sparse, directed & undirected)
@pytest.mark.parametrize("is_directed", [False, True], ids=["undirected", "directed"])
@pytest.mark.parametrize("use_sparse", [False, True], ids=["sparse", "dense"])
@needs.igraph
def test_modularity_sample_structure(use_sparse, is_directed):
    # 4 node adjacency matrix with two separate 2-node communities
    mat = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
        ]
    )
    labels = ["A", "A", "B", "B"]
    adj = csr_matrix(mat) if use_sparse else mat
    score = modularity(adj, labels, is_directed=is_directed)

    # Modularity should be between 0 and 1
    assert 0 <= score <= 1


# Test 2: Edge case when all nodes belong to the same community/cluster
@needs.igraph
def test_modularity_single_community():
    # fully connected graph sample
    adj = np.ones((4, 4)) - np.eye(4)
    labels = ["A", "A", "A", "A"]
    score = modularity(adj, labels, is_directed=False)

    # modularity should be 0
    assert score == pytest.approx(0.0, rel=1e-6)


# Test 3: Invalad input, labels length does not match adjacency matrix size
@needs.igraph
def test_modularity_invalid_labels():
    from igraph._igraph import InternalError

    adj = np.eye(4)
    labels = ["A", "A", "B"]
    with pytest.raises(
        InternalError,
        match="Membership vector size differs",
    ):
        modularity(adj, labels, is_directed=False)


# Test 4: Pass both Louvain and Leiden clustering algorithms
@pytest.mark.parametrize("cluster_method", ["louvain", "leiden"])
@needs.igraph
@needs.louvain
@needs.leidenalg
def test_modularity_adata_multiple_clusterings(cluster_method):
    # Loading PBMC Data and compute PCA and neighbors graph
    adata = sc.datasets.pbmc3k()
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    # Compute modularity using both Louvain and Leiden clustering
    if cluster_method == "louvain":
        sc.tl.louvain(adata)
        score_louvain = modularity_adata(
            adata, labels="louvain", obsp="connectivities", is_directed=False
        )
        # Score should be between 0 and 1
        assert 0 <= score_louvain <= 1
    if cluster_method == "leiden":
        sc.tl.leiden(adata)
        score_leiden = modularity_adata(
            adata, labels="leiden", obsp="connectivities", is_directed=False
        )
        # Score should be between 0 and 1
        assert 0 <= score_leiden <= 1


# Test 5: Modularity should be the same no matter the order of the labels
@needs.igraph
def test_modularity_order():
    adj = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
        ]
    )
    labels1 = ["A", "A", "B", "B"]
    labels2 = ["B", "B", "A", "A"]
    score_1 = modularity(adj, labels1, is_directed=False)
    score_2 = modularity(adj, labels2, is_directed=False)
    assert score_1 == score_2


# Test 6: Modularity on disconnected graph lke edge-case behavior in some algorithms
@needs.igraph
def test_modularity_disconnected_graph():
    adj = np.zeros((4, 4))
    labels = ["A", "B", "C", "D"]
    score = modularity(adj, labels, is_directed=False)

    # Modularity should be undefined for disconnected graphs
    assert np.isnan(score)

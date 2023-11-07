import warnings
from functools import partial
from operator import eq
from string import ascii_letters

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import pytest

from scanpy._compat import DaskArray
from scanpy.testing._helpers.data import pbmc68k_reduced


mark_flaky = pytest.mark.xfail(
    strict=False,
    reason="This used to work reliably, but doesn’t anymore",
)


@pytest.fixture(scope="session", params=[sc.metrics.gearys_c, sc.metrics.morans_i])
def metric(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        pytest.param(eq, marks=[mark_flaky]),
        pytest.param(partial(np.testing.assert_allclose, rtol=1e-15), id="allclose"),
    ],
)
def assert_equal(request: pytest.FixtureRequest):
    return request.param


def test_consistency(metric, assert_equal):
    pbmc = pbmc68k_reduced()
    pbmc.layers["raw"] = pbmc.raw.X.copy()
    g = pbmc.obsp["connectivities"]

    assert_equal(
        metric(g, pbmc.obs["percent_mito"]),
        metric(pbmc, vals=pbmc.obs["percent_mito"]),
    )

    assert_equal(  # Test that series and vectors return same value
        metric(g, pbmc.obs["percent_mito"]),
        metric(g, pbmc.obs["percent_mito"].values),
    )

    np.testing.assert_almost_equal(
        metric(pbmc, obsm="X_pca"),
        metric(g, pbmc.obsm["X_pca"].T),
    )

    all_genes = metric(pbmc, layer="raw")
    first_gene = metric(pbmc, vals=pbmc.obs_vector(pbmc.var_names[0], layer="raw"))

    np.testing.assert_allclose(all_genes[0], first_gene)

    # Test that results are similar for sparse and dense reps of same data
    np.testing.assert_allclose(
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
def test_correctness(metric, size, expected, assert_equal):
    # Test case with perfectly seperated groups
    connected = np.zeros(100)
    connected[np.random.choice(100, size=size, replace=False)] = 1
    graph = np.zeros((100, 100))
    graph[np.ix_(connected.astype(bool), connected.astype(bool))] = 1
    graph[np.ix_(~connected.astype(bool), ~connected.astype(bool))] = 1
    graph = sparse.csr_matrix(graph)

    assert metric(graph, connected) == expected
    assert_equal(
        metric(graph, connected),
        metric(graph, sparse.csr_matrix(connected)),
    )
    # Checking that obsp works
    adata = sc.AnnData(sparse.csr_matrix((100, 100)), obsp={"connectivities": graph})
    assert metric(adata, vals=connected) == expected


def test_graph_metrics_w_constant_values(metric, array_type, assert_equal):
    # https://github.com/scverse/scanpy/issues/1806
    pbmc = pbmc68k_reduced()
    XT = array_type(pbmc.raw.X.T.copy())
    g = pbmc.obsp["connectivities"].copy()

    if isinstance(XT, DaskArray):
        pytest.skip("DaskArray yet not supported")

    const_inds = np.random.choice(XT.shape[0], 10, replace=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
        XT_zero_vals = XT.copy()
        XT_zero_vals[const_inds, :] = 0
        XT_const_vals = XT.copy()
        XT_const_vals[const_inds, :] = 42

    results_full = metric(g, XT)
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
    assert_equal(results_const_zeros, results_const_vals)
    np.testing.assert_array_equal(np.nan, results_const_zeros[const_inds])
    np.testing.assert_array_equal(np.nan, results_const_vals[const_inds])

    non_const_mask = ~np.isin(np.arange(XT.shape[0]), const_inds)
    assert_equal(results_full[non_const_mask], results_const_zeros[non_const_mask])


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

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

import scanpy as sc
from scanpy import Neighbors
from scanpy._compat import CSBase
from testing.scanpy._helpers import anndata_v0_8_constructor_compat

if TYPE_CHECKING:
    from typing import Literal

    from pytest_mock import MockerFixture

# the input data
X = [[1, 0], [3, 0], [5, 6], [0, 4]]
n_neighbors = 3  # includes data points themselves

# distances
distances_euclidean = [
    [0.0, 2.0, 0.0, 4.123105525970459],
    [2.0, 0.0, 0.0, 5.0],
    [0.0, 6.324555397033691, 0.0, 5.385164737701416],
    [4.123105525970459, 5.0, 0.0, 0.0],
]

distances_euclidean_all = [
    [0.0, 2.0, 7.211102485656738, 4.123105525970459],
    [2.0, 0.0, 6.324555397033691, 5.0],
    [7.211102485656738, 6.324555397033691, 0.0, 5.385164737701416],
    [4.123105525970459, 5.0, 5.385164737701416, 0.0],
]


# umap "kernel" – only knn results
connectivities_umap = [
    [0.0, 1.0, 0.0, 1.0],
    [1.0, 0.0, 0.5849691143165735, 0.8277419907567016],
    [0.0, 0.5849691143165735, 0.0, 1.0],
    [1.0, 0.8277419907567016, 1.0, 0.0],
]

transitions_sym_umap = [
    [0.0, 0.4818987107873648, 0.0, 0.3951883393150153],
    [0.48189871078736474, 0.0, 0.3594582764005241, 0.24216345431293487],
    [0.0, 0.3594582764005241, 0.0, 0.5039226836320637],
    [0.39518833931501524, 0.24216345431293487, 0.5039226836320637, 0.0],
]

transitions_umap = [
    [0.0, 0.5395987596963403, 0.0, 0.4604012403036599],
    [0.430368608684738, 0.0, 0.3176747629691457, 0.2519566283461165],
    [0.0, 0.40673754271561435, 0.0, 0.5932624572843856],
    [0.33921243006981133, 0.23275092618009624, 0.42803664375009237, 0.0],
]


# gauss kernel [diffmap, dpt] – knn and dense results
connectivities_gauss_knn = [
    [0.0, 0.8466368913650513, 0.0, 0.5660185813903809],
    [0.8466368913650513, 0.0, 0.4223647117614746, 0.4902938902378082],
    [0.0, 0.4223647117614746, 0.0, 0.5840492248535156],
    [0.5660185813903809, 0.4902938902378082, 0.5840492248535156, 0.0],
]

connectivities_gauss_noknn = [
    [1.0, 0.676927387714386, 0.024883469566702843, 0.1962655782699585],
    [0.676927387714386, 1.0, 0.08414449542760849, 0.1353352814912796],
    [0.024883469566702843, 0.08414449542760849, 1.0, 0.16558068990707397],
    [0.1962655782699585, 0.1353352814912796, 0.16558068990707397, 1.0],
]

transitions_sym_gauss_knn = [
    [0.0, 0.5146393179893494, 0.0, 0.36445462703704834],
    [0.5146393179893494, 0.0, 0.3581143319606781, 0.2239987552165985],
    [0.0, 0.3581143319606781, 0.0, 0.5245543718338013],
    [0.36445462703704834, 0.2239987552165985, 0.5245543718338013, 0.0],
]

transitions_sym_gauss_noknn = [
    [
        0.5093212127685547,
        0.34393802285194397,
        0.016115963459014893,
        0.11607448011636734,
    ],
    [0.34393805265426636, 0.506855845451355, 0.054364752024412155, 0.07984541356563568],
    [
        0.016115965321660042,
        0.054364752024412155,
        0.8235670328140259,
        0.12452481687068939,
    ],
    [0.11607448011636734, 0.07984541356563568, 0.1245248094201088, 0.6867417693138123],
]

transitions_gauss_knn = [
    [0.0, 0.5824036598205566, 0.0, 0.4175964295864105],
    [0.4547595679759979, 0.0, 0.3184431493282318, 0.22679725289344788],
    [0.0, 0.4027276933193207, 0.0, 0.5972723364830017],
    [0.3180755078792572, 0.22123482823371887, 0.46068981289863586, 0.0],
]

transitions_gauss_noknn = [
    [0.5093212127685547, 0.3450769782066345, 0.01887294091284275, 0.12672874331474304],
    [0.34280285239219666, 0.506855845451355, 0.06345486640930176, 0.08688655495643616],
    [0.01376173086464405, 0.04657683148980141, 0.8235670328140259, 0.11609435081481934],
    [0.10631592571735382, 0.07337487488985062, 0.13356748223304749, 0.6867417693138123],
]


def get_neighbors() -> Neighbors:
    return Neighbors(anndata_v0_8_constructor_compat(np.array(X)))


@pytest.fixture
def neigh() -> Neighbors:
    return get_neighbors()


@pytest.mark.parametrize("method", ["umap", "gauss"])
def test_distances_euclidean(
    mocker: MockerFixture, neigh: Neighbors, method: Literal["umap", "gauss"]
):
    """Umap and gauss behave the same for distances.

    They call pynndescent for large data.
    """
    from pynndescent import NNDescent

    # When trying to compress a too-small index, pynndescent complains
    mocker.patch.object(NNDescent, "compress_index", return_val=None)

    neigh.compute_neighbors(n_neighbors, method=method)
    np.testing.assert_allclose(neigh.distances.toarray(), distances_euclidean)


@pytest.mark.parametrize(
    ("transformer", "knn"),
    [
        # knn=False trivially returns all distances
        pytest.param(None, False, id="knn=False"),
        # pynndescent returns all distances when data is so small
        pytest.param("pynndescent", True, id="pynndescent"),
        # Explicit brute force also returns all distances
        pytest.param(
            KNeighborsTransformer(n_neighbors=n_neighbors, algorithm="brute"),
            True,
            id="sklearn",
        ),
    ],
)
def test_distances_all(neigh: Neighbors, transformer, knn):
    neigh.compute_neighbors(
        n_neighbors, transformer=transformer, method="gauss", knn=knn
    )
    dists = (
        neigh.distances.toarray()
        if isinstance(neigh.distances, CSBase)
        else neigh.distances
    )
    np.testing.assert_allclose(dists, distances_euclidean_all)


@pytest.mark.parametrize(
    ("method", "conn", "trans", "trans_sym"),
    [
        pytest.param(
            "umap",
            connectivities_umap,
            transitions_umap,
            transitions_sym_umap,
            id="umap",
        ),
        pytest.param(
            "gauss",
            connectivities_gauss_knn,
            transitions_gauss_knn,
            transitions_sym_gauss_knn,
            id="gauss",
        ),
    ],
)
def test_connectivities_euclidean(neigh: Neighbors, method, conn, trans, trans_sym):
    neigh.compute_neighbors(n_neighbors, method=method)
    np.testing.assert_allclose(neigh.connectivities.toarray(), conn)
    neigh.compute_transitions()
    np.testing.assert_allclose(neigh.transitions_sym.toarray(), trans_sym, rtol=1e-5)
    np.testing.assert_allclose(neigh.transitions.toarray(), trans, rtol=1e-5)


def test_gauss_noknn_connectivities_euclidean(neigh):
    neigh.compute_neighbors(n_neighbors, method="gauss", knn=False)
    np.testing.assert_allclose(neigh.connectivities, connectivities_gauss_noknn)
    neigh.compute_transitions()
    np.testing.assert_allclose(
        neigh.transitions_sym, transitions_sym_gauss_noknn, rtol=1e-5
    )
    np.testing.assert_allclose(neigh.transitions, transitions_gauss_noknn, rtol=1e-5)


def test_metrics_argument():
    no_knn_euclidean = get_neighbors()
    no_knn_euclidean.compute_neighbors(
        n_neighbors, method="gauss", knn=False, metric="euclidean"
    )
    no_knn_manhattan = get_neighbors()
    no_knn_manhattan.compute_neighbors(
        n_neighbors, method="gauss", knn=False, metric="manhattan"
    )
    assert not np.allclose(no_knn_euclidean.distances, no_knn_manhattan.distances)


def test_use_rep_argument():
    adata = AnnData(np.random.randn(30, 300))
    sc.pp.pca(adata)
    neigh_pca = Neighbors(adata)
    neigh_pca.compute_neighbors(n_pcs=5, use_rep="X_pca")
    neigh_none = Neighbors(adata)
    neigh_none.compute_neighbors(n_pcs=5, use_rep=None)
    np.testing.assert_allclose(
        neigh_pca.distances.toarray(), neigh_none.distances.toarray()
    )


@pytest.mark.parametrize("conv", [sparse.csr_matrix.toarray, sparse.csr_matrix])  # noqa: TID251
def test_restore_n_neighbors(neigh, conv):
    neigh.compute_neighbors(n_neighbors, method="gauss")

    ad = AnnData(np.array(X))
    # Allow deprecated usage for now
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")
        ad.uns["neighbors"] = dict(connectivities=conv(neigh.connectivities))
    neigh_restored = Neighbors(ad)
    assert neigh_restored.n_neighbors == 1

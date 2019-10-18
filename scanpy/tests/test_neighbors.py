import numpy as np
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from scanpy import Neighbors

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
    [0.5660185813903809, 0.4902938902378082, 0.5840492248535156, 0.0]]

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
    [0.5093212127685547, 0.34393802285194397, 0.016115963459014893, 0.11607448011636734],
    [0.34393805265426636, 0.506855845451355, 0.054364752024412155, 0.07984541356563568],
    [0.016115965321660042, 0.054364752024412155, 0.8235670328140259, 0.12452481687068939],
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
    return Neighbors(AnnData(np.array(X)))


@pytest.fixture
def neigh() -> Neighbors:
    return get_neighbors()


def test_umap_connectivities_euclidean(neigh):
    neigh.compute_neighbors(method='umap', n_neighbors=n_neighbors)
    assert np.allclose(
        neigh.distances.toarray(), distances_euclidean)
    assert np.allclose(
        neigh.connectivities.toarray(), connectivities_umap)
    neigh.compute_transitions()
    assert np.allclose(neigh.transitions_sym.toarray(), transitions_sym_umap)
    assert np.allclose(neigh.transitions.toarray(), transitions_umap)


def test_gauss_noknn_connectivities_euclidean(neigh):
    neigh.compute_neighbors(method='gauss', knn=False, n_neighbors=3)
    assert np.allclose(
        neigh.distances, distances_euclidean_all)
    assert np.allclose(
        neigh.connectivities, connectivities_gauss_noknn)
    neigh.compute_transitions()
    assert np.allclose(neigh.transitions_sym, transitions_sym_gauss_noknn)
    assert np.allclose(neigh.transitions, transitions_gauss_noknn)


def test_gauss_connectivities_euclidean(neigh):
    neigh.compute_neighbors(method='gauss', n_neighbors=n_neighbors)
    assert np.allclose(
        neigh.distances.toarray(), distances_euclidean)
    assert np.allclose(
        neigh.connectivities.toarray(), connectivities_gauss_knn)
    neigh.compute_transitions()
    assert np.allclose(neigh.transitions_sym.toarray(), transitions_sym_gauss_knn)
    assert np.allclose(neigh.transitions.toarray(), transitions_gauss_knn)


def test_metrics_argument():
    no_knn_euclidean = get_neighbors()
    no_knn_euclidean.compute_neighbors(method="gauss", knn=False,
        n_neighbors=n_neighbors, metric="euclidean")
    no_knn_manhattan = get_neighbors()
    no_knn_manhattan.compute_neighbors(method="gauss", knn=False,
        n_neighbors=n_neighbors, metric="manhattan")
    assert not np.allclose(no_knn_euclidean.distances, no_knn_manhattan.distances)


@pytest.mark.parametrize('conv', [csr_matrix.toarray, csr_matrix])
def test_restore_n_neighbors(neigh, conv):
    neigh.compute_neighbors(method='gauss', n_neighbors=n_neighbors)

    ad = AnnData(np.array(X))
    ad.uns['neighbors'] = dict(connectivities=conv(neigh.connectivities))
    neigh_restored = Neighbors(ad)
    assert neigh_restored.n_neighbors == 1

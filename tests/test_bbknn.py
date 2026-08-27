from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from anndata.acc import A
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

import scanpy as sc
from scanpy.get.get import _rep_from_json
from scanpy.neighbors._bbknn import (
    _compute_batch_balanced_knn,
    _handle_transformer,
    _trim,
)
from testing.scanpy._pytest.params import ARRAY_TYPES_MEM

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from scanpy._compat import CSRBase
    from scanpy.neighbors._types import _Metric


N_PER_BATCH = [60, 40, 30]
BATCHES = np.repeat(["a", "b", "c"], N_PER_BATCH)
N_OBS = len(BATCHES)


@pytest.fixture(
    scope="module", params=[10, 2 * sc.settings.N_PCS], ids=["narrow", "wide"]
)
def rep(request: pytest.FixtureRequest) -> NDArray[np.float32]:
    """Create a representation where each batch is shifted by a constant.

    Either narrow enough for `pp.bbknn` to use `.X`,
    or wide enough that it uses the PCA.
    """
    rng = np.random.default_rng(0)
    x = rng.normal(size=(N_OBS, request.param)).astype(np.float32)
    for i, batch in enumerate(np.unique(BATCHES)):
        x[batch == BATCHES] += 5 * i
    return x


@pytest.fixture
def adata(rep: NDArray[np.float32]) -> AnnData:
    adata = AnnData(rep.copy(), obs=dict(batch=BATCHES.copy()))
    # fewer PCs than `.X` has columns, so the two differ
    sc.pp.pca(adata, n_comps=5, key_added="pca")
    return adata


def n_per_row(a: CSRBase) -> NDArray[np.int64]:
    return np.diff(a.indptr)


def test_bbknn(adata: AnnData) -> None:
    assert sc.pp.bbknn(adata, 3, batches="obs.batch") is None

    params = adata.uns["neighbors"]["params"]
    assert params == dict(
        n_neighbors=9,  # 3 neighbors × 3 batches
        method="umap",
        metric="euclidean",
        batches=["obs", "batch"],
        neighbors_within_batch=3,
        trim=90,
    )
    dists, conns = adata.obsp["distances"], adata.obsp["connectivities"]
    assert dists.shape == conns.shape == (N_OBS, N_OBS)
    # the cell itself is not part of its own neighborhood
    assert (n_per_row(dists) == 8).all()
    assert dists.diagonal().sum() == 0
    assert (conns != conns.T).nnz == 0


def test_bbknn_representation(adata: AnnData) -> None:
    dists = sc.pp.bbknn(adata, 3, batches="obs.batch", copy=True).obsp["distances"]
    # like `pp.neighbors`, we use the PCA – except for data narrower than `N_PCS`
    used, unused = ("X", "pca") if adata.n_vars <= sc.settings.N_PCS else ("pca", "X")

    reps = dict(X=A.X, pca=A.obsm["pca"])
    for name in (used, unused):
        sc.pp.bbknn(adata, 3, batches="obs.batch", use_rep=reps[name], key_added=name)

    # accessors round trip through `.uns` for readers like `tl.umap`
    assert _rep_from_json(adata.uns["pca"]["params"]["use_rep"]) == A.obsm["pca"]

    np.testing.assert_allclose(
        dists.toarray(), adata.obsp[f"{used}_distances"].toarray()
    )
    assert (dists != adata.obsp[f"{unused}_distances"]).nnz > 0


def test_bbknn_is_batch_balanced(adata: AnnData) -> None:
    """Each cell has `neighbors_within_batch` neighbors in each batch, including itself."""
    sc.pp.bbknn(adata, 3, batches="obs.batch")

    dists = adata.obsp["distances"].tolil()
    for i, row in enumerate(dists.rows):
        neighbors = np.asarray([*row, i])  # add back the cell itself
        _, counts = np.unique(BATCHES[neighbors], return_counts=True)
        assert (counts == 3).all()


def test_bbknn_connects_batches(adata: AnnData) -> None:
    """Unlike `pp.neighbors`, `pp.bbknn` connects the (strongly separated) batches."""
    sc.pp.neighbors(adata, n_neighbors=9, key_added="knn")
    sc.pp.bbknn(adata, 3, batches="obs.batch", key_added="bbknn")

    def n_cross_batch(key: str) -> int:
        i, j = adata.obsp[f"{key}_connectivities"].nonzero()
        return int((BATCHES[i] != BATCHES[j]).sum())

    assert n_cross_batch("knn") == 0
    assert n_cross_batch("bbknn") > 0


@pytest.mark.parametrize(
    "transformer",
    [
        pytest.param(None, id="none"),
        pytest.param("sklearn", id="sklearn"),
        pytest.param("pynndescent", id="pynndescent"),
        pytest.param(
            KNeighborsTransformer(n_neighbors=3, algorithm="kd_tree"), id="instance"
        ),
    ],
)
def test_bbknn_transformer(adata: AnnData, transformer) -> None:
    sc.pp.bbknn(adata, 3, batches="obs.batch", transformer=transformer)
    assert (n_per_row(adata.obsp["distances"]) == 8).all()


@pytest.mark.parametrize(
    ("n_obs", "max_batch_size", "metric", "brute"),
    [
        # brute force costs `n_obs` × index size, so both matter
        pytest.param(20_000, 2_000, "euclidean", True, id="many_small_batches"),
        pytest.param(100_000, 50_000, "euclidean", False, id="few_big_batches"),
        pytest.param(300_000, 30_000, "euclidean", False, id="big_data"),
        pytest.param(1_000, 500, "cosine", True, id="small_batch_other_metric"),
        pytest.param(300_000, 30_000, "cosine", False, id="big_data_other_metric"),
    ],
)
def test_bbknn_transformer_choice(
    *, n_obs: int, max_batch_size: int, metric: _Metric, brute: bool
) -> None:
    """`transformer=None` picks a backend based on how big the per-batch indices are."""
    from sklearn.neighbors import KNeighborsTransformer

    transformer, shortcut = _handle_transformer(
        None,
        n_obs=n_obs,
        max_batch_size=max_batch_size,
        n_neighbors=3,
        metric=metric,
        metric_kwds={},
        rng=np.random.default_rng(0),
    )
    assert shortcut is brute
    assert isinstance(transformer, KNeighborsTransformer) is brute


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_bbknn_array_types(rep: NDArray[np.float32], array_type: Callable) -> None:
    adata = AnnData(array_type(np.abs(rep)), obs=dict(batch=BATCHES.copy()))
    sc.pp.bbknn(adata, 3, 0, batches="obs.batch")
    assert (n_per_row(adata.obsp["distances"]) == 8).all()


@pytest.mark.parametrize("trim", [None, 0, 5, 12])
def test_bbknn_trim(adata: AnnData, trim: int | None) -> None:
    sc.pp.bbknn(adata, 3, batches="obs.batch", trim=trim)
    conns = adata.obsp["connectivities"]

    assert adata.uns["neighbors"]["params"]["trim"] == (90 if trim is None else trim)
    if trim:
        # ties are kept, so cells can end up with slightly more than `trim` neighbors
        assert n_per_row(conns).max() >= trim
        assert (conns != conns.T).nnz == 0
    # trimming only ever removes edges
    sc.pp.bbknn(adata, 3, batches="obs.batch", trim=0, key_added="untrimmed")
    untrimmed = adata.obsp["untrimmed_connectivities"]
    assert conns.nnz <= untrimmed.nnz
    assert (conns != conns.multiply(untrimmed != 0)).nnz == 0


def test_trim() -> None:
    """`_trim` cuts each row at its `trim`-th largest value, but keeps the graph symmetric."""
    dense = [
        [0.0, 0.9, 0.8, 0.7],
        [0.9, 0.0, 0.1, 0.0],
        [0.8, 0.1, 0.0, 0.0],
        [0.7, 0.0, 0.0, 0.0],
    ]
    conns = sparse.csr_matrix(dense)  # noqa: TID251
    trimmed = _trim(conns.copy(), 2).toarray()
    # row 0 keeps its top 2 (0.9, 0.8); 0.7 is dropped in both directions
    np.testing.assert_allclose(
        trimmed,
        [
            [0.0, 0.9, 0.8, 0.0],
            [0.9, 0.0, 0.1, 0.0],
            [0.8, 0.1, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ],
    )
    # rows with at most `trim` entries are left alone
    np.testing.assert_allclose(_trim(conns.copy(), 4).toarray(), conns.toarray())


def test_bbknn_knn_is_normalized(rep: NDArray[np.float32]) -> None:
    """The merged neighbors are sorted, and each cell is its own first neighbor.

    kNN backends can report a tiny non-zero distance of a cell to itself,
    which `umap` would mistake for the radius of that cell’s local neighborhood.
    """
    knn_indices, knn_distances = _compute_batch_balanced_knn(
        rep,
        batches=BATCHES,
        unique_batches=np.unique(BATCHES),
        batch_sizes=np.asarray(N_PER_BATCH),
        neighbors_within_batch=3,
        transformer=None,
        metric="euclidean",
        metric_kwds={},
        rng=np.random.default_rng(0),
    )
    np.testing.assert_array_equal(knn_indices[:, 0], np.arange(N_OBS))
    np.testing.assert_array_equal(knn_distances[:, 0], 0.0)
    assert (np.diff(knn_distances, axis=1) >= 0).all()


def test_bbknn_duplicate_cells() -> None:
    """Duplicate cells are at distance 0, but a cell is still its own first neighbor."""
    rng = np.random.default_rng(0)
    x = rng.normal(size=(20, 5))
    x[10:] = x[:10]  # each cell has an exact duplicate
    batches = np.repeat(["a", "b"], 10)

    knn_indices, knn_distances = _compute_batch_balanced_knn(
        x,
        batches=batches,
        unique_batches=np.unique(batches),
        batch_sizes=np.asarray([10, 10]),
        neighbors_within_batch=2,
        transformer=None,
        metric="euclidean",
        metric_kwds={},
        rng=np.random.default_rng(0),
    )
    np.testing.assert_array_equal(knn_indices[:, 0], np.arange(20))
    np.testing.assert_array_equal(knn_distances[:, 0], 0.0)

    adata = AnnData(x, obs=dict(batch=batches))
    sc.pp.bbknn(adata, 2, 0, batches="obs.batch")
    assert adata.obsp["distances"].diagonal().sum() == 0


def test_bbknn_key_added(adata: AnnData) -> None:
    sc.pp.bbknn(adata, 3, batches="obs.batch")
    sc.pp.bbknn(adata, 3, batches="obs.batch", key_added="test")

    assert adata.uns["neighbors"]["params"] == adata.uns["test"]["params"]
    assert adata.uns["test"]["distances_key"] == "test_distances"
    assert adata.uns["test"]["connectivities_key"] == "test_connectivities"
    for key in ("distances", "connectivities"):
        np.testing.assert_allclose(
            adata.obsp[key].toarray(), adata.obsp[f"test_{key}"].toarray()
        )


def test_bbknn_copy(adata: AnnData) -> None:
    copied = sc.pp.bbknn(adata, 3, batches="obs.batch", copy=True)
    assert not adata.obsp
    assert "neighbors" not in adata.uns  # `.uns['pca']` is from the fixture
    assert set(copied.obsp) == {"distances", "connectivities"}


@pytest.mark.parametrize(
    ("kwargs", "error", "pattern"),
    [
        pytest.param(
            dict(batches="obs.nope"),
            KeyError,
            r"Batch key A.obs\['nope'\] not found",
            id="key",
        ),
        pytest.param(
            dict(neighbors_within_batch=0),
            ValueError,
            r"needs to be greater than 0",
            id="n_neighbors",
        ),
        pytest.param(
            dict(neighbors_within_batch=40),
            ValueError,
            r"Not all batches have at least .* \['c'\]",
            id="batch_too_small",
        ),
        pytest.param(
            dict(transformer="nope"),
            ValueError,
            r"Unknown transformer",
            id="transformer",
        ),
    ],
)
def test_bbknn_errors(
    adata: AnnData, kwargs: dict, error: type[Exception], pattern: str
) -> None:
    with pytest.raises(error, match=pattern):
        sc.pp.bbknn(adata, **{"batches": "obs.batch", **kwargs})

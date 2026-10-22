from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from sklearn.neighbors import KNeighborsTransformer, NearestNeighbors

import scanpy as sc
from scanpy.neighbors import _wnn
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray


def _ref_jaccard(nn: NDArray[np.integer], prune: float = 0.0) -> NDArray[np.float64]:
    """Seurat’s ComputeSNN."""
    a = np.zeros((len(nn), len(nn)))
    a[np.arange(len(nn))[:, None], nn] = 1
    shared = a @ a.T
    j = shared / (2 * nn.shape[1] - shared)
    return np.where(j < prune, 0, j)


def _ref_bandwidth(
    d: NDArray[np.float64], nn: NDArray[np.integer], nearest: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Seurat’s ComputeSNNwidth + SNN_SmallestNonzero_Dist."""
    sigma = np.empty(len(nn))
    for i, row in enumerate(_ref_jaccard(nn)):
        partners = np.flatnonzero(row)
        n_i = min(nn.shape[1], partners.size)
        thresh = np.sort(row[partners])[n_i - 1]
        dist = np.maximum(d[i, partners[row[partners] <= thresh]] - nearest[i], 0)
        sigma[i] = np.sort(dist)[::-1][:n_i].mean()
    return sigma


def _seurat_wnn_reference(
    reps: list[NDArray[np.float64]],
    *,
    k: int,
    knn_range: int,
    cross: float = 1e-4,
    smooth: bool = False,
):
    """Straight NumPy port of Seurat’s FindMultiModalNeighbors (exact kNN)."""
    e = [x / np.linalg.norm(x, axis=1, keepdims=True) for x in reps]
    n, n_mod = e[0].shape[0], len(e)
    rows = np.arange(n)
    d = [np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2) for x in e]
    order = [np.argsort(dm, axis=1, kind="stable") for dm in d]
    nn = [o[:, :k] for o in order]
    nearest = [d[m][rows, nn[m][:, 1]] for m in range(n_mod)]
    sigma = [_ref_bandwidth(d[m], nn[m], nearest[m]) for m in range(n_mod)]

    def impute(m, r):  # PredictAssay + impute_dist
        pred = e[m][nn[r][:, 1:]].mean(axis=1)
        return np.maximum(np.linalg.norm(e[m] - pred, axis=1) - nearest[m], 0)

    exp_scores = np.zeros((n_mod, n))
    for m in range(n_mod):
        within = np.exp(-impute(m, m) / sigma[m])
        for r in range(n_mod):
            if r == m:
                continue
            score = np.clip(within / (np.exp(-impute(m, r) / sigma[m]) + cross), 0, 200)
            if smooth:
                score = score[nn[m][:, 1:]].mean(axis=1)
            exp_scores[m] += np.exp(score)
    weights = exp_scores / exp_scores.sum(axis=0)

    idx = np.empty((n, k), dtype=int)
    sim = np.empty((n, k))
    for i in range(n):  # MultiModalNN
        cand = np.concatenate([o[i, 1:knn_range] for o in order])
        cand = np.array(list(dict.fromkeys(cand)))
        s = sum(
            weights[m, i]
            * np.exp(-np.maximum(d[m][i, cand] - nearest[m][i], 0) / sigma[m][i])
            for m in range(n_mod)
        )
        best = np.argsort(-s, kind="stable")[:k]
        idx[i], sim[i] = cand[best], s[best]
    snn = _ref_jaccard(idx, prune=1 / 15)
    return weights, idx, np.sqrt(np.clip((1 - sim) / 2, 0, 1)), snn


def _assert_matches(adata: AnnData, names: list[str], ref) -> None:
    weights, idx, dist, snn = ref
    w = adata.obs[[f"{name}_weight" for name in names]].to_numpy()
    np.testing.assert_allclose(w, np.transpose(weights), rtol=1e-6, atol=1e-8)
    d = adata.obsp["distances"]
    # rows keep the neighbors in decreasing similarity, as Seurat does
    np.testing.assert_array_equal(d.indices.reshape(idx.shape), idx)
    np.testing.assert_allclose(d.data.reshape(idx.shape), dist, atol=1e-6)
    np.testing.assert_allclose(adata.obsp["snn"].toarray(), snn, atol=1e-7)


def _multimodal_adata(n: int = 400) -> AnnData:
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 4, n)
    adata = AnnData(np.zeros((n, 1), dtype=np.float32))
    for key, d, n_groups in (("X_pca", 12, 4), ("X_apca", 8, 2), ("X_lsi", 6, 3)):
        centers = rng.normal(size=(4, d)) * 3
        x = centers[labels % n_groups] + rng.normal(size=(n, d))
        adata.obsm[key] = x.astype(np.float32)
    return adata


@pytest.mark.parametrize("smooth", [False, True])
@pytest.mark.parametrize("n_mod", [2, 3])
def test_wnn_matches_seurat_reference(*, smooth: bool, n_mod: int) -> None:
    adata = _multimodal_adata()
    keys = ["X_pca", "X_apca", "X_lsi"][:n_mod]
    # NumPy integers, e.g. from a parameter sweep, are accepted
    sc.pp.wnn(adata, keys, n_neighbors=np.int64(10), knn_range=40, smooth=smooth)
    reps = [adata.obsm[key].astype(np.float64) for key in keys]
    ref = _seurat_wnn_reference(reps, k=10, knn_range=40, smooth=smooth)
    _assert_matches(adata, keys, ref)


def test_wnn_matches_seurat_fixture(data_dir: Path) -> None:
    """Seurat 5.5.1 with exact RANN search, see `tests/_scripts/generate_wnn_seurat.py`."""
    ref = np.load(data_dir / "wnn_seurat_5.5.1.npz")
    adata = AnnData(np.zeros((200, 1), dtype=np.float32))
    adata.obsm["X_pca"], adata.obsm["X_apca"] = ref["pca"], ref["apca"]
    sc.pp.wnn(adata, {"RNA": "X_pca", "ADT": "X_apca"}, knn_range=40)
    ref = ref["weights"].T, ref["indices"], ref["distances"], ref["snn"]
    _assert_matches(adata, ["RNA", "ADT"], ref)


@pytest.mark.parametrize("k", [5, 30, 70])
def test_wnn_snn(k: int) -> None:
    """SNN bandwidths and graph match Seurat’s definitions."""
    e = _multimodal_adata(n=300).obsm["X_pca"].astype(np.float64)
    e /= np.linalg.norm(e, axis=1, keepdims=True)
    d = np.linalg.norm(e[:, None] - e[None], axis=2)
    nn = np.argsort(d, axis=1, kind="stable")[:, :k].astype(np.int32)
    nearest = d[np.arange(len(d)), nn[:, 1]]

    post_off, postings = _wnn._postings(nn, k)
    sigma = _wnn._snn_bandwidth(nn, k, post_off, postings, e, nearest)
    np.testing.assert_allclose(sigma, _ref_bandwidth(d, nn, nearest), rtol=1e-12)
    snn = _wnn._snn_graph(nn, post_off, postings, prune=1 / 15)
    assert snn.has_canonical_format
    np.testing.assert_allclose(snn.toarray(), _ref_jaccard(nn, 1 / 15), atol=1e-7)


@needs.igraph
def test_wnn_outputs_h5ad_and_downstream(tmp_path: Path) -> None:
    adata = _multimodal_adata()
    adata.X = None  # the layout must use an existing representation
    use_rep = {"RNA": "X_pca", "ADT": "X_apca"}
    out = sc.pp.wnn(
        adata,
        use_rep,
        n_pcs=[10, None],
        n_neighbors=15,
        knn_range=50,
        key_added="wnn",
        copy=True,
    )
    assert out is not None
    assert "wnn" not in adata.uns
    params = out.uns["wnn"]["params"]
    assert (params["use_rep"], params["n_pcs"]) == ("X_pca", 10)
    # scanpy’s convention counts the cell itself
    assert params["n_neighbors"] == 16
    assert params["random_state"] == 0
    assert np.all(np.diff(out.obsp["wnn_distances"].indptr) == 15)
    assert params["modalities"] == {
        "RNA": {"use_rep": "X_pca", "n_pcs": 10},
        "ADT": {"use_rep": "X_apca", "n_pcs": 8},
    }
    w = out.obs[["RNA_weight", "ADT_weight"]].to_numpy()
    np.testing.assert_allclose(w.sum(axis=1), 1)
    # the RNA-like modality separates more groups, so it should dominate
    assert w[:, 0].mean() > w[:, 1].mean()
    c = out.obsp["wnn_connectivities"]
    assert abs(c - c.T).max() < 1e-6
    assert c.diagonal().max() == 0
    out.write_h5ad(tmp_path / "wnn.h5ad")
    loaded = read_h5ad(tmp_path / "wnn.h5ad")
    assert loaded.uns["wnn"]["params"]["modalities"] == params["modalities"]
    for key in ("wnn_distances", "wnn_connectivities", "wnn_snn"):
        assert (loaded.obsp[key] != out.obsp[key]).nnz == 0
    sc.tl.leiden(out, neighbors_key="wnn", flavor="igraph")
    sc.tl.leiden(out, obsp="wnn_snn", key_added="leiden_snn", flavor="igraph")
    sc.tl.umap(loaded, neighbors_key="wnn", maxiter=10)
    assert np.isfinite(loaded.obsm["X_umap"]).all()
    # ingest would silently map through the first modality only
    with pytest.raises(NotImplementedError, match="multimodal"):
        sc.tl.ingest(out[:50].copy(), out, obs="leiden", neighbors_key="wnn")


@pytest.mark.parametrize("method", ["gauss", "jaccard"])
def test_wnn_connectivity_methods(method: str) -> None:
    adata = _multimodal_adata()
    keys = ["X_pca", "X_apca"]
    sc.pp.wnn(adata, keys, n_pcs=5, knn_range=50, method=method, weight_key=None)
    params = adata.uns["neighbors"]["params"]
    assert (params["method"], params["modalities"]["X_apca"]["n_pcs"]) == (method, 5)
    assert adata.obs.columns.empty
    c = adata.obsp["connectivities"]
    assert c.nnz > 0
    assert abs(c - c.T).max() < 1e-6


def test_wnn_rng() -> None:
    """The approximate search is reproducible, also with the legacy `random_state`."""
    adata = _multimodal_adata()
    kw = dict(use_rep=["X_pca", "X_apca"], knn_range=50, transformer="pynndescent")
    a, b = (sc.pp.wnn(adata, **kw, copy=True) for _ in range(2))
    c = sc.pp.wnn(adata, **kw, random_state=0, copy=True)
    assert a is not None
    assert b is not None
    assert c is not None
    assert (a.obsp["distances"] != b.obsp["distances"]).nnz == 0
    assert (a.obsp["distances"] != c.obsp["distances"]).nnz == 0
    with pytest.raises(TypeError, match="at most one"):
        sc.pp.wnn(adata, **kw, rng=0, random_state=0)


def test_wnn_transformer() -> None:
    """Approximate and custom kNN searches give (nearly) the exact result."""
    adata = _multimodal_adata()
    keys = ["X_pca", "X_apca"]
    exact = sc.pp.wnn(adata, keys, knn_range=50, copy=True)
    approx = sc.pp.wnn(
        adata, keys, knn_range=50, transformer="pynndescent", rng=0, copy=True
    )
    custom = sc.pp.wnn(
        adata,
        keys,
        knn_range=50,
        transformer=KNeighborsTransformer(n_neighbors=60, algorithm="kd_tree"),
        copy=True,
    )
    assert exact is not None
    assert approx is not None
    assert custom is not None
    assert (custom.obsp["distances"] != exact.obsp["distances"]).nnz == 0
    same = approx.obsp["distances"].indices == exact.obsp["distances"].indices
    assert same.mean() > 0.95


@needs.mudata
def test_wnn_mudata() -> None:
    import mudata as md

    adata = _multimodal_adata()
    mods = {
        name: AnnData(
            adata.X,
            obs=adata.obs,
            var=pd.DataFrame(index=[name]),
            obsm=dict(X_pca=adata.obsm[key]),
        )
        for name, key in (("rna", "X_pca"), ("prot", "X_apca"))
    }
    mdata = md.MuData(mods)
    with pytest.raises(TypeError, match="mapping"):
        sc.pp.wnn(mdata, ["X_pca", "X_pca"])
    sc.pp.wnn(mdata, {"rna": "X_pca", "prot": "X_pca"}, n_pcs=[12, None], knn_range=50)
    sc.pp.wnn(adata, ["X_pca", "X_apca"], knn_range=50)
    np.testing.assert_allclose(mdata.obs["rna_weight"], adata.obs["X_pca_weight"])
    assert (mdata.obsp["distances"] != adata.obsp["distances"]).nnz == 0
    params = mdata.uns["neighbors"]["params"]
    assert params["modalities"] == {
        "rna": {"use_rep": "X_pca", "n_pcs": 12},
        "prot": {"use_rep": "X_pca", "n_pcs": 8},
    }
    # `tl.umap` can’t read per-modality representations
    assert "use_rep" not in params
    sc.tl.umap(mdata, maxiter=10)
    assert np.isfinite(mdata.obsm["X_umap"]).all()


def test_wnn_duplicated_cells() -> None:
    # the cell moves to the front of its neighbors, or replaces the last one
    knn = np.array([[0, 1, 2], [0, 1, 2], [0, 1, 3], [1, 2, 0]], dtype=np.int32)
    expected = [[0, 1, 2], [1, 0, 2], [2, 0, 1], [3, 1, 2]]
    np.testing.assert_array_equal(_wnn._ensure_self_first(knn), expected)
    adata = _multimodal_adata(n=200)
    for key in ("X_pca", "X_apca"):
        adata.obsm[key][100:115] = adata.obsm[key][100]
    sc.pp.wnn(adata, ["X_pca", "X_apca"], knn_range=50)
    assert np.isfinite(adata.obs[["X_pca_weight", "X_apca_weight"]]).all().all()
    d = adata.obsp["distances"]
    assert np.all(d.indices.reshape(adata.n_obs, -1) != np.arange(adata.n_obs)[:, None])
    # many zero distances make `gauss`’ kernel widths zero
    sc.pp.wnn(adata, ["X_pca", "X_apca"], knn_range=50, method="gauss")
    assert np.isfinite(adata.obsp["connectivities"].data).all()
    # isolated cells with the same (e.g. empty) embedding: zero WNN kernel widths
    adata.obsm["X_apca"][:30] = 0
    with pytest.warns(UserWarning, match=r"\(30 in 'X_apca'\).*unreliable"):
        sc.pp.wnn(adata, ["X_pca", "X_apca"], knn_range=50)


def test_wnn_distances() -> None:
    """Distances keep Seurat’s order, and are exactly 0 for every modality’s nearest neighbor."""
    adata = _multimodal_adata()
    keys = ["X_pca", "X_apca"]
    umap, gauss = (
        sc.pp.wnn(adata, keys, knn_range=50, method=m, copy=True)
        for m in ("umap", "gauss")
    )
    assert umap is not None
    assert gauss is not None
    d = umap.obsp["distances"]
    assert (gauss.obsp["distances"] != d).nnz == 0
    np.testing.assert_array_equal(gauss.obsp["distances"].indices, d.indices)
    assert np.all(np.diff(d.data.reshape(adata.n_obs, -1), axis=1) >= 0)
    # no rounding error for the nearest neighbor in every modality
    nn1 = [
        NearestNeighbors(n_neighbors=2)
        .fit(x := _wnn._l2_normalize(adata.obsm[key].astype(np.float64)))
        .kneighbors(x, return_distance=False)[:, 1]
        for key in keys
    ]
    first = d.indices.reshape(adata.n_obs, -1)[:, 0]
    first_dist = d.data.reshape(adata.n_obs, -1)[:, 0]
    in_all = (first == nn1[0]) & (first == nn1[1])
    assert in_all.any()
    np.testing.assert_array_equal(first_dist[in_all], 0)
    assert (first_dist[~in_all] > 0).all()


@pytest.mark.filterwarnings("ignore:You’re trying to run this on:UserWarning")
def test_wnn_copy_default_rep() -> None:
    """The default representation is computed on the copy."""
    rng = np.random.default_rng(0)
    adata = AnnData(rng.normal(size=(200, 60)).astype(np.float32))
    adata.obsm["X_apca"] = rng.normal(size=(200, 5))
    out = sc.pp.wnn(adata, {"rna": None, "adt": "X_apca"}, knn_range=50, copy=True)
    assert out is not None
    assert "X_pca" in out.obsm
    assert "X_pca" not in adata.obsm
    assert "pca" not in adata.uns


def test_wnn_invalid_neighbors() -> None:
    class Broken(KNeighborsTransformer):
        def fit_transform(self, x, y=None):
            d = super().fit_transform(x, y)
            d.indices[d.indptr[1]] = -1  # e.g. IVF-style padding
            return d

    adata = _multimodal_adata()
    with pytest.raises(ValueError, match="invalid neighbors"):
        sc.pp.wnn(
            adata,
            ["X_pca", "X_apca"],
            knn_range=50,
            transformer=Broken(n_neighbors=50, algorithm="brute"),
        )


def test_wnn_too_few_neighbors() -> None:
    adata = _multimodal_adata()
    transformer = KNeighborsTransformer(n_neighbors=10)
    with pytest.raises(ValueError, match="too few neighbors"):
        sc.pp.wnn(adata, ["X_pca", "X_apca"], knn_range=50, transformer=transformer)
    assert not adata.obsp
    assert "neighbors" not in adata.uns
    view = adata[:300]
    with pytest.raises(ValueError, match="too few neighbors"):
        sc.pp.wnn(view, ["X_pca", "X_apca"], knn_range=50, transformer=transformer)
    assert view.is_view


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param({"use_rep": ["X_pca"]}, "two modalities", id="one_mod"),
        pytest.param({"use_rep": "X_pca"}, "one representation per", id="str"),
        pytest.param({"knn_range": 10, "n_neighbors": 20}, "larger than", id="range"),
        pytest.param({"knn_range": 1000}, "number of cells", id="range_n_obs"),
        pytest.param({"n_neighbors": 1}, "at least 2", id="k"),
        pytest.param({"n_pcs": [5]}, "one entry per modality", id="n_pcs"),
        pytest.param({"cross_constant": [1e-4]}, "one entry per", id="cross"),
        pytest.param({"method": "foo"}, "`method` needs", id="method"),
        pytest.param({"use_rep": ["X_pca", "X_pca"]}, "unique strings", id="dup"),
        pytest.param({"weight_key": "weight"}, "distinct column", id="weight_key"),
        pytest.param({"sd_scale": 0}, "positive", id="sd_scale"),
        pytest.param({"prune_snn": 1.5}, "between 0 and 1", id="prune"),
        pytest.param({"cross_constant": 0}, "positive", id="cross_zero"),
        pytest.param({"transformer": "foo"}, "`transformer` needs", id="transformer"),
        pytest.param(
            {"transformer": KNeighborsTransformer(n_neighbors=60, metric="cosine")},
            "euclidean",
            id="metric",
        ),
        pytest.param(
            {"use_rep": "acc"},
            "one representation per",
            id="acc",
            marks=needs.anndata_acc,
        ),
    ],
)
def test_wnn_errors(kwargs: dict, match: str) -> None:
    adata = _multimodal_adata()
    kwargs = {"use_rep": ["X_pca", "X_apca"], "knn_range": 50} | kwargs
    if kwargs["use_rep"] == "acc":
        from anndata.acc import A

        kwargs["use_rep"] = A.obsm["X_pca"]
    with pytest.raises((ValueError, TypeError), match=match):
        sc.pp.wnn(adata, **kwargs)
    assert not adata.obsp
    assert "neighbors" not in adata.uns
    view = adata[:300]
    with pytest.raises((ValueError, TypeError), match=match):
        sc.pp.wnn(view, **kwargs)
    assert view.is_view

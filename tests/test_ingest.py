from __future__ import annotations

from inspect import signature
from typing import TYPE_CHECKING

import anndata
import numpy as np
import pytest
from scipy import sparse
from sklearn.neighbors import KDTree
from umap import UMAP

import scanpy as sc
from scanpy import settings
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING or hasattr(anndata, "acc"):
    from anndata.acc import A

X = np.array(
    [
        [1.0, 2.5, 3.0, 5.0, 8.7],
        [4.2, 7.0, 9.0, 11.0, 7.0],
        [5.1, 2.0, 9.0, 4.0, 9.0],
        [7.0, 9.4, 6.8, 9.1, 8.0],
        [8.9, 8.6, 9.6, 1.0, 2.0],
        [6.5, 8.9, 2.2, 4.5, 8.9],
    ],
    dtype=np.float32,
)

T = np.array([[2.0, 3.5, 4.0, 1.0, 4.7], [3.2, 2.0, 5.0, 5.0, 8.0]], dtype=np.float32)


@pytest.fixture
def adatas():
    pbmc = pbmc68k_reduced()
    n_split = 500
    adata_ref = sc.AnnData(pbmc.X[:n_split, :], obs=pbmc.obs.iloc[:n_split])
    adata_new = sc.AnnData(pbmc.X[n_split:, :])

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)

    return adata_ref, adata_new


def test_representation(adatas):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)

    assert ing._use_rep == "X_pca"
    assert ing._obsm["rep"].shape == (adata_new.n_obs, settings.N_PCS)
    assert ing._pca_centered

    sc.pp.pca(adata_ref, n_comps=30, zero_center=False)
    sc.pp.neighbors(adata_ref)

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)

    assert ing._use_rep == "X_pca"
    assert ing._obsm["rep"].shape == (adata_new.n_obs, 30)
    assert not ing._pca_centered

    with pytest.warns(FutureWarning, match=r"argument use_rep is deprecated"):
        sc.pp.neighbors(adata_ref, use_rep="X")

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)

    assert ing._use_rep == "X"
    assert ing._obsm["rep"] is adata_new.X


@needs.anndata_acc
def test_representation_acc(adatas) -> None:
    """`use` is stored in `.uns` like `use_rep`, and used for the new data."""
    adata_ref, adata_new = (a.copy() for a in adatas)
    adata_new.obsm["X_pca"] = adata_ref.obsm["X_pca"][: adata_new.n_obs]
    sc.pp.neighbors(adata_ref, use=A.obsm["X_pca"])

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)

    assert ing._use_rep == "X_pca"
    np.testing.assert_array_equal(ing._obsm["rep"], adata_new.obsm["X_pca"])


@needs.anndata_acc
@pytest.mark.parametrize("as_sparse", [False, True])
def test_pca_transform_uses_reference_mean(
    as_sparse, monkeypatch: pytest.MonkeyPatch
) -> None:
    mask = np.tile([True, False, True, True, True], 12)
    ref_values = np.tile(X, (1, 12))
    query_values = np.tile(T, (1, 12))
    ref_values[:, ~mask] = np.nan
    query_values[:, ~mask] = np.inf
    if as_sparse:
        ref_x = sparse.csr_matrix(ref_values)  # noqa: TID251
        query_x = sparse.csr_matrix(query_values)  # noqa: TID251
    else:
        ref_x = ref_values
        query_x = query_values
    adata_ref = sc.AnnData(ref_x)
    adata_new = sc.AnnData(query_x)
    adata_ref.var["selected"] = mask
    sc.pp.pca(adata_ref, n_comps=3, mask="var.selected")
    sc.pp.neighbors(adata_ref, n_neighbors=3, n_pcs=2)
    ing = sc.tl.Ingest(adata_ref)

    if as_sparse:

        def fail_toarray(*args, **kwargs) -> None:
            pytest.fail("PCA ingest must not densify the query matrix")

        monkeypatch.setattr(type(adata_new.X), "toarray", fail_toarray)

    expected = (query_values[:, mask] - ref_values[:, mask].mean(axis=0)) @ (
        adata_ref.varm["PCs"][mask]
    )
    ing.fit(adata_new)
    np.testing.assert_allclose(ing._obsm["rep"], expected[:, :2], rtol=1e-5, atol=1e-5)
    ing.map_embedding("pca")
    np.testing.assert_allclose(ing._obsm["X_pca"], expected, rtol=1e-5, atol=1e-5)

    individual = []
    for i in range(adata_new.n_obs):
        ing.fit(adata_new[[i]].copy())
        ing.map_embedding("pca")
        individual.append(ing._obsm["X_pca"])
    individual = np.vstack(individual)
    np.testing.assert_allclose(individual, expected, rtol=1e-5, atol=1e-5)


def test_neighbors(adatas):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)
    ing.neighbors(k=10)
    indices = ing._indices

    tree = KDTree(adata_ref.obsm["X_pca"])
    true_indices = tree.query(ing._obsm["rep"], 10, return_distance=False)

    num_correct = 0.0
    for i in range(adata_new.n_obs):
        num_correct += np.sum(np.isin(true_indices[i], indices[i]))
    percent_correct = num_correct / (adata_new.n_obs * 10)

    assert percent_correct > 0.99


@pytest.mark.parametrize("n", [8, 12])
def test_neighbors_defaults(adatas, n):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    sc.pp.neighbors(adata_ref, n_neighbors=n)

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)
    ing.neighbors()
    assert ing._indices.shape[1] == n


# https://github.com/lmcinnes/umap/issues/1174
@pytest.mark.filterwarnings("ignore:.*renamed to.*ensure_all_finite:FutureWarning")
def test_ingest_function(adatas: tuple[sc.AnnData, sc.AnnData]) -> None:
    adata_ref, adata_new = (ad.copy() for ad in adatas)

    sc.tl.ingest(
        adata_new,
        adata_ref,
        obs="bulk_labels",
        embedding_method=["umap", "pca"],
    )

    assert "bulk_labels" in adata_new.obs
    assert "X_umap" in adata_new.obsm
    assert "X_pca" in adata_new.obsm

    with pytest.warns(FutureWarning, match=r"argument inplace is deprecated"):
        ad = sc.tl.ingest(
            adata_new,
            adata_ref,
            obs="bulk_labels",
            embedding_method=["umap", "pca"],
            inplace=False,
        )

    assert "bulk_labels" in ad.obs
    assert "X_umap" in ad.obsm
    assert "X_pca" in ad.obsm


# https://github.com/lmcinnes/umap/issues/1174
@pytest.mark.filterwarnings("ignore:.*renamed to.*ensure_all_finite:FutureWarning")
@pytest.mark.filterwarnings("ignore:.*argument use_rep is deprecated:FutureWarning")
@pytest.mark.parametrize("n_components", [2, 3])
@pytest.mark.parametrize("as_sparse", [False, True])
@pytest.mark.parametrize("metric_kwds", [{}, {"p": 3}])
def test_ingest_map_embedding_umap(
    n_components, as_sparse, metric_kwds, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sc.settings, "n_jobs", 1)
    metric = "minkowski" if metric_kwds else "euclidean"
    ref_x = sparse.csr_matrix(X) if as_sparse else X  # noqa: TID251
    query_x = sparse.csr_matrix(T) if as_sparse else T  # noqa: TID251
    adata_ref = sc.AnnData(ref_x)
    adata_new = sc.AnnData(query_x)

    sc.pp.neighbors(
        adata_ref,
        method="umap",
        use_rep="X",
        n_neighbors=4,
        metric=metric,
        metric_kwds=metric_kwds,
        transformer="pynndescent",
        random_state=0,
    )
    sc.tl.umap(adata_ref, n_components=n_components, a=1.2, b=0.8, random_state=0)
    reference_embedding = adata_ref.obsm["X_umap"].copy()

    # Use the modern transform when available, retaining support for UMAP 0.5.
    compat = (
        dict(compatibility_layout=False)
        if "compatibility_layout" in signature(UMAP).parameters
        else {}
    )
    reducer = UMAP(
        n_components=n_components,
        metric=metric,
        metric_kwds=metric_kwds,
        a=1.2,
        b=0.8,
        random_state=0,
        n_neighbors=4,
        n_jobs=1,
        **compat,
    ).fit(ref_x)
    # Compare transforms into the same stored layout, independent of how it was fit.
    reducer.embedding_ = reference_embedding.copy()
    umap_transformed_t = reducer.transform(query_x)

    def fail_fit(*args, **kwargs):
        pytest.fail("Ingest must reuse the reference embedding without fitting UMAP")

    monkeypatch.setattr(UMAP, "fit", fail_fit)
    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)
    ing.map_embedding(method="umap")

    assert ing._obsm["X_umap"].shape == (adata_new.n_obs, n_components)
    assert np.isfinite(ing._obsm["X_umap"]).all()
    assert np.allclose(ing._obsm["X_umap"], umap_transformed_t)
    np.testing.assert_array_equal(adata_ref.obsm["X_umap"], reference_embedding)


@pytest.mark.parametrize(
    ("n_ref", "as_sparse", "use_rep"),
    [
        (48, False, "X"),
        (48, True, "X"),
        (4096, False, "X"),
        (4096, True, "X"),
        (48, False, "X_pca"),
        (48, False, "custom"),
    ],
)
def test_ingest_umap_reference(
    n_ref, as_sparse, use_rep, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sc.settings, "n_jobs", 1)
    rng = np.random.default_rng(42)
    n_vars = 60 if use_rep == "X_pca" else 8
    x = rng.normal(scale=0.1, size=(n_ref, n_vars)).astype(np.float32)
    x[n_ref // 2 :] += 10
    query = x[[0, -1]] + 0.01
    adata_ref = sc.AnnData(sparse.csr_matrix(x) if as_sparse else x)  # noqa: TID251
    adata_new = sc.AnnData(sparse.csr_matrix(query) if as_sparse else query)  # noqa: TID251
    adata_ref.obs["label"] = np.repeat(["a", "b"], n_ref // 2)
    adata_ref.obs["label"] = adata_ref.obs["label"].astype("category")
    expected_rep = query
    if use_rep == "X_pca":
        sc.pp.pca(adata_ref, n_comps=3)
        expected_rep = (query - x.mean(axis=0)) @ adata_ref.varm["PCs"]
    elif use_rep == "custom":
        adata_ref.obsm["custom"] = x * 2
        adata_new.obsm["custom"] = expected_rep = query * 2
    rep_args = dict(n_pcs=3) if use_rep == "X_pca" else dict(use_rep=use_rep)
    with pytest.warns(
        FutureWarning, match=r"argument (use_rep|key_added) is deprecated"
    ):
        sc.pp.neighbors(adata_ref, **rep_args, n_neighbors=4, key_added="reference")
    with pytest.warns(FutureWarning, match=r"argument neighbors_key is deprecated"):
        sc.tl.umap(adata_ref, neighbors_key="reference", n_components=3, maxiter=20)
    reference_embedding = adata_ref.obsm["X_umap"].copy()

    sc.tl.ingest(
        adata_new,
        adata_ref,
        obs="label",
        embedding_method="umap",
        neighbors_key="reference",
    )

    assert adata_new.obsm["X_umap"].shape == (2, 3)
    assert np.isfinite(adata_new.obsm["X_umap"]).all()
    actual_rep = adata_new.obsm["rep"]
    if as_sparse:
        actual_rep = actual_rep.toarray()
    np.testing.assert_allclose(actual_rep, expected_rep, rtol=1e-6, atol=1e-5)
    assert adata_new.obs["label"].tolist() == ["a", "b"]
    assert adata_new.obs["label"].dtype == adata_ref.obs["label"].dtype
    np.testing.assert_array_equal(adata_ref.obsm["X_umap"], reference_embedding)


def test_ingest_backed(adatas, tmp_path):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    adata_new.write_h5ad(f"{tmp_path}/new.h5ad")

    adata_new = anndata.read_h5ad(f"{tmp_path}/new.h5ad", backed="r")

    ing = sc.tl.Ingest(adata_ref)
    with pytest.raises(
        NotImplementedError,
        match=f"Ingest.fit is not implemented for matrices of type {type(adata_new.X)}",
    ):
        ing.fit(adata_new)

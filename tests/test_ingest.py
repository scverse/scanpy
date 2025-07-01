from __future__ import annotations

import anndata
import numpy as np
import pytest
from sklearn.neighbors import KDTree
from umap import UMAP

import scanpy as sc
from scanpy import settings
from scanpy._compat import pkg_version
from testing.scanpy._helpers.data import pbmc68k_reduced

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

    sc.pp.neighbors(adata_ref, use_rep="X")

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)

    assert ing._use_rep == "X"
    assert ing._obsm["rep"] is adata_new.X


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


@pytest.mark.parametrize("n", [3, 4])
def test_neighbors_defaults(adatas, n):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    sc.pp.neighbors(adata_ref, n_neighbors=n)

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)
    ing.neighbors()
    assert ing._indices.shape[1] == n


@pytest.mark.skipif(
    pkg_version("anndata") < sc.tl._ingest.ANNDATA_MIN_VERSION,
    reason="`AnnData.concatenate` does not concatenate `.obsm` in old anndata versions",
)
def test_ingest_function(adatas):
    adata_ref = adatas[0].copy()
    adata_new = adatas[1].copy()

    sc.tl.ingest(
        adata_new,
        adata_ref,
        obs="bulk_labels",
        embedding_method=["umap", "pca"],
        inplace=True,
    )

    assert "bulk_labels" in adata_new.obs
    assert "X_umap" in adata_new.obsm
    assert "X_pca" in adata_new.obsm

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


def test_ingest_map_embedding_umap():
    adata_ref = sc.AnnData(X)
    adata_new = sc.AnnData(T)

    sc.pp.neighbors(
        adata_ref, method="umap", use_rep="X", n_neighbors=4, random_state=0
    )
    sc.tl.umap(adata_ref, random_state=0)

    ing = sc.tl.Ingest(adata_ref)
    ing.fit(adata_new)
    ing.map_embedding(method="umap")

    reducer = UMAP(min_dist=0.5, random_state=0, n_neighbors=4)
    reducer.fit(X)
    umap_transformed_t = reducer.transform(T)

    assert np.allclose(ing._obsm["X_umap"], umap_transformed_t)


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

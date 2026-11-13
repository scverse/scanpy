from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from anndata.tests.helpers import assert_equal
from sklearn.metrics import silhouette_score

import scanpy as sc
from scanpy.preprocessing._combat import _design_matrix, _standardize_data
from testing.scanpy._pytest.marks import needs


def test_norm():
    # this test trivially checks whether mean normalisation worked

    # load in data
    adata = sc.datasets.blobs()
    key = "blobs"
    data = pd.DataFrame(data=adata.X.T, index=adata.var_names, columns=adata.obs_names)

    # construct a pandas series of the batch annotation
    batch = pd.Series(adata.obs[key])
    model = pd.DataFrame({"batch": batch})

    # standardize the data
    s_data, _design, _var_pooled, _stand_mean = _standardize_data(model, data, "batch")

    assert np.allclose(s_data.mean(axis=1), np.zeros(s_data.shape[0]))


def test_covariates():
    rng = np.random.default_rng()
    adata = sc.datasets.blobs()
    key = "blobs"

    x1 = sc.pp.combat(adata, key, out=True)

    adata.obs["cat1"] = rng.binomial(3, 0.5, size=adata.n_obs)
    adata.obs["cat2"] = rng.binomial(2, 0.1, size=adata.n_obs)
    adata.obs["num1"] = rng.standard_normal(adata.n_obs)

    x2 = sc.pp.combat(adata, key, covariates=["cat1", "cat2", "num1"], out=True)
    sc.pp.combat(adata, key, covariates=["cat1", "cat2", "num1"])

    assert x1.shape == x2.shape

    df = adata.obs[["cat1", "cat2", "num1", key]]
    batch_cats = adata.obs[key].cat.categories
    design = _design_matrix(df, key, batch_cats)

    assert len(design.columns) == 4 + len(batch_cats) - 1


@pytest.mark.parametrize("covariates", [["cond"], ["cond", "donor"]])
def test_categorical_covariates(covariates: list[str]) -> None:
    """Categorical covariates are modeled, so their effect survives batch correction."""
    rng = np.random.default_rng(0)
    batch = rng.integers(0, 3, 600)
    treated = rng.random(600) < np.where(batch == 0, 0.8, 0.2)
    cond_effect = rng.normal(0, 2, 10)
    x = (
        5
        + rng.normal(0, 2, (3, 10))[batch]
        + np.outer(treated, cond_effect)
        + rng.normal(0, 0.5, (600, 10))
    )
    adata = AnnData(x)
    adata.obs = pd.DataFrame(
        dict(
            batch=batch.astype(str),
            cond=np.where(treated, "treated", "ctrl"),
            donor=rng.choice(["d0", "d1", "d2"], 600),
        ),
        index=adata.obs_names,
    ).astype("category")

    corrected = sc.pp.combat(adata, "batch", covariates=covariates, out=True)

    np.testing.assert_allclose(
        corrected[treated].mean(0) - corrected[~treated].mean(0), cond_effect, atol=0.2
    )


def test_combat_obs_names():
    """Regression test for <https://github.com/scverse/scanpy/issues/1170>."""
    rng = np.random.default_rng()
    x = rng.random((200, 100))
    obs = pd.DataFrame(
        {"batch": pd.Categorical(rng.integers(0, 2, (200,)))},
        index=np.repeat(np.arange(100), 2).astype(str),  # Non-unique index
    )
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        a = sc.AnnData(x, obs)
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        b = a.copy()
    b.obs_names_make_unique()

    sc.pp.combat(a, "batch")
    sc.pp.combat(b, "batch")

    assert_equal(a.X, b.X)

    a.obs_names_make_unique()
    assert_equal(a, b)


def test_combat_single_cell_batch():
    """Test that combat raises an error when a batch has fewer than 2 cells.

    Regression test for https://github.com/scverse/scanpy/issues/1175
    """
    adata = sc.datasets.blobs()
    # Create a batch where one category has only 1 cell
    batch = pd.Categorical(["single"] + ["other"] * (adata.n_obs - 1))
    adata.obs["batch"] = batch

    with pytest.raises(ValueError, match="fewer than 2 cells"):
        sc.pp.combat(adata, batch="batch")


def test_silhouette():
    # this test checks wether combat can align data from several gaussians
    # it checks this by computing the silhouette coefficient in a pca embedding

    # load in data
    adata = sc.datasets.blobs()

    # apply combat
    sc.pp.combat(adata, "blobs")

    # compute pca
    sc.pp.pca(adata)
    x_pca = adata.obsm["X_pca"]

    # compute silhouette coefficient in pca
    sh = silhouette_score(x_pca[:, :2], adata.obs["blobs"].values)

    assert sh < 0.1


@needs.anndata_acc
def test_use_out() -> None:
    """`use`/`out` pick where to read and write, leaving `X` alone."""
    from anndata.acc import A

    adata = sc.datasets.blobs()
    adata.layers["counts"] = adata.X.copy()
    expected = sc.pp.combat(adata, "blobs", out=True)
    sc.pp.combat(adata, "blobs", use=A.layers["counts"], out=A.layers["p"])
    np.testing.assert_array_equal(adata.X, adata.layers["counts"])
    np.testing.assert_allclose(adata.layers["p"], expected)


def test_inplace_deprecated() -> None:
    adata = sc.datasets.blobs()
    with pytest.warns(FutureWarning, match=r"argument inplace is deprecated"):
        x = sc.pp.combat(adata, "blobs", inplace=False)
    np.testing.assert_allclose(x, sc.pp.combat(adata, "blobs", out=True))


@needs.anndata_acc
def test_key_ref() -> None:
    """`key` and `covariates` can also be references."""
    from anndata.acc import A

    adata = sc.datasets.blobs()
    adata.obs["cov"] = np.random.default_rng(0).normal(size=adata.n_obs)
    expected = sc.pp.combat(adata, "blobs", covariates=["cov"], out=True)
    res = sc.pp.combat(adata, A.obs["blobs"], covariates=[A.obs["cov"]], out=True)
    np.testing.assert_allclose(res, expected)


def test_key_deprecated() -> None:
    adata = sc.datasets.blobs()
    expected = sc.pp.combat(adata, "blobs", out=True)
    with pytest.warns(FutureWarning, match=r"argument key is deprecated"):
        res = sc.pp.combat(adata, key="blobs", out=True)
    np.testing.assert_allclose(res, expected)

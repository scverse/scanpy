from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix  # noqa: TID251

import scanpy as sc
from scanpy.preprocessing._combat import _design_matrix


def _categorical_example(n_covariates):
    combinations = list(product(["batch_a", "batch_b"], ["a", "b", "c"], ["x", "y"]))
    obs = pd.DataFrame(combinations * 4, columns=["batch", "condition", "donor"])
    obs.index = obs.index.astype(str)
    for key in obs:
        obs[key] = obs[key].astype("category")
    rng = np.random.default_rng(4454)
    obs["age"] = rng.normal(size=len(obs))
    x = rng.normal(size=(len(obs), 24))
    x += obs["batch"].cat.codes.to_numpy()[:, None] * np.linspace(0.5, 2, 24)
    x += obs["condition"].cat.codes.to_numpy()[:, None] * np.linspace(0.3, 1, 24)
    categories = ["condition", "donor"][:n_covariates]
    numeric = pd.concat(
        [
            obs[["batch", "age"]],
            pd.get_dummies(obs[categories], drop_first=True, dtype=float),
        ],
        axis=1,
    )
    return AnnData(x, obs=obs), numeric, categories


@pytest.mark.parametrize("n_covariates", [1, 2])
def test_categorical_design_is_full_rank(n_covariates):
    adata, numeric, categories = _categorical_example(n_covariates)
    model = adata.obs[["batch", *categories, "age"]]
    levels = adata.obs["batch"].cat.categories
    design = _design_matrix(model, "batch", levels)
    expected = _design_matrix(numeric, "batch", levels)
    assert design.shape == expected.shape
    assert np.linalg.matrix_rank(design) == design.shape[1]
    expected_columns = [*expected.columns[:2], *numeric.columns[2:], "age"]
    np.testing.assert_allclose(design.to_numpy(), expected[expected_columns].to_numpy())


@pytest.mark.parametrize("n_covariates", [1, 2])
@pytest.mark.parametrize("sparse", [False, True])
def test_categorical_combat_matches_reference_coding(n_covariates, sparse):
    adata, numeric, categories = _categorical_example(n_covariates)
    x = adata.X.copy()
    if sparse:
        adata.X = csr_matrix(adata.X)
    reference = AnnData(x.copy(), obs=numeric)
    expected = sc.pp.combat(
        reference, key="batch", covariates=list(numeric.columns[1:]), inplace=False
    )
    actual = sc.pp.combat(
        adata, key="batch", covariates=[*categories, "age"], inplace=False
    )
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)
    np.testing.assert_array_equal(adata.X.toarray() if sparse else adata.X, x)
    sc.pp.combat(adata, key="batch", covariates=[*categories, "age"])
    np.testing.assert_allclose(adata.X, expected, rtol=1e-10, atol=1e-10)


def test_single_level_categorical_covariate_matches_batch_only():
    adata, _, _ = _categorical_example(1)
    adata.obs["constant"] = pd.Categorical(["same"] * adata.n_obs)
    expected = sc.pp.combat(adata, key="batch", inplace=False)
    actual = sc.pp.combat(adata, key="batch", covariates=["constant"], inplace=False)
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)

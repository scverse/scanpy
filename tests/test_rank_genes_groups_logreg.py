from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scanpy as sc


@pytest.mark.filterwarnings("ignore:invalid value encountered in log2:RuntimeWarning")
@pytest.mark.parametrize("method", ["t-test", "logreg"])
def test_rank_genes_groups_with_renamed_categories(method):
    adata = sc.datasets.blobs(n_variables=4, n_centers=3, n_observations=200)
    assert np.allclose(adata.X[1], [9.214668, -2.6487126, 4.2020774, 0.51076424])

    # for method in ['logreg', 't-test']:

    sc.tl.rank_genes_groups(adata, "blobs", method=method)
    assert adata.uns["rank_genes_groups"]["names"].dtype.names == ("0", "1", "2")
    assert adata.uns["rank_genes_groups"]["names"][0].tolist() == ("1", "3", "0")

    adata.rename_categories("blobs", ["Zero", "One", "Two"])
    assert adata.uns["rank_genes_groups"]["names"][0].tolist() == ("1", "3", "0")

    sc.tl.rank_genes_groups(adata, "blobs", method=method)
    assert adata.uns["rank_genes_groups"]["names"][0].tolist() == ("1", "3", "0")
    assert adata.uns["rank_genes_groups"]["names"].dtype.names == ("Zero", "One", "Two")


def test_rank_genes_groups_with_renamed_categories_use_rep():
    adata = sc.datasets.blobs(n_variables=4, n_centers=3, n_observations=200)
    assert np.allclose(adata.X[1], [9.214668, -2.6487126, 4.2020774, 0.51076424])

    adata.layers["to_test"] = adata.X.copy()
    adata.X = adata.X[::-1, :]

    sc.tl.rank_genes_groups(
        adata, "blobs", method="logreg", layer="to_test", use_raw=False
    )
    assert adata.uns["rank_genes_groups"]["names"].dtype.names == ("0", "1", "2")
    assert adata.uns["rank_genes_groups"]["names"][0].tolist() == ("1", "3", "0")

    sc.tl.rank_genes_groups(adata, "blobs", method="logreg")
    assert adata.uns["rank_genes_groups"]["names"][0].tolist() != ("3", "1", "0")


def test_rank_genes_groups_with_unsorted_groups():
    adata = sc.datasets.blobs(n_variables=10, n_centers=5, n_observations=200)
    adata._sanitize()
    adata.rename_categories("blobs", ["Zero", "One", "Two", "Three", "Four"])
    bdata = adata.copy()
    sc.tl.rank_genes_groups(
        adata, "blobs", groups=["Zero", "One", "Three"], method="logreg"
    )
    sc.tl.rank_genes_groups(
        bdata, "blobs", groups=["One", "Three", "Zero"], method="logreg"
    )
    array_ad = pd.DataFrame(
        adata.uns["rank_genes_groups"]["scores"]["Three"]
    ).to_numpy()
    array_bd = pd.DataFrame(
        bdata.uns["rank_genes_groups"]["scores"]["Three"]
    ).to_numpy()
    np.testing.assert_equal(array_ad, array_bd)


@pytest.mark.parametrize("categories", [["A", "unused", "B"], ["B", "unused", "A"]])
@pytest.mark.parametrize("target", ["A", "B"])
@pytest.mark.parametrize("representation", ["dense", "sparse", "raw", "layer"])
def test_binary_logreg_scores_point_toward_requested_group(
    categories, target, representation
):
    from anndata import AnnData
    from scipy.sparse import csr_matrix  # noqa: TID251

    x = np.vstack([
        np.tile([10.0, 0.0, 1.0], (30, 1)),
        np.tile([0.0, 10.0, 1.0], (30, 1)),
    ])
    adata = AnnData(
        x,
        obs=pd.DataFrame(
            {"group": pd.Categorical(["A"] * 30 + ["B"] * 30, categories=categories)},
            index=[f"cell_{i}" for i in range(60)],
        ),
        var=pd.DataFrame(index=["A_marker", "B_marker", "shared"]),
    )
    kwargs = {"use_raw": False}
    if representation == "sparse":
        adata.X = csr_matrix(adata.X)
    elif representation == "raw":
        adata.raw = adata.copy()
        adata.X = np.zeros_like(x)
        kwargs = {"use_raw": True}
    elif representation == "layer":
        adata.layers["expression"] = adata.X.copy()
        adata.X = np.zeros_like(x)
        kwargs["layer"] = "expression"
    reference = "B" if target == "A" else "A"
    sc.tl.rank_genes_groups(
        adata, "group", groups=[target], reference=reference, method="logreg", **kwargs
    )
    result = adata.uns["rank_genes_groups"]
    names = result["names"][target]
    scores = dict(zip(names, result["scores"][target], strict=True))
    assert names[0] == f"{target}_marker"
    assert scores[f"{target}_marker"] > 0
    assert scores[f"{reference}_marker"] < 0


@pytest.mark.parametrize("categories", [["A", "B"], ["B", "A"]])
def test_binary_logreg_default_group_direction(categories):
    from anndata import AnnData

    adata = AnnData(
        np.vstack([np.tile([10.0, 0.0], (20, 1)), np.tile([0.0, 10.0], (20, 1))]),
        obs=pd.DataFrame(
            {"group": pd.Categorical(["A"] * 20 + ["B"] * 20, categories=categories)},
            index=[f"cell_{i}" for i in range(40)],
        ),
        var=pd.DataFrame(index=["A_marker", "B_marker"]),
    )
    sc.tl.rank_genes_groups(adata, "group", method="logreg", use_raw=False)
    result = adata.uns["rank_genes_groups"]
    target = categories[0]
    assert result["names"][target][0] == f"{target}_marker"
    assert result["scores"][target][0] > 0

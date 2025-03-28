from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scanpy as sc


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

from __future__ import annotations

import pytest

import scanpy as sc
from scanpy.testing._helpers.data import pbmc68k_reduced

n_neighbors = 5
key = "test"


@pytest.mark.parametrize(
    ("groupby", "dendrogram_key"),
    [
        pytest.param(None, "dendrogram_obs", id="obs"),
        pytest.param("bulk_labels", "dendrogram_bulk_labels", id="bulk_labels"),
        pytest.param(
            ["bulk_labels", "phase"],
            "dendrogram_bulk_labels_phase",
            id="bulk_labels+phase",
        ),
    ],
)
@pytest.mark.parametrize(
    "key_added", [None, "custom_key"], ids=["derived", "custom_key"]
)
def test_dendrogram_key_added(groupby, dendrogram_key, key_added):
    adata = pbmc68k_reduced()
    sc.tl.dendrogram(adata, groupby=groupby, key_added=key_added, use_rep="X_pca")

    if key_added is None:
        key_added = dendrogram_key
    assert key_added in adata.uns


def test_dendrogram_var_names():
    adata = pbmc68k_reduced()
    sc.tl.dendrogram(adata, var_names=adata.var_names[:10], axis=1)
    assert len(adata.uns["dendrogram_var"]["categories_ordered"]) == 10

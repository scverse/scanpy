from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.testing import assert_array_equal

import scanpy as sc


@pytest.fixture
def adata():
    return AnnData(
        np.arange(12, dtype=float).reshape(4, 3),
        obs=pd.DataFrame(
            {"group": pd.Categorical(["a", "b", "a", "b"]), "score": [3, 1, 4, 2]},
            index=["a", "b", "c", "d"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


def _coordinates(ax):
    ax.figure.canvas.draw()
    collection = ax.collections[0]
    if ax.name == "3d":
        return np.column_stack(collection._offsets3d)
    return collection.get_offsets()


@pytest.mark.parametrize(
    ("basis", "keys", "expected_key", "label"),
    [
        ("umap", ["umap"], "umap", "UMAP"),
        ("X_umap", ["X_umap", "X_X_umap"], "X_umap", "UMAP"),
        ("umap", ["X_umap"], "X_umap", "UMAP"),
        ("umap", ["umap", "X_umap"], "umap", "UMAP"),
        ("custom", ["custom", "X_custom"], "custom", "custom"),
        ("X_X_umap", ["X_X_umap", "X_umap"], "X_X_umap", "X_umap"),
        ("X_umap", ["X_X_umap"], "X_X_umap", "UMAP"),
    ],
)
def test_scatter_basis_lookup(adata, basis, keys, expected_key, label):
    for i, key in enumerate(keys):
        adata.obsm[key] = np.arange(16).reshape(4, 4) + 100 * i

    ax = sc.pl.scatter(adata, basis=basis, show=False)

    assert_array_equal(_coordinates(ax), adata.obsm[expected_key][:, :2])
    assert (ax.get_xlabel(), ax.get_ylabel()) == (f"{label}1", f"{label}2")


@pytest.mark.parametrize(
    ("projection", "components", "columns"),
    [("2d", "3,1", [2, 0]), ("3d", ["4", "2", "1"], [3, 1, 0])],
)
def test_scatter_basis_components(adata, projection, components, columns):
    adata.obsm["X_custom"] = np.arange(16).reshape(4, 4)

    ax = sc.pl.scatter(
        adata,
        basis="X_custom",
        projection=projection,
        components=components,
        show=False,
    )

    assert_array_equal(_coordinates(ax), adata.obsm["X_custom"][:, columns])
    labels = [ax.get_xlabel(), ax.get_ylabel()]
    if projection == "3d":
        labels.append(ax.get_zlabel())
    assert labels == [f"custom{c + 1}" for c in columns]


@pytest.mark.parametrize(
    ("basis", "keys", "projection", "components", "columns"),
    [
        ("diffmap", ["X_diffmap"], "2d", None, [1, 2]),
        ("diffmap", ["diffmap", "X_diffmap"], "2d", "3,1", [3, 1]),
        ("X_diffmap", ["X_diffmap"], "3d", None, [1, 2, 3]),
        ("X_diffmap", ["X_diffmap"], "3d", "3,1,2", [3, 1, 2]),
    ],
)
def test_scatter_basis_diffmap(adata, *, basis, keys, projection, components, columns):
    for i, key in enumerate(keys):
        adata.obsm[key] = np.arange(16).reshape(4, 4) + 100 * i

    ax = sc.pl.scatter(
        adata, basis=basis, projection=projection, components=components, show=False
    )

    assert_array_equal(_coordinates(ax), adata.obsm[keys[0]][:, columns])
    labels = [ax.get_xlabel(), ax.get_ylabel()]
    if projection == "3d":
        labels.append(ax.get_zlabel())
    assert labels == [f"DC{c}" for c in columns]


def test_scatter_basis_missing_prefix(adata):
    adata.obsm["umap"] = np.arange(8).reshape(4, 2)
    with pytest.raises(
        KeyError, match=r"Could not find 'X_umap' or 'X_X_umap' in \.obsm"
    ):
        sc.pl.scatter(adata, basis="X_umap", show=False)


@pytest.mark.parametrize("source", ["raw", "X", "layer"])
def test_scatter_basis_expression(adata, source):
    adata.raw = adata.copy()
    adata.X = adata.X + 100
    adata.layers["counts"] = adata.X + 100
    adata.obsm["custom"] = np.arange(8).reshape(4, 2)
    kwargs = {"use_raw": source == "raw"}
    if source == "layer":
        kwargs["layers"] = "counts"
    expected = {"raw": adata.raw.X, "X": adata.X, "layer": adata.layers["counts"]}

    ax = sc.pl.scatter(adata, basis="custom", color="g1", show=False, **kwargs)

    assert_array_equal(_coordinates(ax), adata.obsm["custom"])
    assert_array_equal(ax.collections[0].get_array(), expected[source][:, 0])


def test_scatter_basis_categorical_axes(adata, plt):
    adata.obsm["custom"] = np.arange(8).reshape(4, 2)
    _, ax = plt.subplots()
    result = sc.pl.scatter(
        adata,
        basis="custom",
        color="group",
        groups="a",
        palette=["red", "blue"],
        title="Selected cells",
        size=17,
        alpha=0.5,
        frameon=False,
        ax=ax,
        show=False,
    )

    assert result is ax
    assert_array_equal(_coordinates(ax), adata.obsm["custom"][[0, 2]])
    assert_array_equal(ax.collections[1].get_offsets(), adata.obsm["custom"][[1, 3]])
    assert_array_equal(ax.collections[0].get_sizes(), [17])
    assert ax.collections[0].get_alpha() == 0.5
    assert ax.get_title() == "Selected cells"
    assert not ax.get_frame_on()
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["a"]


def test_scatter_xy_unchanged(adata):
    ax = sc.pl.scatter(adata, x="g1", y="g2", color="score", show=False)

    order = np.argsort(adata.obs["score"])
    assert_array_equal(_coordinates(ax), adata.X[order, :2])
    assert_array_equal(ax.collections[0].get_array(), adata.obs["score"].iloc[order])
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("g1", "g2")

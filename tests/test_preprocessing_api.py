from __future__ import annotations

import pytest

import scanpy as sc


@pytest.mark.parametrize("namespace", ["pp", "preprocessing"])
@pytest.mark.parametrize(
    "name", ["recipe_seurat", "recipe_zheng17", "recipe_weinreb17"]
)
def test_legacy_recipes_removed(namespace: str, name: str) -> None:
    module = getattr(sc, namespace)
    assert not hasattr(module, name)
    assert name not in module.__all__


def test_experimental_recipe_available() -> None:
    assert "recipe_pearson_residuals" in sc.experimental.pp.__all__
    assert callable(sc.experimental.pp.recipe_pearson_residuals)

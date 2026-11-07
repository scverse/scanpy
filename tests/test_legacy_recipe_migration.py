from __future__ import annotations

from contextlib import nullcontext
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from scanpy._compat import CSBase

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize("sparse_input", [False, True], ids=["dense", "csr"])
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_zheng_migration(project_dir: Path, dtype: str, *, sparse_input: bool) -> None:
    """Execute the guide with no accessor requirement and check each transformation."""
    rng = np.random.default_rng(0)
    counts = rng.poisson(rng.lognormal(size=(48, 1201))).astype(dtype)
    counts[:, -1] = 0  # Exercise filtering before normalization and HVG selection.
    adata = AnnData(sparse.csr_matrix(counts) if sparse_input else counts.copy())  # noqa: TID251

    totals = counts.sum(axis=1)
    factors = totals / np.median(totals)
    logged = np.log1p(counts[:, :-1] / factors[:, None])
    expected_hvg = sc.pp.highly_variable_genes(
        AnnData(logged), flavor="cell_ranger", n_top_genes=1000, out=True
    )["highly_variable"].to_numpy()
    assert 0 < expected_hvg.sum() < logged.shape[1]

    guide = project_dir / "docs/how-to/legacy-recipes.md"
    section = guide.read_text().split("## Zheng et al. (2017)", 1)[1]
    snippet = section.split("```python\n", 1)[1].split("```", 1)[0]
    before_scale, scale, after_scale = snippet.partition("sc.pp.scale(")
    assert scale
    namespace = {"adata": adata}
    exec(compile(before_scale, str(guide), "exec"), namespace)
    migrated = namespace["adata"]
    actual_logged = (
        migrated.X.toarray() if isinstance(migrated.X, CSBase) else migrated.X
    )
    np.testing.assert_array_equal(namespace["hvg"]["highly_variable"], expected_hvg)
    np.testing.assert_array_equal(migrated.var_names, adata.var_names[expected_hvg])
    np.testing.assert_allclose(migrated.obs["n_counts_all"], factors, rtol=1e-6)
    np.testing.assert_allclose(actual_logged, logged[:, expected_hvg], rtol=1e-6)
    assert migrated.uns["log1p"] == {"base": None}

    warn_ctx = (
        pytest.warns(UserWarning, match=r"zero-centering.*densifies")
        if sparse_input
        else nullcontext()
    )
    with warn_ctx:
        exec(compile(scale + after_scale, str(guide), "exec"), namespace)
    selected = logged[:, expected_hvg].astype(np.float64)
    std = selected.std(axis=0, ddof=1)
    std[std == 0] = 1
    expected_scaled = (selected - selected.mean(axis=0)) / std
    np.testing.assert_allclose(migrated.X, expected_scaled, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(migrated.obs["n_counts_all"], factors, rtol=1e-6)

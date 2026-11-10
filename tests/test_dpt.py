from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


@pytest.fixture
def pbmc():
    adata = pbmc68k_reduced()
    adata.uns["iroot"] = 0
    return adata


@needs.anndata_acc
def test_dpt_use(pbmc) -> None:
    """`use` finds a diffmap written with the same `out`, including its eigenvalues."""
    expected = pbmc.copy()
    sc.tl.diffmap(expected)
    sc.tl.dpt(expected)

    sc.tl.diffmap(pbmc, out="obsm.p")  # also writes `.uns['p']`
    sc.tl.dpt(pbmc, use="obsm.p")
    np.testing.assert_allclose(
        pbmc.obs["dpt_pseudotime"], expected.obs["dpt_pseudotime"]
    )


@needs.anndata_acc
def test_dpt_out(pbmc) -> None:
    """`out=True` returns the outputs, a mapping names the columns, `.uns` follows the first."""
    sc.tl.diffmap(pbmc)
    obs_before, uns_before = set(pbmc.obs), set(pbmc.uns)
    res = sc.tl.dpt(pbmc, n_branchings=1, out=True)
    assert set(pbmc.obs) == obs_before
    assert set(pbmc.uns) == uns_before
    assert all(v is not None for v in res)

    sc.tl.dpt(
        pbmc, n_branchings=1, out=dict(pseudotime="obs.pt", groups="obs.pt_groups")
    )
    np.testing.assert_allclose(pbmc.obs["pt"], res.pseudotime)
    assert "pt_groups" in pbmc.obs
    assert {"pt_changepoints", "pt_grouptips"} <= set(pbmc.uns)
    assert "dpt_pseudotime" not in pbmc.obs

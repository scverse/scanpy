from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced


@pytest.fixture
def pbmc():
    adata = pbmc68k_reduced()
    adata.uns["iroot"] = 0
    return adata


def test_dpt_use(pbmc) -> None:
    """`use` finds a diffmap written with the same `out`, including its eigenvalues."""
    expected = pbmc.copy()
    sc.tl.diffmap(expected)
    sc.tl.dpt(expected)

    sc.tl.diffmap(pbmc, out=sc.Prefix("p"))
    sc.tl.dpt(pbmc, use=sc.Prefix("p"))
    np.testing.assert_allclose(
        pbmc.obs["dpt_pseudotime"], expected.obs["dpt_pseudotime"]
    )


def test_dpt_out(pbmc) -> None:
    """`out=True` returns the outputs, `sc.Prefix` names the columns and `.uns` entries."""
    sc.tl.diffmap(pbmc)
    obs_before, uns_before = set(pbmc.obs), set(pbmc.uns)
    res = sc.tl.dpt(pbmc, n_branchings=1, out=True)
    assert set(pbmc.obs) == obs_before
    assert set(pbmc.uns) == uns_before
    assert all(v is not None for v in res)

    sc.tl.dpt(pbmc, n_branchings=1, out=sc.Prefix("p"))
    np.testing.assert_allclose(pbmc.obs["p_dpt_pseudotime"], res.pseudotime)
    assert {"p_dpt_changepoints", "p_dpt_grouptips"} <= set(pbmc.uns)
    assert "dpt_pseudotime" not in pbmc.obs

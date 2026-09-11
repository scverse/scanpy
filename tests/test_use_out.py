"""Tests for the `use`/`out` parameters replacing `layer`/`obsm`/`use_rep`/`inplace`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData

import scanpy as sc
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable

pytestmark = needs.anndata_acc


@pytest.fixture
def adata() -> AnnData:
    x = np.random.default_rng(0).poisson(3, (20, 5)).astype("float32")
    adata = AnnData(x.copy())
    adata.layers["counts"] = x.copy()
    adata.obs["batch"] = np.tile([0.0, 1.0], 10)
    return adata


def _orig(adata: AnnData) -> np.ndarray:
    return np.asarray(adata.layers["counts"])


# `use`/`out` on the same-shape transforms.
# Each entry: id -> (call, name of the `.X`-derived result check)
TRANSFORMS: dict[str, Callable[..., object]] = {
    "scale": lambda adata, **kw: sc.pp.scale(adata, zero_center=True, **kw),
    "log1p": sc.pp.log1p,
    "normalize_total": lambda adata, **kw: sc.pp.normalize_total(
        adata, target_sum=1, **kw
    ),
    "regress_out": lambda adata, **kw: sc.pp.regress_out(adata, "batch", **kw),
}


@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_writes_elsewhere(adata: AnnData, name: str) -> None:
    """`out` writes to its own destination and leaves `use`’s source alone."""
    from anndata.acc import A

    before = _orig(adata).copy()
    TRANSFORMS[name](adata, use=A.layers["counts"], out=A.layers["result"])
    assert "result" in adata.layers
    np.testing.assert_array_equal(_orig(adata), before)
    assert not np.allclose(adata.layers["result"], before)


@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_none_returns(adata: AnnData, name: str) -> None:
    """`out=None` returns the result without touching `adata`."""
    from anndata.acc import A

    before = _orig(adata).copy()
    res = TRANSFORMS[name](adata, use=A.layers["counts"], out=None)
    if name == "normalize_total":
        res = res["X"]
    assert res is not None
    assert np.shape(res) == adata.shape
    np.testing.assert_array_equal(_orig(adata), before)
    assert "result" not in adata.layers


@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_default_writes_back(adata: AnnData, name: str) -> None:
    """Without `out`, the result goes back where `use` read it from."""
    from anndata.acc import A

    before = _orig(adata).copy()
    x_before = np.asarray(adata.X).copy()
    assert TRANSFORMS[name](adata, use=A.layers["counts"]) is None
    assert not np.allclose(_orig(adata), before)
    np.testing.assert_array_equal(np.asarray(adata.X), x_before)


@pytest.mark.parametrize("name", TRANSFORMS)
def test_use_str_is_resolved(adata: AnnData, name: str) -> None:
    """A `use` string is always an `anndata.acc` spec, regardless of preset."""
    from anndata.acc import A

    expected = adata.copy()
    TRANSFORMS[name](expected, use=A.layers["counts"], out=A.layers["result"])
    TRANSFORMS[name](adata, use="layers.counts", out="layers.result")
    np.testing.assert_allclose(adata.layers["result"], expected.layers["result"])


def test_use_bare_str_rejected(adata: AnnData) -> None:
    """A bare (non-period-separated) `use` string is not a legacy key."""
    with pytest.raises(ValueError, match=r"Cannot parse accessor"):
        sc.pp.log1p(adata, use="counts")


@pytest.mark.parametrize(
    ("call", "legacy"),
    [
        pytest.param(sc.pp.log1p, "layer", id="log1p"),
        pytest.param(
            lambda a, **kw: sc.pp.scale(a, zero_center=True, **kw), "layer", id="scale"
        ),
        pytest.param(lambda a, **kw: sc.pp.pca(a, n_comps=2, **kw), "layer", id="pca"),
    ],
)
def test_legacy_layer_deprecated(
    adata: AnnData, call: Callable[..., object], legacy: str
) -> None:
    """The legacy `layer` argument still works, but warns."""
    with pytest.warns(FutureWarning, match=rf"argument {legacy} is deprecated"):
        call(adata, **{legacy: "counts"})


def test_use_and_layer_conflict(adata: AnnData) -> None:
    from anndata.acc import A

    with (
        pytest.warns(FutureWarning, match=r"argument layer is deprecated"),
        pytest.raises(TypeError, match=r"`use` cannot be combined with"),
    ):
        sc.pp.log1p(adata, use=A.layers["counts"], layer="counts")


def test_score_genes_out(adata: AnnData) -> None:
    """`score_genes`’ `out` references an `.obs` column; `score_name` is its deprecated version."""
    from anndata.acc import A

    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    genes = list(adata.var_names[:2])

    sc.tl.score_genes(adata, genes, out=A.obs["mine"])
    assert "mine" in adata.obs

    scores = sc.tl.score_genes(adata, genes, out=None)
    assert scores.shape == (adata.n_obs,)
    assert "score" not in adata.obs

    with pytest.warns(FutureWarning, match=r"argument score_name is deprecated"):
        sc.tl.score_genes(adata, genes, score_name="legacy")
    assert "legacy" in adata.obs

    # `score_name` names a column, `out` takes an accessor: not interchangeable
    with pytest.raises(ValueError, match=r"Cannot parse accessor"):
        sc.tl.score_genes(adata, genes, out="legacy")


def test_score_genes_out_and_score_name(adata: AnnData) -> None:
    from anndata.acc import A

    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    with (
        pytest.warns(FutureWarning, match=r"argument score_name is deprecated"),
        pytest.raises(TypeError, match=r"Pass either `out` or `score_name`"),
    ):
        sc.tl.score_genes(
            adata, list(adata.var_names[:2]), out=A.obs["a"], score_name="b"
        )


def test_normalize_total_inplace_deprecated(adata: AnnData) -> None:
    """`inplace=False` is the deprecated way to say `out=None`."""
    with pytest.warns(FutureWarning, match=r"argument inplace is deprecated"):
        res = sc.pp.normalize_total(adata, target_sum=1, inplace=False)
    assert set(res) == {"X", "norm_factor"}
    np.testing.assert_allclose(np.asarray(res["X"]).sum(axis=1), 1, rtol=1e-6)


def test_pca_use(adata: AnnData) -> None:
    """`pca` takes `use` but keeps `key_added` (it writes obsm, varm and uns)."""
    from anndata.acc import A

    sc.pp.pca(adata, n_comps=2, use=A.layers["counts"], key_added="mine")
    assert "mine" in adata.obsm
    assert "mine" in adata.varm
    assert "mine" in adata.uns

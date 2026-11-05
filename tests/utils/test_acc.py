from __future__ import annotations

from contextlib import suppress
from typing import TYPE_CHECKING

import pytest

from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable

with suppress(ImportError):
    from hv_anndata import A

    from scanpy._utils._acc import _resolve, _resolve_all, _resolve_some


pytestmark = [needs.scanpy2]


def test_resolve_vec() -> None:
    assert _resolve("obs.a") == A.obs["a"] == _resolve(A.obs["a"])


def test_resolve_container() -> None:
    from anndata.acc import LayerAcc
    from hv_anndata import AdDim

    counts = LayerAcc("counts", ref_class=AdDim)
    assert _resolve("layers.counts", vec=False) == counts
    assert _resolve(A.layers["counts"], vec=False) == counts


def test_resolve_all() -> None:
    assert _resolve_all(["obs.a", A.var["b"]]) == [A.obs["a"], A.var["b"]]
    assert _resolve_all(["X", A.obsm["pca"]], vec=False) == [A.X, A.obsm["pca"]]


def test_resolve_some() -> None:
    assert _resolve_some("obs.a") == A.obs["a"]
    assert _resolve_some(["obs.a"]) == [A.obs["a"]]


@pytest.mark.parametrize(
    ("call", "match"),
    [
        pytest.param(
            lambda: _resolve(A.X),  # type: ignore
            "Expected a string or an AdDim",
            id="container-as-dim",
        ),
        pytest.param(
            lambda: _resolve(A.obs["a"], vec=False),  # type: ignore
            "Expected a string or a container accessor",
            id="dim-as-container",
        ),
        pytest.param(
            lambda: _resolve_all("obs.a"),
            "Expected a collection of specs",
            id="str-as-collection",
        ),
        pytest.param(
            lambda: _resolve_all(A.X),  # type: ignore
            "Expected a collection of specs",
            id="acc-as-collection",
        ),
    ],
)
def test_resolve_errors(call: Callable[..., object], match: str) -> None:
    with pytest.raises(TypeError, match=match):
        call()

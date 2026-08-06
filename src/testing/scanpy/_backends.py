"""Backend conformance helpers for Scanpy-compatible accelerators."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from anndata import AnnData
from scverse_backends.testing import run_conformance

import scanpy as sc
from scanpy._backends import get_backend
from scanpy._compat import SpBase

if TYPE_CHECKING:
    from collections.abc import Sequence


def _as_numpy(x):
    if hasattr(x, "get"):
        x = x.get()
    if isinstance(x, SpBase):
        return x.toarray()
    return np.asarray(x)


def _conformance_adata() -> AnnData:
    x = np.array(
        [
            [10, 0, 0, 8],
            [9, 1, 0, 7],
            [0, 8, 9, 0],
            [1, 7, 8, 0],
            [5, 5, 1, 1],
            [4, 4, 2, 2],
        ],
        dtype=np.float32,
    )
    return AnnData(x)


def _test_normalize_total(backend_name: str) -> None:
    expected = _conformance_adata()
    actual = expected.copy()

    sc.pp.normalize_total(expected, target_sum=10, backend="cpu")
    sc.pp.normalize_total(actual, target_sum=10, backend=backend_name)

    np.testing.assert_allclose(_as_numpy(actual.X), _as_numpy(expected.X))


def _test_neighbors(backend_name: str) -> None:
    expected = _conformance_adata()
    actual = expected.copy()

    sc.pp.neighbors(expected, n_neighbors=3, use_rep="X", backend="cpu")
    sc.pp.neighbors(actual, n_neighbors=3, use_rep="X", backend=backend_name)

    for key in ["distances", "connectivities"]:
        np.testing.assert_allclose(
            _as_numpy(actual.obsp[key]),
            _as_numpy(expected.obsp[key]),
            rtol=1e-5,
            atol=1e-6,
        )


def _test_umap(backend_name: str) -> None:
    adata = _conformance_adata()
    sc.pp.neighbors(adata, n_neighbors=3, use_rep="X", backend="cpu")

    sc.tl.umap(adata, rng=0, backend=backend_name)

    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert np.isfinite(_as_numpy(adata.obsm["X_umap"])).all()


_TESTS = {
    "normalize_total": _test_normalize_total,
    "neighbors": _test_neighbors,
    "umap": _test_umap,
}


def validate_backend(
    backend_name: str,
    *,
    functions: Sequence[str] | None = None,
    raise_on_failure: bool = True,
) -> dict[str, str]:
    """Run Scanpy's backend conformance checks against an installed backend."""
    return run_conformance(
        backend_name=backend_name,
        tests=_TESTS,
        get_backend=get_backend,
        functions=functions,
        raise_on_failure=raise_on_failure,
    )


__all__ = ["validate_backend"]

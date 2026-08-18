from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse, stats

import scanpy as sc
import scanpy.external as sce

if TYPE_CHECKING:
    from collections.abc import Callable

HASHES = [f"Hash{i}" for i in range(10)]


@pytest.fixture
def counts() -> pd.DataFrame:
    """Hashing counts: cell `i` is a singlet for barcode `i % 10`.

    Except for the first 10 cells (doublets) and the last 10 (negatives).
    """
    rng = np.random.default_rng(0)
    signal = stats.poisson.rvs(1000, 1, 990, random_state=rng)
    doublet_signal = stats.poisson.rvs(1000, 1, 10, random_state=rng)
    x = np.reshape(stats.poisson.rvs(500, 1, 10000, random_state=rng), (1000, 10))
    for idx, signal_count in enumerate(signal):
        x[idx, idx % 10] = signal_count
    for idx, signal_count in enumerate(doublet_signal):
        x[idx, (idx % 10) - 1] = signal_count
    return pd.DataFrame(x, columns=HASHES, index=[f"cell{i}" for i in range(len(x))])


def _in_obs(counts: pd.DataFrame) -> tuple[AnnData, list]:
    from anndata.acc import A

    rng = np.random.default_rng(0)
    adata = AnnData(rng.integers(0, 100, size=counts.shape), obs=counts.copy())
    return adata, A.obs[HASHES]


def _in_x(counts: pd.DataFrame) -> tuple[AnnData, list]:
    from anndata.acc import A

    adata = AnnData(counts.to_numpy(), var=pd.DataFrame(index=HASHES))
    return adata, A.X[:, HASHES]


def _in_x_sparse(counts: pd.DataFrame) -> tuple[AnnData, list]:
    adata, refs = _in_x(counts)
    adata.X = sparse.csr_matrix(adata.X)  # noqa: TID251
    return adata, refs


@pytest.mark.parametrize(
    "make", [_in_obs, _in_x, _in_x_sparse], ids=["obs", "X", "X-sparse"]
)
def test_cell_demultiplexing(counts: pd.DataFrame, make: Callable) -> None:
    adata, refs = make(counts)
    sc.pp.hashsolo(adata, refs)

    expected = pd.array(
        ["Doublet"] * 10
        + np.repeat(HASHES, 98).reshape(98, 10, order="F").ravel().tolist()
        + ["Negative"] * 10,
        dtype="string",
    )
    classification = adata.obs["hashsolo"].array.astype("string")
    # This is a bit flaky, so allow some mismatches:
    # (Series() because of https://github.com/pandas-dev/pandas/issues/63458)
    if pd.Series(expected != classification).sum() > 3:
        # Compare lists for better error message
        assert classification.tolist() == expected.tolist()

    probs = adata.obsm["hashsolo"]
    assert list(probs.columns) == ["negative", "singlet", "doublet"]
    np.testing.assert_allclose(probs.to_numpy().sum(axis=1), 1)


def test_copy(counts: pd.DataFrame) -> None:
    adata, refs = _in_obs(counts)
    copied = sc.pp.hashsolo(adata, refs, copy=True)
    assert "hashsolo" not in adata.obs
    assert "hashsolo" in copied.obs


def test_legacy_api(counts: pd.DataFrame) -> None:
    adata, _ = _in_obs(counts)
    with pytest.warns(FutureWarning, match=r"scanpy\.pp\.hashsolo"):
        sce.pp.hashsolo(adata, HASHES)

    expected_new = adata.copy()
    sc.pp.hashsolo(expected_new, HASHES)
    pd.testing.assert_series_equal(
        adata.obs["Classification"].astype("string"),
        expected_new.obs["hashsolo"].astype("string"),
        check_names=False,
    )
    np.testing.assert_array_equal(
        adata.obs["most_likely_hypothesis"],
        np.argmax(expected_new.obsm["hashsolo"].to_numpy(), axis=1),
    )

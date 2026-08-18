from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse, stats

import scanpy as sc
import scanpy.external as sce
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata.acc import AdRef

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


@pytest.fixture
def adata(counts: pd.DataFrame) -> AnnData:
    """Hashing counts in `.obs`, with unrelated expression data in `X`."""
    rng = np.random.default_rng(0)
    return AnnData(rng.integers(0, 100, size=counts.shape), obs=counts.copy())


@pytest.fixture(
    params=[
        "obs-str",
        pytest.param("obs", marks=needs.anndata_acc),
        pytest.param("X", marks=needs.anndata_acc),
        pytest.param("X-sparse", marks=needs.anndata_acc),
    ]
)
def hashed(
    request: pytest.FixtureRequest, adata: AnnData, counts: pd.DataFrame
) -> tuple[AnnData, Sequence[AdRef | str]]:
    """Build an `AnnData` holding the hashing counts, plus references to them."""
    match request.param:
        case "obs-str":  # plain `.obs` column names work without `anndata.acc`
            return adata, HASHES
        case "obs":
            from anndata.acc import A

            return adata, A.obs[HASHES]
        case "X" | "X-sparse":
            from anndata.acc import A

            x = counts.to_numpy()
            if request.param == "X-sparse":
                x = sparse.csr_matrix(x)  # noqa: TID251
            return AnnData(x, var=pd.DataFrame(index=HASHES)), A.X[:, HASHES]
        case _:
            pytest.fail(f"Unknown param {request.param!r}")


def test_cell_demultiplexing(hashed: tuple[AnnData, Sequence[AdRef | str]]) -> None:
    adata, refs = hashed
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
    assert isinstance(probs, pd.DataFrame)
    assert list(probs.columns) == ["negative", "singlet", "doublet"]
    np.testing.assert_allclose(probs.to_numpy().sum(axis=1), 1)


def test_copy(adata: AnnData) -> None:
    copied = sc.pp.hashsolo(adata, HASHES, copy=True)
    assert "hashsolo" not in adata.obs
    assert "hashsolo" in copied.obs


def test_legacy_api(adata: AnnData) -> None:
    """The deprecated `sce.pp.hashsolo` matches `sc.pp.hashsolo` and needs no `anndata.acc`."""
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

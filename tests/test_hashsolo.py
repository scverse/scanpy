from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

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

    from anndata.acc import AdRef, MultiAcc

HASHES = [f"Hash{i}" for i in range(10)]


class Hashed(NamedTuple):
    """An `AnnData` with hashing counts, how to reference them, and their names."""

    adata: AnnData
    refs: MultiAcc | Sequence[AdRef | str]
    names: Sequence[str]


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


@pytest.fixture(params=["obs", "X", "X-sparse", "obsm-df", "obsm-array"])
def hashed(
    request: pytest.FixtureRequest, adata: AnnData, counts: pd.DataFrame
) -> Hashed:
    """Build an `AnnData` holding the hashing counts, plus references to them."""
    match request.param:
        case "obs":
            from anndata.acc import A

            return Hashed(adata, A.obs[HASHES], HASHES)
        case "X" | "X-sparse":
            from anndata.acc import A

            x = counts.to_numpy()
            if request.param == "X-sparse":
                x = sparse.csr_matrix(x)  # noqa: TID251
            adata = AnnData(x, var=pd.DataFrame(index=HASHES))
            return Hashed(adata, A.X[:, HASHES], HASHES)
        case "obsm-df" | "obsm-array":  # a `MultiAcc` means “all of its columns”
            from anndata.acc import A

            df = request.param == "obsm-df"
            adata.obsm["hto"] = counts.copy() if df else counts.to_numpy()
            # a plain array has no column names, so they fall back to positions
            names = HASHES if df else [str(i) for i in range(len(HASHES))]
            return Hashed(adata, A.obsm["hto"], names)
        case _:
            pytest.fail(f"Unknown param {request.param!r}")


@pytest.mark.parametrize(
    ("n_hashes", "n_noise"),
    [
        pytest.param(1, None, id="1-hash"),
        pytest.param(2, None, id="2-hashes-default"),
        pytest.param(len(HASHES), 0, id="explicit-0"),
    ],
)
@needs.anndata_acc
def test_too_few_noise_barcodes(
    adata: AnnData, n_hashes: int, n_noise: int | None
) -> None:
    """Without a noise barcode, the noise distribution is undefined."""
    from anndata.acc import A

    with pytest.raises(ValueError, match=r"noise barcodes?"):
        sc.pp.hashsolo(adata, A.obs[HASHES[:n_hashes]], n_noise_barcodes=n_noise)


@needs.anndata_acc
def test_cell_demultiplexing(hashed: Hashed) -> None:
    adata, refs, names = hashed
    sc.pp.hashsolo(adata, refs)

    expected = pd.array(
        ["Doublet"] * 10
        + np.repeat(names, 98).reshape(98, 10, order="F").ravel().tolist()
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


@needs.anndata_acc
def test_copy(adata: AnnData) -> None:
    from anndata.acc import A

    copied = sc.pp.hashsolo(adata, A.obs[HASHES], copy=True)
    assert "hashsolo" not in adata.obs
    assert "hashsolo" in copied.obs


def test_legacy_api_needs_no_acc(subtests: pytest.Subtests, adata: AnnData) -> None:
    with pytest.warns(FutureWarning, match=r"scanpy\.pp\.hashsolo"):
        sce.pp.hashsolo(adata, HASHES)


@needs.anndata_acc
def test_legacy_api_matches(subtests: pytest.Subtests, adata: AnnData) -> None:
    from anndata.acc import A

    with pytest.warns(FutureWarning, match=r"scanpy\.pp\.hashsolo"):
        sce.pp.hashsolo(adata, HASHES)

    expected_new = adata.copy()
    sc.pp.hashsolo(expected_new, A.obs[HASHES])
    pd.testing.assert_series_equal(
        adata.obs["Classification"].astype("string"),
        expected_new.obs["hashsolo"].astype("string"),
        check_names=False,
    )
    np.testing.assert_array_equal(
        adata.obs["most_likely_hypothesis"],
        np.argmax(expected_new.obsm["hashsolo"].to_numpy(), axis=1),
    )

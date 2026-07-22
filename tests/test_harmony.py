from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pooch
import pytest
from anndata import AnnData, read_h5ad
from scipy.stats import pearsonr

from scanpy.preprocessing import harmony_integrate
from scanpy.preprocessing._harmony.core import (
    _SUPPRESS_PENALTY,
    _compute_lambda_kb,
    _correction_multi,
    _get_batch_codes,
    _get_theta_array,
)
from testing.scanpy._helpers.data import pbmc68k_reduced

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike


_HARMONY_DATA_BASE = (
    "https://scverse-exampledata.s3.amazonaws.com/rapids-singlecell/harmony_data"
)
_IRCOLITIS_HARMONYPY2_H5AD = (
    "ircolitis_blood_cd8_2048_harmonypy2_2_0_0.h5ad",
    "sha256:c52c4a916fc6b811134dbfb1dc105d83f53195ea3092f277f30c8dbb987d641a",
)
_HARMONYPY2_MULTIKEY_H5AD = (
    "harmonypy2_two_covariates_2_0_0.h5ad",
    "sha256:1ac1542ee31b0660175ed307c29077c5621225af062a0e14d251322abcc1ac46",
)

DATA = dict(
    pca=("pbmc_3500_pcs.tsv.gz", "md5:27e319b3ddcc0c00d98e70aa8e677b10"),
    pca_harmonized=(
        "pbmc_3500_pcs_harmonized.tsv.gz",
        "md5:a7c4ce4b98c390997c66d63d48e09221",
    ),
    meta=("pbmc_3500_meta.tsv.gz", "md5:8c7ca20e926513da7cf0def1211baecb"),
)


def _get_measure(
    x: np.ndarray, base: np.ndarray, norm: Literal["r", "L2"]
) -> np.ndarray:
    """Compute correlation or L2 distance between arrays."""
    if norm == "r":  # Compute per-column correlation
        if x.ndim == 1:
            corr, _ = pearsonr(x, base)
            return corr
        return np.array([pearsonr(x[:, i], base[:, i])[0] for i in range(x.shape[1])])
    if norm == "L2":
        # L2 distance normalized by base norm
        if x.ndim == 1:
            return np.linalg.norm(x - base) / np.linalg.norm(base)
        return np.array([
            np.linalg.norm(x[:, i] - base[:, i]) / np.linalg.norm(base[:, i])
            for i in range(x.shape[1])
        ])
    pytest.fail(f"Unknown {norm=!r}")


@pytest.fixture(scope="module")
def adata_reference_module() -> AnnData:
    """Load reference data once per module (avoids re-reading CSV)."""
    paths = {
        f: pooch.retrieve(
            f"{_HARMONY_DATA_BASE}/{name}",
            known_hash=hash_,
        )
        for f, (name, hash_) in DATA.items()
    }
    dfs = {f: pd.read_csv(path, delimiter="\t") for f, path in paths.items()}
    # Create unique index using row number + cell name
    dfs["meta"].index = [f"{i}_{cell}" for i, cell in enumerate(dfs["meta"]["cell"])]
    return AnnData(
        X=None,
        obs=dfs["meta"],
        obsm={
            "X_pca": dfs["pca"].to_numpy(),
            "harmony_org": dfs["pca_harmonized"].to_numpy(),
        },
    )


@pytest.fixture
def adata_reference(adata_reference_module: AnnData) -> AnnData:
    """Return a fresh copy per test so tests don't mutate shared state."""
    return adata_reference_module.copy()


@pytest.fixture(scope="module")
def adata_harmonypy2_multikey() -> AnnData:
    filename, known_hash = _HARMONYPY2_MULTIKEY_H5AD
    reference_file = pooch.retrieve(
        f"{_HARMONY_DATA_BASE}/{filename}", known_hash=known_hash
    )
    return read_h5ad(reference_file)


@pytest.fixture(scope="module")
def adata_ircolitis_harmonypy2() -> AnnData:
    """Stratified 2,048-cell IRcolitis harmonypy 2.0.0 reference."""
    filename, known_hash = _IRCOLITIS_HARMONYPY2_H5AD
    reference_file = pooch.retrieve(
        f"{_HARMONY_DATA_BASE}/{filename}", known_hash=known_hash
    )
    return read_h5ad(reference_file)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("flavor", ["harmony1", "harmony2"])
def test_harmony_integrate(
    dtype: DTypeLike,
    flavor: Literal["harmony1", "harmony2"],
) -> None:
    """Test that Harmony integrate works."""
    adata = pbmc68k_reduced()
    harmony_integrate(
        adata,
        "bulk_labels",
        correction_method="fast",
        dtype=dtype,
        flavor=flavor,
    )
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_original_correction_rejected() -> None:
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="correction_method must be 'fast'"):
        harmony_integrate(adata, "bulk_labels", correction_method="original")


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate_reference(
    *,
    subtests: pytest.Subtests,
    adata_reference: AnnData,
    dtype: DTypeLike,
) -> None:
    """Test that Harmony1 produces results similar to the reference implementation."""
    harmony_integrate(
        adata_reference,
        "donor",
        dtype=dtype,
        max_iter_harmony=20,
        flavor="harmony1",
    )
    x, base = adata_reference.obsm["harmony_org"], adata_reference.obsm["X_pca_harmony"]

    with subtests.test("r"):
        assert _get_measure(x, base, "r").min() > 0.95
    with subtests.test("L2"):
        assert _get_measure(x, base, "L2").max() < 0.1


def test_harmony_multiple_keys() -> None:
    """Test Harmony models multiple batch keys separately."""
    adata = pbmc68k_reduced()
    # Create a second batch key
    rng = np.random.default_rng(0)
    adata.obs["batch2"] = rng.choice(["A", "B", "C"], size=adata.n_obs)

    harmony_integrate(
        adata,
        ["bulk_labels", "batch2"],
        theta=[2.0, 0.1],
        max_iter_harmony=1,
        rng=0,
    )

    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_multikey_marginal_codes_and_theta() -> None:
    obs = pd.DataFrame({
        "batch": pd.Categorical(
            ["b0", "b1", "b2", "b3", "b4", "b0"],
            categories=["b0", "b1", "b2", "b3", "b4"],
        ),
        "sex": pd.Categorical(["f", "m", "f", "m", "f", "m"], categories=["f", "m"]),
    })

    codes, n_levels = _get_batch_codes(obs, ["batch", "sex"])

    np.testing.assert_array_equal(n_levels, np.array([5, 2], dtype=np.int32))
    np.testing.assert_array_equal(codes[:, 0], np.array([0, 1, 2, 3, 4, 0]))
    np.testing.assert_array_equal(codes[:, 1], np.array([5, 6, 5, 6, 5, 6]))
    np.testing.assert_array_equal(
        _get_theta_array([2.0, 0.1], n_levels, np.dtype(np.float32)),
        np.array([[2, 2, 2, 2, 2, 0.1, 0.1]], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        _get_theta_array(2.0, n_levels, np.dtype(np.float32)),
        np.full((1, 7), 2.0, dtype=np.float32),
    )


@pytest.mark.parametrize("size", [1, 6])
def test_harmony_multikey_theta_rejects_invalid_length(size: int) -> None:
    with pytest.raises(
        ValueError,
        match=r"batch variables \(2\) or categorical levels \(7\)",
    ):
        _get_theta_array([2.0] * size, np.array([5, 2]), np.dtype(np.float32))


def test_harmony_theta_rejects_unsupported_type() -> None:
    with pytest.raises(ValueError, match="theta must be a scalar or an array-like"):
        _get_theta_array({"batch": 2.0}, np.array([5, 2]), np.dtype(np.float32))


def test_harmony_batch_keys_are_nonempty_and_complete() -> None:
    obs = pd.DataFrame({"batch": ["a", None]})
    with pytest.raises(ValueError, match="contains missing values"):
        _get_batch_codes(obs, "batch")
    with pytest.raises(ValueError, match="at least one column"):
        _get_batch_codes(obs, [])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_multikey_correction_matches_dense_design(
    dtype: type[np.floating],
) -> None:
    rng = np.random.default_rng(734)
    levels = np.array([2, 3], dtype=np.int32)
    offsets = np.array([0, 2], dtype=np.int32)
    n_batches = int(levels.sum())
    n_cells, n_pcs, n_clusters = 31, 5, 3

    local_codes = np.column_stack([
        rng.integers(0, level, size=n_cells) for level in levels
    ]).astype(np.int32)
    codes = local_codes + offsets
    x = rng.normal(size=(n_cells, n_pcs)).astype(dtype)
    r = rng.random(size=(n_cells, n_clusters)).astype(dtype)
    r /= r.sum(axis=1, keepdims=True)
    lambda_kb = rng.uniform(0.2, 1.0, size=(n_batches, n_clusters)).astype(dtype)
    lambda_kb[1, 1] = dtype(_SUPPRESS_PENALTY)

    result = _correction_multi(x, codes, n_batches, r, lambda_kb=lambda_kb)

    design = np.zeros((n_cells, n_batches + 1), dtype=dtype)
    design[:, 0] = 1
    for covariate in range(codes.shape[1]):
        design[np.arange(n_cells), codes[:, covariate] + 1] = 1

    expected = x.copy()
    for cluster, r_k in enumerate(r.T):
        active = lambda_kb[:, cluster] < dtype(_SUPPRESS_PENALTY)
        retained = np.concatenate(([True], active))
        weighted_design = r_k[:, None] * design
        gram = design.T @ weighted_design
        active_indices = np.flatnonzero(active) + 1
        gram[active_indices, active_indices] += lambda_kb[active, cluster]
        rhs = weighted_design.T @ x
        w = np.zeros((n_batches + 1, n_pcs), dtype=dtype)
        retained_indices = np.flatnonzero(retained)
        w[retained] = np.linalg.lstsq(
            gram[np.ix_(retained_indices, retained_indices)],
            rhs[retained],
            rcond=None,
        )[0]
        w[0] = 0
        expected -= r_k[:, None] * (design @ w)

    atol = 2e-5 if dtype == np.float32 else 1e-11
    np.testing.assert_allclose(result, expected, atol=atol, rtol=atol)


def test_harmony_multikey_singular_correction_is_finite() -> None:
    rng = np.random.default_rng(734)
    batch = np.resize(["a", "b", "c"], 60)
    adata = AnnData(
        X=None,
        obs=pd.DataFrame(
            {"batch": batch, "duplicate_batch": batch},
            index=[f"cell_{index}" for index in range(60)],
        ),
        obsm={"X_pca": rng.normal(size=(60, 6)).astype(np.float32)},
    )

    harmony_integrate(
        adata,
        ["batch", "duplicate_batch"],
        flavor="harmony1",
        ridge_lambda=0.0,
        n_clusters=3,
        max_iter_harmony=1,
        max_iter_clustering=2,
        block_proportion=1.0,
        rng=734,
    )

    assert np.isfinite(adata.obsm["X_pca_harmony"]).all()


@pytest.mark.parametrize(
    ("case", "theta", "n_clusters", "max_iter_harmony"),
    [
        ("nclust1", [2.0, 0.1], 1, 1),
        ("nclust4", [2.0, 0.1], 4, 10),
        ("nclust4_batch_only", [2.0, 0.0], 4, 10),
        ("nclust4_sex_only", [0.0, 2.0], 4, 10),
    ],
)
def test_harmony2_multikey_reference(
    adata_harmonypy2_multikey: AnnData,
    case: str,
    theta: list[float],
    n_clusters: int,
    max_iter_harmony: int,
) -> None:
    adata = adata_harmonypy2_multikey.copy()
    reference = adata.obsm[f"harmony2_ref_{case}"].copy()

    harmony_integrate(
        adata,
        ["batch", "sex"],
        theta=theta,
        flavor="harmony2",
        dtype=np.float64,
        sigma=0.1,
        n_clusters=n_clusters,
        max_iter_harmony=max_iter_harmony,
        max_iter_clustering=4,
        tol_clustering=1e-3,
        tol_harmony=1e-2,
        block_proportion=0.05,
        rng=734,
        alpha=0.2,
        batch_prune_threshold=1e-5,
    )

    result = adata.obsm["X_pca_harmony"]
    assert _get_measure(reference, result, "r").min() > 0.95
    assert _get_measure(reference, result, "L2").max() < 0.1


def test_harmony_custom_parameters() -> None:
    """Test Harmony with custom parameters."""
    adata = pbmc68k_reduced()
    harmony_integrate(
        adata,
        "bulk_labels",
        theta=1.5,
        sigma=0.15,
        n_clusters=50,
        max_iter_harmony=5,
        ridge_lambda=0.5,
        flavor="harmony1",
    )
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_no_nan_output() -> None:
    """Test that Harmony output contains no NaN values."""
    adata = pbmc68k_reduced()
    harmony_integrate(adata, "bulk_labels")
    assert not np.isnan(adata.obsm["X_pca_harmony"]).any()


def test_harmony_input_validation(subtests) -> None:
    """Test that Harmony raises errors for invalid inputs."""
    adata = pbmc68k_reduced()

    with subtests.test("no basis"), pytest.raises(ValueError, match="not available"):
        harmony_integrate(adata, "bulk_labels", basis="nonexistent")
    with subtests.test("no key"), pytest.raises(KeyError):
        harmony_integrate(adata, "nonexistent_key")


def test_harmony_invalid_flavor() -> None:
    """Test that invalid flavor raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="flavor must be"):
        harmony_integrate(adata, "bulk_labels", flavor="harmony3")


@pytest.mark.parametrize("bad_alpha", [-0.1, 0.0, float("inf"), float("nan")])
def test_harmony_integrate_bad_alpha(bad_alpha: float) -> None:
    """Non-positive or non-finite alpha with flavor='harmony2' raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="alpha must be a finite positive"):
        harmony_integrate(adata, "bulk_labels", alpha=bad_alpha)


@pytest.mark.parametrize("bad_threshold", [-0.1, 1.5, 2.0])
def test_harmony_integrate_bad_prune_threshold(bad_threshold: float) -> None:
    """batch_prune_threshold outside [0, 1] raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="batch_prune_threshold must be in"):
        harmony_integrate(adata, "bulk_labels", batch_prune_threshold=bad_threshold)


def test_harmony_flavor_warnings() -> None:
    """Test that flavor-incompatible parameter warnings are raised."""
    adata = pbmc68k_reduced()

    # harmony2 with ridge_lambda should warn
    with pytest.warns(UserWarning, match="ridge_lambda is ignored"):
        harmony_integrate(
            adata,
            "bulk_labels",
            flavor="harmony2",
            ridge_lambda=2.0,
            max_iter_harmony=1,
        )

    # harmony1 with alpha should warn
    with pytest.warns(UserWarning, match="alpha is ignored"):
        harmony_integrate(
            adata, "bulk_labels", flavor="harmony1", alpha=0.5, max_iter_harmony=1
        )

    # harmony1 with batch_prune_threshold should warn
    with pytest.warns(UserWarning, match="batch_prune_threshold is ignored"):
        harmony_integrate(
            adata,
            "bulk_labels",
            flavor="harmony1",
            batch_prune_threshold=0.01,
            max_iter_harmony=1,
        )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_pruning(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb suppresses correction for N_b==0 and below-threshold pairs."""
    n_batches, n_clusters = 4, 3
    alpha = 0.2
    threshold = 1e-5
    sentinel = dtype(_SUPPRESS_PENALTY)

    n_b = np.array([0, 100, 1, 50], dtype=dtype)
    o = np.array(
        [[0, 0, 0], [30, 40, 30], [0, 0, 1], [20, 15, 15]],
        dtype=dtype,
    )
    e = np.ones((n_batches, n_clusters), dtype=dtype) * 10

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=alpha,
        threshold=threshold,
        ridge_lambda=1.0,
        dynamic_lambda=True,
    )

    # batch 0 (N_b==0): all clusters must be sentinel
    assert np.all(result[0] == sentinel)
    # batch 1 (well-represented): should be alpha * E = 2.0
    np.testing.assert_allclose(result[1], np.full(n_clusters, alpha * 10, dtype=dtype))
    # batch 2, clusters 0,1 (O/N_b = 0/1 < threshold): sentinel
    assert result[2, 0] == sentinel
    assert result[2, 1] == sentinel
    # batch 2, cluster 2 (O/N_b = 1/1 = 1.0 >= threshold): alpha * E
    np.testing.assert_allclose(result[2, 2], dtype(alpha * 10))
    # batch 3: all alpha * E
    np.testing.assert_allclose(result[3], np.full(n_clusters, alpha * 10, dtype=dtype))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_dynamic_false(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb returns uniform ridge_lambda when dynamic_lambda=False."""
    n_batches, n_clusters = 3, 5
    e = np.ones((n_batches, n_clusters), dtype=dtype)
    o = np.ones((n_batches, n_clusters), dtype=dtype)
    n_b = np.ones(n_batches, dtype=dtype)

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=0.5,
        threshold=1e-5,
        ridge_lambda=1.0,
        dynamic_lambda=False,
    )
    np.testing.assert_array_equal(result, np.full_like(e, 1.0))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_fixed_ridge_zero(dtype: type[np.floating]) -> None:
    """Fixed zero ridge still guards a zero correction denominator."""
    sentinel = dtype(_SUPPRESS_PENALTY)
    e = np.array([[5.0, 0.0]], dtype=dtype)
    o = np.array([[10.0, 0.0]], dtype=dtype)
    n_b = np.array([100.0], dtype=dtype)

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=0.2,
        threshold=None,
        ridge_lambda=0.0,
        dynamic_lambda=False,
    )

    assert result[0, 0] == dtype(0.0)
    assert result[0, 1] == sentinel


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_zero_denom(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb guards against O==0 and E==0 (zero-denominator)."""
    sentinel = dtype(_SUPPRESS_PENALTY)
    e = np.array([[0.0, 5.0]], dtype=dtype)
    o = np.array([[0.0, 10.0]], dtype=dtype)
    n_b = np.array([100.0], dtype=dtype)

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=0.2,
        threshold=None,
        ridge_lambda=1.0,
        dynamic_lambda=True,
    )
    # (0,0): O+lambda_kb = 0+0 = 0 -> sentinel
    assert result[0, 0] == sentinel
    # (0,1): normal -> alpha * E = 1.0
    np.testing.assert_allclose(result[0, 1], dtype(1.0))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony2_ircolitis_reference(
    adata_ircolitis_harmonypy2: AnnData,
    dtype: type[np.floating],
) -> None:
    """Harmony2 on a real 11-batch subset matches harmonypy 2.0.0."""
    adata = adata_ircolitis_harmonypy2.copy()
    harmony_integrate(
        adata,
        "batch",
        theta=2.0,
        flavor="harmony2",
        dtype=dtype,
        sigma=0.1,
        n_clusters=2,
        max_iter_harmony=10,
        max_iter_clustering=4,
        tol_clustering=1e-3,
        tol_harmony=1e-2,
        block_proportion=0.05,
        rng=734,
        alpha=0.2,
        batch_prune_threshold=1e-5,
    )

    reference = adata.obsm["harmony2_ref"]
    result = adata.obsm["X_pca_harmony"]
    assert _get_measure(reference, result, "r").min() > 0.95
    assert _get_measure(reference, result, "L2").max() < 0.1

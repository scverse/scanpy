from __future__ import annotations

from functools import partial
from itertools import pairwise
from typing import TYPE_CHECKING

import numba
import numpy as np
import pandas as pd
import pooch
import pytest
from anndata import AnnData, read_h5ad
from scipy.stats import pearsonr

from scanpy.preprocessing import harmony_integrate
from scanpy.preprocessing._harmony import core
from scanpy.preprocessing._harmony._kernels import (
    _apply_correction,
    _cluster_pass,
    _joint_stats,
)
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

# These settings were used to generate the stored harmonypy 2.0.0 reference outputs.
# Keep them fixed here so tests only vary the parameters relevant to each fixture.
_integrate_harmony2_reference = partial(
    harmony_integrate,
    flavor="harmony2",
    sigma=0.1,
    max_iter_clustering=4,
    tol_clustering=1e-3,
    tol_harmony=1e-2,
    block_proportion=0.05,
    rng=734,
    alpha=0.2,
    batch_prune_threshold=1e-5,
    seeding="kmeans",  # harmonypy seeds with k-means
    stop_rule="reference",  # and records a round's objective only if clustering converged
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


def _dense_correction(
    x: np.ndarray,
    codes: np.ndarray,
    n_batches: int,
    r: np.ndarray,
    lambda_kb: np.ndarray,
) -> np.ndarray:
    """Ridge correction by an explicit one-hot design and a least-squares solve per cluster."""
    dtype = x.dtype.type
    n_cells, n_pcs = x.shape
    design = np.zeros((n_cells, n_batches + 1), dtype=dtype)
    design[:, 0] = 1
    for covariate in range(codes.shape[1]):
        design[np.arange(n_cells), codes[:, covariate] + 1] = 1

    expected = x.copy()
    for cluster, r_k in enumerate(r.T):
        active = lambda_kb[:, cluster] < dtype(_SUPPRESS_PENALTY)
        retained = np.concatenate(([True], active))
        sqrt_r_k = np.sqrt(r_k)
        weighted_design = sqrt_r_k[:, None] * design[:, retained]
        weighted_x = sqrt_r_k[:, None] * x

        penalty_design = np.zeros(
            (active.sum(), retained.sum()),
            dtype=dtype,
        )
        penalty_design[:, 1:] = np.diag(np.sqrt(lambda_kb[active, cluster]))
        augmented_design = np.vstack((weighted_design, penalty_design))
        augmented_x = np.vstack((
            weighted_x,
            np.zeros((active.sum(), n_pcs), dtype=dtype),
        ))

        w = np.zeros((n_batches + 1, n_pcs), dtype=dtype)
        w[retained] = np.linalg.lstsq(
            augmented_design,
            augmented_x,
            rcond=None,
        )[0]
        w[0] = 0
        expected -= r_k[:, None] * (design @ w)
    return expected


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_multikey_correction_matches_dense_design(
    dtype: type[np.floating],
) -> None:
    """Compare the optimized correction with a direct dense ridge-regression solve.

    Production avoids materializing the one-hot design and exploits its marginal
    structure. This independent oracle explicitly builds that design and solves
    the weighted, regularized least-squares problem, catching errors in category
    indexing, cross-covariate terms, and suppressed penalties.
    """
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
    expected = _dense_correction(x, codes, n_batches, r, lambda_kb)

    atol = 2e-5 if dtype == np.float32 else 1e-11
    np.testing.assert_allclose(result, expected, atol=atol, rtol=atol)


def test_harmony_multikey_correction_singular_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A system the solver rejects as singular falls back to least squares."""
    rng = np.random.default_rng(734)
    x = rng.normal(size=(12, 3))
    codes = np.column_stack((
        np.resize([0, 1], len(x)),
        np.resize([2, 3, 4], len(x)),
    )).astype(np.int32)
    r = rng.random(size=(len(x), 2))
    r /= r.sum(axis=1, keepdims=True)
    lambda_kb = np.ones((5, 2))

    expected = _correction_multi(x, codes, 5, r, lambda_kb=lambda_kb)

    def singular(*_args, **_kwargs):
        raise np.linalg.LinAlgError

    monkeypatch.setattr(np.linalg, "solve", singular)
    result = _correction_multi(x, codes, 5, r, lambda_kb=lambda_kb)

    np.testing.assert_allclose(result, expected, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("ridge_lambda", [0.0, -0.1, float("inf"), float("nan")])
@pytest.mark.parametrize("multikey", [False, True])
def test_harmony1_rejects_bad_ridge_lambda(
    ridge_lambda: float, *, multikey: bool
) -> None:
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
    key = ["batch", "duplicate_batch"] if multikey else "batch"

    with pytest.raises(ValueError, match="ridge_lambda must be a finite positive"):
        harmony_integrate(
            adata,
            key,
            flavor="harmony1",
            ridge_lambda=ridge_lambda,
            n_clusters=3,
            max_iter_harmony=1,
            max_iter_clustering=2,
            block_proportion=1.0,
            rng=734,
        )


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

    _integrate_harmony2_reference(
        adata,
        ["batch", "sex"],
        theta=theta,
        dtype=np.float64,
        n_clusters=n_clusters,
        max_iter_harmony=max_iter_harmony,
    )

    result = adata.obsm["X_pca_harmony"]
    assert _get_measure(reference, result, "r").min() > 0.95
    assert _get_measure(reference, result, "L2").max() < 0.1


def _within_group_r2(x: np.ndarray, group: np.ndarray, covariate: np.ndarray) -> float:
    """Share of the within-group variance of ``x`` explained by ``covariate``."""
    total = explained = 0.0
    for g in np.unique(group):
        xg, cg = x[group == g], covariate[group == g]
        total += ((xg - xg.mean(0)) ** 2).sum()
        explained += sum(
            (cg == c).sum() * ((xg[cg == c].mean(0) - xg.mean(0)) ** 2).sum()
            for c in np.unique(cg)
        )
    return explained / total


def test_harmony2_multikey_default_seeding(adata_harmonypy2_multikey: AnnData) -> None:
    """The default seeding integrates at least as well as harmonypy on its multikey fixture.

    It starts further from the k-means fixed point than harmonypy, so it runs more rounds
    and lands farther from the stored output: compare the biology, not the coordinates.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    adata = adata_harmonypy2_multikey.copy()
    x = adata.obsm["X_pca"]
    group = KMeans(4, n_init=10, random_state=0).fit_predict(
        x / np.linalg.norm(x, axis=1, keepdims=True)
    )
    assert np.all(
        np.bincount(group) == 20
    )  # four groups of 20, each balanced over batch x sex
    reference = adata.obsm["harmony2_ref_nclust4"]

    harmony_integrate(
        adata,
        ["batch", "sex"],
        theta=[2.0, 0.1],
        n_clusters=4,
        max_iter_clustering=4,
        tol_clustering=1e-3,
        tol_harmony=1e-2,
        rng=734,
    )

    result = adata.obsm["X_pca_harmony"]
    batch = adata.obs["batch"].to_numpy()
    assert _within_group_r2(result, group, batch) <= _within_group_r2(
        reference, group, batch
    )
    assert silhouette_score(result, group) >= silhouette_score(reference, group)


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
    _integrate_harmony2_reference(
        adata,
        "batch",
        theta=2.0,
        dtype=dtype,
        n_clusters=2,
        max_iter_harmony=10,
    )

    reference = adata.obsm["harmony2_ref"]
    result = adata.obsm["X_pca_harmony"]
    assert _get_measure(reference, result, "r").min() > 0.95
    assert _get_measure(reference, result, "L2").max() < 0.1


# Kernels and driver of the numba implementation

_N_BATCH, _K, _SIGMA = 3, 5, 0.1
_DTYPES = [np.float32, np.float64]


def _toy(
    n: int = 600, d: int = 8, seed: int = 0, dtype: DTypeLike = np.float32
) -> tuple[np.ndarray, pd.DataFrame]:
    """Four cell types, three batches with a shift each, and a second key without effect."""
    rng = np.random.default_rng(seed)
    cell_type, batch = rng.integers(4, size=n), rng.integers(_N_BATCH, size=n)
    x = (
        3 * rng.normal(size=(4, d))[cell_type]
        + rng.normal(size=(_N_BATCH, d))[batch]
        + 0.5 * rng.normal(size=(n, d))
    )
    obs = pd.DataFrame({
        "batch": pd.Categorical(batch.astype(str)),
        "half": pd.Categorical((np.arange(n) % 2).astype(str)),
    })
    return x.astype(dtype), obs


def _state(dtype: DTypeLike, keys: str | list[str] = "batch", seed: int = 1):
    """Build a random soft assignment with consistent counts."""
    x, obs = _toy(dtype=dtype)
    z = core._normalize_rows_l2(x)
    codes, n_levels = _get_batch_codes(obs, keys)
    r = np.random.default_rng(seed).random((len(x), _K)).astype(dtype) ** 8
    r /= r.sum(axis=1, keepdims=True)
    o = np.zeros((int(n_levels.sum()), _K), dtype=dtype)
    for covariate in codes.T:
        np.add.at(o, covariate, r)
    n_b = np.bincount(codes.ravel()).astype(dtype)
    return x, z, codes, n_levels, r, o, n_b, n_b / dtype(len(x))


def _bounds(n: int, n_blocks: int = 19) -> np.ndarray:
    sizes = [len(a) for a in np.array_split(np.arange(n), n_blocks)]
    return np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)


def _clustering_iteration(  # noqa: PLR0917
    z, codes, pr_b, r, o, theta, perm, *, stabilized, sum_penalties
):
    """One clustering iteration in float64 NumPy, block by block, with either combination rule."""
    r, o = r.astype(np.float64), o.astype(np.float64)
    y = core._normalize_rows_l2(r.T @ z).astype(np.float64)
    r_sum = r.sum(axis=0)
    for idx in np.array_split(perm, 19):
        for covariate in codes[idx].T:
            np.add.at(o, covariate, -r[idx])
        r_sum -= r[idx].sum(axis=0)
        e = pr_b[:, None] * r_sum
        den = o + e + 1 if stabilized else o + 1
        log_factor = theta[:, None] * (np.log(e + 1) - np.log(den))
        per_key = log_factor[codes[idx]]  # cells x keys x K
        if sum_penalties:
            log_penalty = np.logaddexp.reduce(per_key, axis=1)
        else:
            log_penalty = per_key.sum(axis=1)
        logit = -2 / _SIGMA * (1 - z[idx] @ y.T) + log_penalty
        new = np.exp(logit - logit.max(axis=1, keepdims=True))
        r[idx] = new / new.sum(axis=1, keepdims=True)
        for covariate in codes[idx].T:
            np.add.at(o, covariate, r[idx])
        r_sum += r[idx].sum(axis=0)
    return r, o


@pytest.mark.parametrize("dtype", _DTYPES)
@pytest.mark.parametrize(
    ("keys", "stabilized", "n_chunks", "scatter", "sum_penalties"),
    [
        ("batch", True, 32, False, False),
        (["batch", "half"], False, 3, True, False),
        ("batch", True, 3, False, False),
        (["batch", "half"], False, 32, True, True),
        (["batch", "half"], True, 32, False, True),
    ],
)
def test_cluster_pass_matches_numpy(
    dtype, keys, *, stabilized, n_chunks, scatter, sum_penalties
) -> None:
    x, z, codes, _, r0, o0, _, pr_b = _state(dtype, keys)
    theta = np.full(o0.shape[0], 2, dtype=dtype)
    perm = np.random.default_rng(7).permutation(len(x))
    kw = dict(stabilized=stabilized, sum_penalties=sum_penalties)
    r_ref, o_ref = _clustering_iteration(z, codes, pr_b, r0, o0, theta, perm, **kw)

    y = core._normalize_rows_l2(r0.T @ z)
    r, o, r_sum = r0.copy(), o0.copy(), r0.sum(axis=0)
    dist, ent, o_fresh, y_t = _cluster_pass(
        *(z, np.ascontiguousarray(y.T), r, codes, o, r_sum, pr_b, theta, perm),
        *(_bounds(len(x)), dtype(-2 / _SIGMA), stabilized, n_chunks, scatter),
        sum_penalties,
    )
    atol = 2e-5 if dtype == np.float32 else 1e-10
    np.testing.assert_allclose(r, r_ref, atol=atol)
    np.testing.assert_allclose(o, o_ref, atol=50 * atol)
    # objective parts, summed on the fly
    assert dist == pytest.approx(np.sum(r_ref * 2 * (1 - z @ y.T)), rel=1e-4)
    assert ent == pytest.approx(np.sum(r_ref * np.log(r_ref + 1e-300)), rel=1e-4)
    # invariants
    assert r.dtype == o.dtype == o_fresh.dtype == y_t.dtype == dtype
    np.testing.assert_allclose(r.sum(axis=1), 1, atol=1e-5)
    np.testing.assert_allclose(r_sum, r.sum(axis=0), rtol=1e-4)
    exact = np.zeros_like(o)
    for covariate in codes.T:
        np.add.at(exact, covariate, r)
    np.testing.assert_allclose(o_fresh, exact, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(y_t.T, core._normalize_rows_l2(r.T @ z), atol=1e-5)
    if isinstance(keys, list):  # the two rules really differ on this input
        other, _ = _clustering_iteration(
            *(z, codes, pr_b, r0, o0, theta, perm),
            stabilized=stabilized,
            sum_penalties=not sum_penalties,
        )
        assert np.abs(r_ref - other).max() > 1e-3


@pytest.mark.parametrize("dtype", _DTYPES)
def test_cluster_pass_sum_rule_handles_large_penalties(dtype) -> None:
    """Summed penalty factors far above exp(80) keep their order between clusters."""
    x, z, codes, _, r0, o0, _, pr_b = _state(dtype, ["batch", "half"])
    r0[codes[:, 0] == 0, :2] = 0  # batch 0 absent from clusters 0 and 1
    r0[:, 1] *= 0.3  # with different expected counts
    r0 /= r0.sum(axis=1, keepdims=True)
    o0[:] = 0
    for covariate in codes.T:
        np.add.at(o0, covariate, r0)
    theta = np.full(o0.shape[0], 30, dtype=dtype)
    perm = np.random.default_rng(7).permutation(len(x))
    r_ref, _ = _clustering_iteration(
        *(z, codes, pr_b, r0, o0, theta, perm), stabilized=False, sum_penalties=True
    )
    y = core._normalize_rows_l2(r0.T @ z)
    r, o, r_sum = r0.copy(), o0.copy(), r0.sum(axis=0)
    _cluster_pass(
        *(z, np.ascontiguousarray(y.T), r, codes, o, r_sum, pr_b, theta, perm),
        *(_bounds(len(x)), dtype(-2 / _SIGMA), False, 32, False, True),
    )
    np.testing.assert_allclose(r, r_ref, atol=5e-5 if dtype == np.float32 else 1e-10)


def test_cluster_pass_without_penalty_is_soft_kmeans() -> None:
    x, z, codes, _, r0, _, _, pr_b = _state(np.float32)
    y = core._normalize_rows_l2(r0.T @ z)
    r = np.zeros_like(r0)
    o, r_sum = np.zeros((_N_BATCH, _K), np.float32), np.zeros(_K, np.float32)
    n = len(x)
    _cluster_pass(
        *(z, np.ascontiguousarray(y.T), r, codes, o, r_sum, pr_b),
        *(np.zeros(_N_BATCH, np.float32), np.arange(n), np.array([0, n])),
        *(np.float32(-20), True, 3, False, False),
    )
    expected = np.exp(-20 * (1 - z @ y.T).astype(np.float64))
    expected /= expected.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(r, expected, atol=1e-5)
    # sum_k o[b, k] = N_b
    np.testing.assert_allclose(o.sum(axis=1), np.bincount(codes[:, 0]), rtol=1e-5)


def _random_design(rng: np.random.Generator, n: int, sizes: tuple[int, ...]):
    """Random levels for several keys and a random soft assignment, in float64."""
    obs = pd.DataFrame({
        f"key{c}": pd.Categorical(
            rng.integers(size, size=n).astype(str),
            categories=[str(v) for v in range(size)],
        )
        for c, size in enumerate(sizes)
    })
    codes, n_levels = _get_batch_codes(obs, list(obs.columns))
    x = rng.normal(size=(n, 6))
    r = rng.random((n, _K)) ** 4
    r /= r.sum(axis=1, keepdims=True)
    return codes, n_levels, x, r


@pytest.mark.parametrize("dtype", _DTYPES)
# big block first / last / a size-1 key / one key
@pytest.mark.parametrize("sizes", [(4, 9, 3), (4, 3, 9), (1, 8), (7,)])
def test_correction_with_pruning_and_an_empty_cluster(dtype, sizes) -> None:
    rng = np.random.default_rng(sum(sizes))
    codes, n_levels, x, r = _random_design(rng, 500, sizes)
    r[:, 3] = 0  # an empty cluster
    r /= r.sum(axis=1, keepdims=True)
    n_batches = int(n_levels.sum())
    lam = rng.uniform(0.2, 2.0, size=(n_batches, _K))
    # pruned pairs in every block, several per cluster
    lam[rng.random((n_batches, _K)) < 0.25] = _SUPPRESS_PENALTY
    lam[:, 1] = _SUPPRESS_PENALTY  # everything pruned in one cluster
    lam[: n_levels[0], 2] = _SUPPRESS_PENALTY  # a whole key pruned in another

    result = _correction_multi(
        x.astype(dtype), codes, n_batches, r.astype(dtype), lambda_kb=lam.astype(dtype)
    )
    expected = _dense_correction(x, codes, n_batches, r, lam)
    assert result.dtype == dtype
    np.testing.assert_allclose(
        result, expected, atol=5e-5 if dtype == np.float32 else 1e-10
    )


def test_correction_with_nested_keys_and_unused_levels() -> None:
    rng = np.random.default_rng(1)
    n = 600
    sample = rng.integers(12, size=n)
    obs = pd.DataFrame({
        # donors 4 and 5 never occur
        "donor": pd.Categorical(
            (sample // 3).astype(str), categories=[str(v) for v in range(6)]
        ),
        "sample": pd.Categorical(sample.astype(str)),
    })
    codes, _ = _get_batch_codes(obs, ["donor", "sample"])
    x = rng.normal(size=(n, 6))
    r = rng.random((n, _K)) ** 4
    r /= r.sum(axis=1, keepdims=True)
    lam = np.full((18, _K), 0.7)
    lam[5] = (
        _SUPPRESS_PENALTY  # one unused level pruned, the other left with a finite penalty
    )
    result = _correction_multi(x, codes, 18, r, lambda_kb=lam)
    np.testing.assert_allclose(
        result, _dense_correction(x, codes, 18, r, lam), atol=1e-10
    )


@pytest.mark.parametrize("dtype", _DTYPES)
def test_apply_correction_normalises_and_returns_centroids(dtype) -> None:
    x, _, codes, n_levels, r, o, n_b, pr_b = _state(dtype, ["batch", "half"])
    lam = _compute_lambda_kb(
        pr_b[:, None] * r.sum(axis=0),
        o=o,
        n_b=n_b,
        alpha=0.2,
        threshold=1e-5,
        ridge_lambda=1.0,
        dynamic_lambda=True,
    )
    joint, levels = core._joint_categories(codes)
    s, o_joint = _joint_stats(
        x, r, joint, np.argsort(joint, kind="stable"), len(levels), 3
    )
    w = core._solve_correction(s, o_joint, levels, n_levels, lam, _SUPPRESS_PENALTY)
    z_hat, z_norm = np.empty_like(x), np.empty_like(x)
    y_t = _apply_correction(x, r, joint, w, z_hat, z_norm, 3)
    assert z_hat.dtype == z_norm.dtype == y_t.dtype == dtype
    np.testing.assert_allclose(
        z_hat,
        _correction_multi(x, codes, int(n_levels.sum()), r, lambda_kb=lam),
        atol=1e-5,
    )
    np.testing.assert_allclose(np.linalg.norm(z_norm, axis=1), 1, atol=1e-5)
    np.testing.assert_allclose(y_t.T, core._normalize_rows_l2(r.T @ z_norm), atol=1e-5)


@pytest.mark.parametrize("dtype", _DTYPES)
# normal, many categories across chunks, more chunks than cells (many empty)
@pytest.mark.parametrize(("n_chunks", "n_joint"), [(3, 7), (32, 300), (1000, 2)])
def test_joint_stats_equal_direct_sums(dtype, n_chunks: int, n_joint: int) -> None:
    rng = np.random.default_rng(n_chunks + n_joint)
    n = 400
    x = rng.normal(size=(n, 5)).astype(dtype)
    r = rng.random((n, _K)).astype(dtype)
    joint = rng.integers(n_joint, size=n).astype(np.int32)
    joint[:150] = 0  # one large category that spans many chunks
    s, o = _joint_stats(
        x, r, joint, np.argsort(joint, kind="stable"), n_joint, n_chunks
    )
    s_ref, o_ref = np.zeros((n_joint, _K, 5)), np.zeros((n_joint, _K))
    np.add.at(o_ref, joint, r.astype(np.float64))
    np.add.at(s_ref, joint, r[:, :, None].astype(np.float64) * x[:, None, :])
    assert s.dtype == o.dtype == dtype
    tol = 1e-4 if dtype == np.float32 else 1e-10
    np.testing.assert_allclose(o, o_ref, rtol=tol, atol=tol)
    np.testing.assert_allclose(s, s_ref, rtol=tol, atol=tol)


def test_joint_categories() -> None:
    rng = np.random.default_rng(0)
    for n_keys, n_levels in ((1, 5), (3, 4), (16, 2), (40, 3)):
        codes = rng.integers(n_levels, size=(3000, n_keys)) + n_levels * np.arange(
            n_keys
        )
        joint, levels = core._joint_categories(codes.astype(np.int32))
        assert joint.dtype == np.int32
        assert joint.flags.c_contiguous
        np.testing.assert_array_equal(levels[joint], codes)
        assert len(levels) == len(np.unique(codes, axis=0)) == joint.max() + 1


@pytest.mark.parametrize("keys", ["batch", ["batch", "half"]])
def test_harmonize_float32_matches_float64(keys) -> None:
    x, obs = _toy(n=1500, dtype=np.float64)
    kw = dict(seeding="kmeans++", rng=0)  # the same seeds in both precisions
    got64, info = core.harmonize(x, obs, keys, **kw)
    got32, _ = core.harmonize(x.astype(np.float32), obs, keys, **kw)
    assert got64.dtype == np.float64
    assert got32.dtype == np.float32
    assert np.abs(got64 - x).max() > 0.1  # it did correct something
    assert len(info["iterations"]) <= 10
    r = [np.corrcoef(got32[:, j], got64[:, j])[0, 1] for j in range(x.shape[1])]
    assert min(r) > 0.9999
    rel = np.linalg.norm(got32 - got64, axis=0) / np.linalg.norm(got64, axis=0)
    assert rel.max() < 1e-3


def test_harmonize_many_small_batches() -> None:
    """A design the rule sends through the scatter path agrees across precisions."""
    assert not core._use_scatter(166, 1, 584_944 // 19)  # atlas-like: accumulators
    assert core._use_scatter(
        4000, 1, 20_000 // 19
    )  # thousands of tiny batches: scatters
    assert core._use_scatter(11, 1, 10 * core.N_CHUNKS)
    assert not core._use_scatter(10, 1, 10 * core.N_CHUNKS)
    x, obs = _toy(n=1500, dtype=np.float64)
    obs = obs.assign(well=pd.Categorical(np.arange(1500) % 250))
    assert core._use_scatter(250, 1, 1500 // 19)
    got64, _ = core.harmonize(x, obs, "well", seeding="kmeans++", rng=0)
    got32, _ = core.harmonize(
        x.astype(np.float32), obs, "well", seeding="kmeans++", rng=0
    )
    rel = np.linalg.norm(got32 - got64, axis=0) / np.linalg.norm(got64, axis=0)
    assert rel.max() < 1e-3


def test_harmonize_stop_rules() -> None:
    assert not core._converged([10.0], 1e-2, window=1, signed=False)
    assert core._converged([10.0, 10.05], 1e-2, window=1, signed=False)
    assert not core._converged([10.0, 11.0], 1e-2, window=1, signed=False)
    # a rise passes the signed test
    assert core._converged([10.0, 11.0], 1e-2, window=1, signed=True)
    # needs window + 1 values
    assert not core._converged([10.0, 9.0, 8.0], 1e-2, window=3, signed=True)
    x, obs = _toy()
    for rule in ("reference_record", "paper", "paper_window3"):
        _, info = core.harmonize(x, obs, "batch", stop_rule=rule, rng=0)
        objs = info["objectives"]
        assert len(objs) == len(info["iterations"]) + 1  # every round recorded
        assert max(info["iterations"]) <= 200
        tol = core.STOP_RULES[rule][2]["harmony2"][2]
        passed = [old - new < tol * abs(old) for old, new in pairwise(objs)]
        assert not any(passed[:-1])  # no earlier round met the stop test
        # the last one did, or the cap ended the run
        assert passed[-1] or len(info["iterations"]) == 10


def test_harmonize_reference_rule_records_only_converged_rounds() -> None:
    x, obs = _toy()
    # a window of 3 can never pass with one iteration per round
    kw = dict(max_iter_clustering=1, rng=0)
    _, info = core.harmonize(x, obs, "batch", stop_rule="reference", **kw)
    assert info["iterations"] == [1] * 10
    assert len(info["objectives"]) == 1  # only the initial objective
    _, info = core.harmonize(x, obs, "batch", **kw)
    assert len(info["objectives"]) == len(info["iterations"]) + 1


@pytest.mark.parametrize("seeding", ["kmeans++", "afkmc2", "paper_nearest"])
def test_harmonize_seedings_are_reproducible(seeding: str) -> None:
    x, obs = _toy()
    out = [core.harmonize(x, obs, "batch", seeding=seeding, rng=3)[0] for _ in range(2)]
    np.testing.assert_array_equal(out[0], out[1])
    assert np.isfinite(out[0]).all()


def test_harmonize_default_seeding_follows_flavor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class SeededError(Exception):
        pass

    def afkmc2(*_args, **_kwargs):
        raise SeededError

    monkeypatch.setattr(core, "_seed_afkmc2", afkmc2)
    x, obs = _toy()
    with pytest.raises(SeededError):
        core.harmonize(x, obs, "batch", max_iter_harmony=1, rng=0)
    # harmony1 keeps the k-means initialisation
    core.harmonize(x, obs, "batch", flavor="harmony1", max_iter_harmony=1, rng=0)


def test_harmonize_edge_inputs() -> None:
    x, obs = _toy()
    for extra in (dict(tau=5), dict(stop_rule="paper_window3"), dict(n_clusters=2)):
        assert np.isfinite(core.harmonize(x, obs, "batch", rng=0, **extra)[0]).all()
    # a one-cell batch
    obs1 = obs.assign(batch=pd.Categorical(["a"] * (len(x) - 1) + ["b"]))
    assert np.isfinite(core.harmonize(x, obs1, "batch", rng=0)[0]).all()
    # fewer cells per block than chunks
    assert np.isfinite(core.harmonize(x[:40], obs[:40], "batch", rng=0)[0]).all()


@pytest.mark.skipif(
    numba.config.NUMBA_NUM_THREADS < 2, reason="needs at least two numba threads"
)
@pytest.mark.parametrize("dtype", _DTYPES)
def test_harmonize_does_not_depend_on_thread_count(dtype) -> None:
    x, obs = _toy(dtype=dtype)
    previous = numba.get_num_threads()
    try:
        out = []
        for n_threads in (1, 2):
            numba.set_num_threads(n_threads)
            # afkmc2 is bit-reproducible; scikit-learn's KMeans is not
            out.append(core.harmonize(x, obs, ["batch", "half"], rng=0))
    finally:
        numba.set_num_threads(previous)
    np.testing.assert_array_equal(out[0][0], out[1][0])
    assert out[0][1]["objectives"] == out[1][1]["objectives"]


@pytest.mark.parametrize("keys", ["batch", ["batch", "half"]])
def test_harmonize_low_memory_is_identical(keys) -> None:
    x, obs = _toy(n=1500)
    full = core.harmonize(x, obs, keys, rng=0)[0]
    lean = core.harmonize(x, obs, keys, low_memory=True, rng=0)[0]
    np.testing.assert_array_equal(full, lean)


def test_harmonize_leaves_out_all_zero_cells() -> None:
    x, obs = _toy()
    x0 = x.copy()
    x0[[3, 77]] = 0
    with pytest.warns(UserWarning, match="2 cells have an all-zero embedding"):
        out, _ = core.harmonize(x0, obs, "batch", rng=0)
    keep = np.ones(len(x), dtype=bool)
    keep[[3, 77]] = False
    np.testing.assert_array_equal(out[~keep], 0)
    # same as running without them
    expected, _ = core.harmonize(x[keep], obs[keep], "batch", rng=0)
    np.testing.assert_array_equal(out[keep], expected)


def test_harmonize_four_keys_and_a_level_emptied_by_the_zero_cell_filter() -> None:
    x, obs = _toy()
    obs = obs.assign(
        third=pd.Categorical(np.arange(len(x)) % 3),
        lone=pd.Categorical(["only"] + ["rest"] * (len(x) - 1)),
    )
    x = x.copy()
    x[0] = 0  # the one cell of level "only" is all-zero and gets filtered
    with pytest.warns(UserWarning, match="1 cells have an all-zero embedding"):
        out, _ = core.harmonize(x, obs, ["batch", "third", "half", "lone"], rng=0)
    assert np.isfinite(out).all()
    assert np.abs(out - x).max() > 0.1
    np.testing.assert_array_equal(out[0], 0)


def test_harmonize_tau_and_recompute_counts() -> None:
    x, obs = _toy()

    def run(**extra):
        return core.harmonize(x, obs, "batch", rng=0, **extra)[0]

    np.testing.assert_array_equal(run(tau=1e-9), run())  # discount factor 1
    np.testing.assert_allclose(run(tau=1e9), run(theta=0.0), atol=1e-4)  # factor 0
    assert not np.array_equal(run(tau=10), run())  # in between it changes the result
    incremental = run(recompute_counts=False)
    assert not np.array_equal(incremental, run())  # the flag takes effect ...
    np.testing.assert_allclose(incremental, run(), atol=1e-2)  # ... at rounding level


def test_harmonize_penalty_combine() -> None:
    x, obs = _toy()
    one_key = [
        core.harmonize(x, obs, "batch", penalty_combine=rule, rng=0)[0]
        for rule in ("product", "sum")
    ]
    np.testing.assert_array_equal(one_key[0], one_key[1])  # one key: nothing to combine
    two_keys = [
        core.harmonize(x, obs, ["batch", "half"], penalty_combine=rule, rng=0)[0]
        for rule in ("product", "sum")
    ]
    assert np.isfinite(two_keys[1]).all()
    assert not np.array_equal(two_keys[0], two_keys[1])
    # a key with theta = 0 is switched off under either rule, so with the other key alone
    # the sum equals the product
    off = [
        core.harmonize(
            x, obs, ["batch", "half"], theta=[2.0, 0.0], penalty_combine=rule, rng=0
        )[0]
        for rule in ("product", "sum")
    ]
    np.testing.assert_allclose(off[0], off[1], atol=1e-5)


def test_harmonize_seeding_options() -> None:
    x, obs = _toy()
    kw = dict(seeding="afkmc2", seeding_lloyd_steps=3, rng=0)
    plain = core.harmonize(x, obs, "batch", seeding="afkmc2", rng=0)[0]
    polished = core.harmonize(x, obs, "batch", **kw)[0]
    assert not np.array_equal(polished, plain)  # the polish takes effect
    # a cap above the number of cells changes nothing; a binding cap does
    uncapped = core.harmonize(x, obs, "batch", seeding_cells_per_cluster=10**6, **kw)[0]
    np.testing.assert_array_equal(uncapped, polished)
    capped = core.harmonize(x, obs, "batch", seeding_cells_per_cluster=10, **kw)[0]
    assert np.isfinite(capped).all()
    assert not np.array_equal(capped, polished)
    # subsampled KMeans
    km = core.harmonize(
        x, obs, "batch", seeding="kmeans", seeding_cells_per_cluster=10, rng=0
    )
    assert np.isfinite(km[0]).all()
    # without Lloyd steps nothing uses the subsample, so a cap changes nothing
    np.testing.assert_array_equal(
        core.harmonize(x, obs, "batch", seeding_cells_per_cluster=10, rng=0)[0],
        core.harmonize(x, obs, "batch", rng=0)[0],
    )


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (dict(alpha=float("nan")), "alpha"),
        (dict(flavor="harmony1", ridge_lambda=0.0), "ridge_lambda"),
        (dict(block_proportion=2.0), "block_proportion"),
        (dict(sigma=0.0), "sigma"),
        (dict(max_iter_harmony=0), "max_iter_harmony"),
        (dict(seeding="nope"), "seeding must be one of"),
        (dict(penalty_combine="mean"), "penalty_combine must be one of"),
        (dict(stop_rule="never"), "stop_rule must be one of"),
        (dict(flavor="harmony3"), "flavor must be one of"),
        (dict(n_clusters=0), "n_clusters must be between 1"),
        (dict(n_clusters=601), "n_clusters must be between 1"),
    ],
)
def test_harmonize_rejects_bad_parameters(kwargs: dict, match: str) -> None:
    x, obs = _toy()
    with pytest.raises(ValueError, match=match):
        core.harmonize(x, obs, "batch", **kwargs)


def test_harmonize_rejects_bad_inputs() -> None:
    x, obs = _toy()
    with pytest.raises(ValueError, match="rows but x has"):
        core.harmonize(x, obs[:500], "batch")
    with pytest.raises(ValueError, match="rows but x has"):
        core.harmonize(x[:500], obs, "batch")
    with pytest.raises(ValueError, match="2-dimensional"):
        core.harmonize(x[:, 0], obs, "batch")
    for dtype in (np.float16, np.int64):
        with pytest.raises(ValueError, match="dtype must be float32 or float64"):
            core.harmonize(x.astype(dtype), obs, "batch")
    x0 = x.copy()
    x0[0, 0] = np.inf
    with pytest.raises(ValueError, match="non-finite"):
        core.harmonize(x0, obs, "batch")
    with (
        pytest.raises(ValueError, match="non-zero embedding"),
        pytest.warns(UserWarning, match="all-zero embedding"),
    ):
        core.harmonize(np.zeros_like(x), obs, "batch")


def test_block_bounds_match_array_split() -> None:
    for n in (7, 40, 600, 1001):
        n_blocks = int(min(n, 1 // 0.05))
        np.testing.assert_array_equal(core._block_bounds(n, 0.05), _bounds(n, n_blocks))


def _closest_pair(y: np.ndarray) -> float:
    d = np.linalg.norm(y[:, None] - y[None], axis=2)
    return d[np.triu_indices(len(y), 1)].min()


# squared-distance weights (k-means++, AFK-MC²) spread seeds more than plain-distance weights
@pytest.mark.parametrize(
    ("seeder", "factor"),
    [("_seed_kmeanspp", 1.5), ("_seed_afkmc2", 1.5), ("_seed_paper_nearest", 1.2)],
)
def test_distance_weighted_seeders_spread_seeds(seeder: str, factor: float) -> None:
    _, z, *_ = _state(np.float32)
    gaps, uniform = [], []
    for seed in range(20):
        rng = np.random.default_rng(seed)
        y = getattr(core, seeder)(z, 10, rng)
        assert y.shape == (10, z.shape[1])
        assert y.dtype == np.float32
        np.testing.assert_allclose(np.linalg.norm(y, axis=1), 1, atol=1e-5)
        gaps.append(_closest_pair(y))
        uniform.append(_closest_pair(z[rng.choice(len(z), 10, replace=False)]))
    # further apart than uniformly drawn cells
    assert np.mean(gaps) > factor * np.mean(uniform)


def test_seeders_survive_identical_cells() -> None:
    z = np.tile(np.array([[1.0, 0.0, 0.0]], dtype=np.float32), (50, 1))
    for seeder in (core._seed_kmeanspp, core._seed_afkmc2, core._seed_paper_nearest):
        assert np.isfinite(seeder(z, 5, np.random.default_rng(0))).all()


@pytest.mark.parametrize("dtype", _DTYPES)
def test_lloyd_polish_lowers_the_kmeans_cost(dtype) -> None:
    _, z, *_ = _state(dtype)
    y = core._seed_afkmc2(z, _K, np.random.default_rng(0))

    def cost(c):
        return float(np.sum(1 - (z @ c.T).max(axis=1)))

    polished = core._lloyd_spherical(z, y, 3)
    assert cost(polished) < cost(y)
    assert cost(core._lloyd_spherical(z, polished, 1)) <= cost(polished) + 1e-4
    np.testing.assert_allclose(np.linalg.norm(polished, axis=1), 1, atol=1e-5)
    assert polished.dtype == dtype

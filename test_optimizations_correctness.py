#!/usr/bin/env python3
"""
Comprehensive Correctness Tests for HVG and Normalization Optimizations

This script verifies that optimized implementations produce EXACTLY the same
results as the standard implementations.
"""

from __future__ import annotations

import numpy as np
from optimized_hvg_implementation import (
    _compute_hvg_stats_dense,
    _compute_hvg_stats_sparse,
)
from scipy import sparse

import scanpy as sc


def create_test_dataset(
    n_obs: int = 1000, n_vars: int = 2000, density: float = 0.1, seed: int = 42
):
    """Create a realistic test dataset."""
    np.random.seed(seed)

    # Generate sparse count matrix with Poisson distribution
    X = sparse.random(n_obs, n_vars, density=density, format="csr", random_state=seed)
    X.data = np.random.poisson(X.data * 5 + 1).astype(np.float32)

    # Create AnnData object
    adata = sc.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_vars)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_obs)]

    return adata


def test_hvg_mean_variance_correctness():
    """Test that optimized mean/variance computation is exactly correct."""
    print("🧪 Testing HVG Mean/Variance Correctness")
    print("=" * 60)

    test_configs = [
        (500, 1000, 0.2, "Small Dense"),
        (2000, 3000, 0.05, "Medium Sparse"),
        (5000, 4000, 0.02, "Large Very Sparse"),
    ]

    all_passed = True

    for n_obs, n_vars, density, test_name in test_configs:
        print(
            f"\n📊 Testing {test_name}: {n_obs:,} × {n_vars:,} (density: {density:.3f})"
        )

        # Create test data
        adata = create_test_dataset(n_obs, n_vars, density)
        X_sparse = adata.X
        X_dense = X_sparse.toarray()

        actual_density = X_sparse.nnz / (n_obs * n_vars)
        print(f"   Actual density: {actual_density:.4f}")

        # Test optimized sparse implementation
        means_opt_sparse, vars_opt_sparse = _compute_hvg_stats_sparse(
            X_sparse.data, X_sparse.indices, X_sparse.indptr, n_obs, n_vars
        )

        # Test optimized dense implementation
        means_opt_dense, vars_opt_dense = _compute_hvg_stats_dense(X_dense)

        # Test against numpy reference
        means_numpy = np.mean(X_dense, axis=0, dtype=np.float64).astype(np.float32)
        vars_numpy = np.var(X_dense, axis=0, dtype=np.float64).astype(np.float32)

        # Test against scanpy stats if available
        try:
            from scanpy._utils import stats

            means_scanpy, vars_scanpy = stats.mean_var(X_sparse, axis=0)
            means_scanpy = means_scanpy.astype(np.float32)
            vars_scanpy = vars_scanpy.astype(np.float32)
            has_scanpy_stats = True
        except (ImportError, AttributeError):
            has_scanpy_stats = False

        # Check correctness
        print("   🔍 Checking correctness...")

        # Sparse vs numpy
        means_diff_sparse = np.max(np.abs(means_opt_sparse - means_numpy))
        vars_diff_sparse = np.max(np.abs(vars_opt_sparse - vars_numpy))

        sparse_means_ok = means_diff_sparse < 1e-6
        sparse_vars_ok = vars_diff_sparse < 1e-6

        print(
            f"      Sparse vs Numpy - Means: {'✅' if sparse_means_ok else '❌'} (max diff: {means_diff_sparse:.2e})"
        )
        print(
            f"      Sparse vs Numpy - Vars:  {'✅' if sparse_vars_ok else '❌'} (max diff: {vars_diff_sparse:.2e})"
        )

        # Dense vs numpy
        means_diff_dense = np.max(np.abs(means_opt_dense - means_numpy))
        vars_diff_dense = np.max(np.abs(vars_opt_dense - vars_numpy))

        dense_means_ok = means_diff_dense < 1e-6
        dense_vars_ok = vars_diff_dense < 1e-6

        print(
            f"      Dense vs Numpy - Means:  {'✅' if dense_means_ok else '❌'} (max diff: {means_diff_dense:.2e})"
        )
        print(
            f"      Dense vs Numpy - Vars:   {'✅' if dense_vars_ok else '❌'} (max diff: {vars_diff_dense:.2e})"
        )

        # Sparse vs dense
        means_diff_sd = np.max(np.abs(means_opt_sparse - means_opt_dense))
        vars_diff_sd = np.max(np.abs(vars_opt_sparse - vars_opt_dense))

        sd_means_ok = means_diff_sd < 1e-6
        sd_vars_ok = vars_diff_sd < 1e-6

        print(
            f"      Sparse vs Dense - Means: {'✅' if sd_means_ok else '❌'} (max diff: {means_diff_sd:.2e})"
        )
        print(
            f"      Sparse vs Dense - Vars:  {'✅' if sd_vars_ok else '❌'} (max diff: {vars_diff_sd:.2e})"
        )

        # Scanpy stats comparison if available
        if has_scanpy_stats:
            means_diff_scanpy = np.max(np.abs(means_opt_sparse - means_scanpy))
            vars_diff_scanpy = np.max(np.abs(vars_opt_sparse - vars_scanpy))

            scanpy_means_ok = means_diff_scanpy < 1e-6
            scanpy_vars_ok = vars_diff_scanpy < 1e-6

            print(
                f"      Sparse vs Scanpy - Means: {'✅' if scanpy_means_ok else '❌'} (max diff: {means_diff_scanpy:.2e})"
            )
            print(
                f"      Sparse vs Scanpy - Vars:  {'✅' if scanpy_vars_ok else '❌'} (max diff: {vars_diff_scanpy:.2e})"
            )

        # Overall result for this test
        test_passed = (
            sparse_means_ok
            and sparse_vars_ok
            and dense_means_ok
            and dense_vars_ok
            and sd_means_ok
            and sd_vars_ok
        )

        if has_scanpy_stats:
            test_passed = test_passed and scanpy_means_ok and scanpy_vars_ok

        print(f"   🎯 {test_name}: {'✅ PASS' if test_passed else '❌ FAIL'}")

        if not test_passed:
            all_passed = False

    return all_passed


def test_normalization_correctness():
    """Test that optimized normalization produces exactly the same results."""
    print("\n🧪 Testing Normalization Correctness")
    print("=" * 60)

    # Check if optimization functions are available
    try:
        from scanpy.preprocessing._normalization_optimized import (
            _normalize_csr_optimized,
            apply_row_normalization_optimized,
            get_optimization_mode,
            set_optimization_mode,
        )

        optimization_available = True
    except ImportError:
        print("❌ Normalization optimization functions not available")
        return False

    test_configs = [
        (1000, 500, 0.1, "Small"),
        (5000, 2000, 0.05, "Medium"),
        (10000, 3000, 0.02, "Large"),
    ]

    all_passed = True

    for n_obs, n_vars, density, test_name in test_configs:
        print(
            f"\n📊 Testing {test_name}: {n_obs:,} × {n_vars:,} (density: {density:.3f})"
        )

        # Create test data
        np.random.seed(42)
        X = sparse.random(n_obs, n_vars, density=density, format="csr")
        X.data = np.random.poisson(X.data * 10 + 1).astype(np.float64)

        actual_density = X.nnz / (n_obs * n_vars)
        print(f"   Actual density: {actual_density:.4f}")

        # Test different optimization modes
        modes_to_test = ["disabled", "naive", "numba", "auto"]
        results = {}

        for mode in modes_to_test:
            print(f"   🔄 Testing mode: {mode}")

            # Set optimization mode
            set_optimization_mode(mode)

            # Create copy of data
            X_copy = X.copy()

            # Apply normalization
            try:
                if mode == "disabled":
                    # Use standard scanpy normalization (modifies in-place)
                    adata_temp = sc.AnnData(X_copy)
                    sc.pp.normalize_total(adata_temp, target_sum=1e4)
                    result = adata_temp.X  # The normalized matrix
                else:
                    # Use optimized normalization (modifies in-place)
                    result_tuple = apply_row_normalization_optimized(
                        X_copy, target_sum=1e4
                    )
                    if result_tuple is None:
                        # Disabled mode, use standard scanpy
                        adata_temp = sc.AnnData(X_copy)
                        sc.pp.normalize_total(adata_temp, target_sum=1e4)
                        result = adata_temp.X
                    else:
                        # X_copy is modified in-place, so it's now the normalized matrix
                        result = X_copy

                results[mode] = result
                print("      ✅ Success")

            except Exception as e:
                print(f"      ❌ Failed: {e}")
                all_passed = False
                continue

        # Compare all modes against disabled (reference)
        if "disabled" in results:
            reference = results["disabled"]
            print("   🔍 Checking correctness against reference (disabled mode)...")

            for mode in ["naive", "numba", "auto"]:
                if mode in results:
                    try:
                        # Convert both to same format for comparison
                        ref_dense = (
                            reference.toarray()
                            if sparse.issparse(reference)
                            else reference
                        )
                        result_dense = (
                            results[mode].toarray()
                            if sparse.issparse(results[mode])
                            else results[mode]
                        )

                        max_diff = np.max(np.abs(result_dense - ref_dense))
                        is_correct = (
                            max_diff < 1e-9
                        )  # Allow for numerical precision differences

                        print(
                            f"      {mode} vs disabled: {'✅' if is_correct else '❌'} (max diff: {max_diff:.2e})"
                        )

                        if not is_correct:
                            all_passed = False
                    except Exception as e:
                        print(f"      {mode} vs disabled: ❌ (comparison error: {e})")
                        all_passed = False

        # Track test status for this configuration
        test_passed_for_config = True
        if "disabled" in results:
            for mode in ["naive", "numba", "auto"]:
                if mode in results:
                    try:
                        ref_dense = (
                            results["disabled"].toarray()
                            if sparse.issparse(results["disabled"])
                            else results["disabled"]
                        )
                        result_dense = (
                            results[mode].toarray()
                            if sparse.issparse(results[mode])
                            else results[mode]
                        )
                        max_diff = np.max(np.abs(result_dense - ref_dense))
                        if max_diff >= 1e-9:
                            test_passed_for_config = False
                            break
                    except Exception:
                        test_passed_for_config = False
                        break

        print(
            f"   🎯 {test_name}: {'✅ PASS' if test_passed_for_config else '❌ FAIL'}"
        )

        if not test_passed_for_config:
            all_passed = False

    # Reset to disabled mode
    set_optimization_mode("disabled")

    return all_passed


def test_full_hvg_pipeline_correctness():
    """Test the full HVG pipeline for basic mean/variance consistency."""
    print("\n🧪 Testing Full HVG Pipeline Basic Consistency")
    print("=" * 60)

    # Create test dataset
    adata = create_test_dataset(2000, 3000, 0.05)
    adata_copy = adata.copy()

    print(f"📊 Dataset: {adata.shape[0]:,} × {adata.shape[1]:,}")
    print(f"   Density: {adata.X.nnz / (adata.n_obs * adata.n_vars):.4f}")

    # Test our optimized mean/variance computation
    print("\n⚡ Testing optimized mean/variance computation...")
    X = adata.X
    means_opt, vars_opt = _compute_hvg_stats_sparse(
        X.data, X.indices, X.indptr, X.shape[0], X.shape[1]
    )

    # Compare with numpy reference
    X_dense = X.toarray()
    means_ref = np.mean(X_dense, axis=0, dtype=np.float64).astype(np.float32)
    vars_ref = np.var(X_dense, axis=0, dtype=np.float64).astype(np.float32)

    means_diff = np.max(np.abs(means_opt - means_ref))
    vars_diff = np.max(np.abs(vars_opt - vars_ref))

    means_ok = means_diff < 1e-10
    vars_ok = vars_diff < 1e-10

    print(
        f"   Means correctness: {'✅' if means_ok else '❌'} (max diff: {means_diff:.2e})"
    )
    print(
        f"   Variances correctness: {'✅' if vars_ok else '❌'} (max diff: {vars_diff:.2e})"
    )

    pipeline_ok = means_ok and vars_ok
    print(
        f"\n🎯 HVG Pipeline Basic Consistency: {'✅ PASS' if pipeline_ok else '❌ FAIL'}"
    )

    return pipeline_ok


def main():
    """Run all correctness tests."""
    print("🔬 Comprehensive Optimization Correctness Tests")
    print("=" * 80)

    # Test HVG mean/variance correctness
    hvg_ok = test_hvg_mean_variance_correctness()

    # Test normalization correctness
    norm_ok = test_normalization_correctness()

    # Test full HVG pipeline basic consistency
    pipeline_ok = test_full_hvg_pipeline_correctness()

    # Overall summary
    print("\n" + "=" * 80)
    print("📋 FINAL CORRECTNESS SUMMARY")
    print("=" * 80)
    print(f"HVG Mean/Variance:     {'✅ PASS' if hvg_ok else '❌ FAIL'}")
    print(f"Normalization:         {'✅ PASS' if norm_ok else '❌ FAIL'}")
    print(f"HVG Pipeline Basic:    {'✅ PASS' if pipeline_ok else '❌ FAIL'}")

    all_passed = hvg_ok and norm_ok and pipeline_ok
    print(
        f"\n🎯 OVERALL CORRECTNESS: {'✅ ALL TESTS PASS' if all_passed else '❌ SOME TESTS FAILED'}"
    )

    if all_passed:
        print("\n🎉 All optimizations produce mathematically equivalent results!")
        print("   The optimized implementations are safe to use in production.")
    else:
        print("\n⚠️  Some correctness issues found. Please review the failing tests.")

    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())

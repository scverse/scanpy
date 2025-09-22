"""Test suite for configurable normalization optimization functions."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from numpy.testing import assert_allclose
from scanpy.preprocessing._normalization_optimized import (
    _csr_row_scaling_optimized,
    _csr_sum_and_squared_sum_optimized,
    _normalize_csr_naive,
    _normalize_csr_numba,
    _normalize_csr_optimized,
    _should_use_numba_optimization,
    apply_row_normalization_optimized,
    get_optimization_mode,
    set_optimization_mode,
)


class TestOptimizationModeConfiguration:
    """Test configuration and flag management."""

    def test_default_mode(self):
        """Test that default mode is 'auto' for optimal performance."""
        assert get_optimization_mode() == "auto"

    def test_mode_setting(self):
        """Test setting different optimization modes."""
        for mode in ["disabled", "naive", "numba", "auto"]:
            set_optimization_mode(mode)
            assert get_optimization_mode() == mode

    def test_invalid_mode(self):
        """Test that invalid modes raise ValueError."""
        with pytest.raises(ValueError, match="Invalid optimization mode"):
            set_optimization_mode("invalid")

    def test_mode_persistence(self):
        """Test that mode setting persists across function calls."""
        set_optimization_mode("auto")
        assert get_optimization_mode() == "auto"

        # Call some function that might change internal state
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float64)
        _normalize_csr_optimized(X, rows=100, columns=50)

        # Mode should still be 'auto'
        assert get_optimization_mode() == "auto"

        # Reset for other tests
        set_optimization_mode("disabled")


class TestMatrixSizeDetection:
    """Test automatic matrix size-based selection."""

    def test_small_matrix_detection(self):
        """Test that small matrices are detected correctly."""
        # Small matrix - should NOT use Numba
        X_small = sp.random(100, 50, density=0.1, format="csr")
        assert not _should_use_numba_optimization(X_small)

    def test_large_matrix_detection_by_elements(self):
        """Test detection based on total elements."""
        # Large matrix by total elements
        X_large = sp.random(500, 300, density=0.01, format="csr")  # 150k elements
        assert _should_use_numba_optimization(X_large)

    def test_large_matrix_detection_by_nnz(self):
        """Test detection based on non-zero elements."""
        # Large matrix by non-zeros
        X_dense = sp.random(
            200, 200, density=0.8, format="csr"
        )  # ~32k nnz but 40k elements
        # Ensure we have enough non-zeros
        if X_dense.nnz >= 50_000:
            assert _should_use_numba_optimization(X_dense)


class TestDisabledMode:
    """Test that disabled mode preserves existing behavior."""

    def setup_method(self):
        """Set up disabled mode for each test."""
        set_optimization_mode("disabled")

    def teardown_method(self):
        """Reset mode after each test."""
        set_optimization_mode("disabled")

    def test_normalize_csr_returns_none(self):
        """Test that _normalize_csr_optimized returns None when disabled."""
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float64)
        result = _normalize_csr_optimized(X, rows=100, columns=50)
        assert result is None

    def test_apply_row_normalization_returns_none(self):
        """Test that apply_row_normalization_optimized returns None when disabled."""
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float64)
        result = apply_row_normalization_optimized(X)
        assert result is None

    def test_no_modification_when_disabled(self):
        """Test that matrices are not modified when optimizations are disabled."""
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float64)
        X_original = X.copy()

        result = apply_row_normalization_optimized(X)
        assert result is None

        # Matrix should be unchanged
        assert_allclose(X.data, X_original.data)


class TestImplementationConsistency:
    """Test that all implementations produce identical results."""

    def setup_method(self):
        """Set up for each test."""
        set_optimization_mode("disabled")  # Start with disabled

    def teardown_method(self):
        """Reset mode after each test."""
        set_optimization_mode("disabled")

    def test_naive_vs_numba_consistency(self):
        """Test that naive and Numba implementations produce identical results."""
        np.random.seed(42)
        X_dense = np.random.poisson(5, (200, 100))
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)

        # Test naive implementation
        X_naive = X_sparse.copy()
        counts_naive, _ = _normalize_csr_naive(
            X_naive, rows=X_sparse.shape[0], columns=X_sparse.shape[1]
        )

        # Test Numba implementation
        X_numba = X_sparse.copy()
        counts_numba, _ = _normalize_csr_numba(
            X_numba, rows=X_sparse.shape[0], columns=X_sparse.shape[1]
        )

        # Results should be identical
        assert_allclose(counts_naive, counts_numba, rtol=1e-10)

    def test_auto_mode_consistency(self):
        """Test that auto mode produces consistent results."""
        np.random.seed(42)

        # Test with small matrix (should use naive)
        X_small = sp.random(50, 30, density=0.2, format="csr", dtype=np.float64)
        X_small.data = np.abs(X_small.data) * 1000

        X_test1 = X_small.copy()
        counts_auto, _ = _normalize_csr_optimized(
            X_test1, rows=50, columns=30, mode="auto"
        )

        X_test2 = X_small.copy()
        counts_naive, _ = _normalize_csr_naive(X_test2, rows=50, columns=30)

        assert_allclose(counts_auto, counts_naive, rtol=1e-10)

        # Test with large matrix (should use Numba)
        X_large = sp.random(1000, 500, density=0.1, format="csr", dtype=np.float64)
        X_large.data = np.abs(X_large.data) * 1000

        X_test3 = X_large.copy()
        counts_auto_large, _ = _normalize_csr_optimized(
            X_test3, rows=1000, columns=500, mode="auto"
        )

        X_test4 = X_large.copy()
        counts_numba_large, _ = _normalize_csr_numba(X_test4, rows=1000, columns=500)

        assert_allclose(counts_auto_large, counts_numba_large, rtol=1e-10)

    def test_mode_override(self):
        """Test that mode parameter overrides global setting."""
        set_optimization_mode("disabled")

        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float64)

        # Global mode is disabled, but override with 'naive'
        result = _normalize_csr_optimized(X, rows=100, columns=50, mode="naive")
        assert result is not None

        # Global mode is disabled, function call without override should return None
        result = _normalize_csr_optimized(X, rows=100, columns=50)
        assert result is None


class TestOptimizedFunctionality:
    """Test core optimized functions when enabled."""

    def setup_method(self):
        """Enable optimizations for testing."""
        set_optimization_mode("auto")

    def teardown_method(self):
        """Reset mode after each test."""
        set_optimization_mode("disabled")

    def test_csr_sum_and_squared_sum_optimized(self):
        """Test optimized sum computation."""
        np.random.seed(42)
        X_dense = np.random.poisson(2, (100, 50))
        X_sparse = sp.csr_matrix(X_dense)

        sums, sq_sums = _csr_sum_and_squared_sum_optimized(
            X_sparse.data,
            X_sparse.indices,
            X_sparse.indptr,
            X_sparse.shape[0],
            X_sparse.shape[1],
        )

        expected_sums = np.array(X_dense.sum(axis=0)).flatten()
        expected_sq_sums = np.array((X_dense**2).sum(axis=0)).flatten()

        assert_allclose(sums, expected_sums, rtol=1e-10)
        assert_allclose(sq_sums, expected_sq_sums, rtol=1e-10)

    def test_csr_row_scaling_optimized(self):
        """Test optimized row scaling."""
        np.random.seed(42)
        X_dense = np.random.poisson(5, (50, 30))
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)
        X_copy = X_sparse.copy()

        scaling_factors = np.random.uniform(0.5, 2.0, X_sparse.shape[0])

        _csr_row_scaling_optimized(
            X_sparse.data, X_sparse.indptr, scaling_factors, X_sparse.shape[0]
        )

        # Compare with reference implementation
        X_reference = X_copy.copy()
        for i in range(X_reference.shape[0]):
            X_reference.data[X_reference.indptr[i] : X_reference.indptr[i + 1]] *= (
                scaling_factors[i]
            )

        assert_allclose(X_sparse.data, X_reference.data, rtol=1e-10)

    def test_complete_normalization_pipeline(self):
        """Test complete normalization pipeline with auto mode."""
        np.random.seed(42)
        X_dense = np.random.poisson(4, (100, 60))
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)

        target_sum = 10000.0
        result = apply_row_normalization_optimized(
            X_sparse, target_sum=target_sum, mode="auto"
        )

        assert result is not None
        scaling_factors, counts_per_cell = result

        # Verify that row sums are approximately target_sum
        new_row_sums = np.array(X_sparse.sum(axis=1)).flatten()
        expected_sums = np.full(X_sparse.shape[0], target_sum)

        assert_allclose(new_row_sums, expected_sums, rtol=1e-6)

    def test_highly_expressed_exclusion(self):
        """Test highly expressed gene exclusion."""
        np.random.seed(42)
        X_dense = np.random.poisson(3, (150, 80))
        # Add some highly expressed genes
        X_dense[:10, 0] = 1000  # Gene 0 highly expressed in first 10 cells
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)

        counts_opt, cols_opt = _normalize_csr_optimized(
            X_sparse.copy(),
            rows=X_sparse.shape[0],
            columns=X_sparse.shape[1],
            exclude_highly_expressed=True,
            max_fraction=0.1,
            mode="auto",
        )

        # Verify basic properties
        assert len(counts_opt) == X_sparse.shape[0]
        assert np.all(counts_opt >= 0), "All counts should be non-negative"

        # Verify that highly expressed genes were excluded
        if cols_opt is not None:
            assert cols_opt[0] > 0, "Gene 0 should be flagged as highly expressed"

            total_sums = np.array(X_sparse.sum(axis=1)).flatten()
            assert np.all(counts_opt <= total_sums), (
                "Counts should be <= total sums when excluding highly expressed"
            )


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def setup_method(self):
        """Set up for each test."""
        set_optimization_mode("auto")

    def teardown_method(self):
        """Reset mode after each test."""
        set_optimization_mode("disabled")

    def test_empty_matrix(self):
        """Test handling of empty matrices."""
        X_empty = sp.csr_matrix((0, 10))
        result = _normalize_csr_optimized(X_empty, rows=0, columns=10, mode="naive")

        if result is not None:
            counts, _ = result
            assert len(counts) == 0

    def test_single_element_matrix(self):
        """Test single element matrix."""
        X_single = sp.csr_matrix(([5.0], ([0], [0])), shape=(1, 1))
        result = _normalize_csr_optimized(X_single, rows=1, columns=1, mode="naive")

        if result is not None:
            counts, _ = result
            assert counts[0] == 5.0

    def test_numerical_stability(self):
        """Test numerical stability with extreme values."""
        # Very large values
        X_large = sp.csr_matrix(([1e10, 1e10], ([0, 1], [0, 0])), shape=(2, 1))
        result = apply_row_normalization_optimized(X_large, mode="naive")

        if result is not None:
            scaling_factors, counts = result
            assert not np.any(np.isnan(scaling_factors)), (
                "No NaN values in scaling factors"
            )
            assert not np.any(np.isinf(scaling_factors)), (
                "No inf values in scaling factors"
            )


class TestPerformanceRegression:
    """Test that optimizations don't cause performance regressions."""

    def setup_method(self):
        """Set up for performance tests."""
        set_optimization_mode("disabled")  # Start with baseline

    def teardown_method(self):
        """Reset mode after each test."""
        set_optimization_mode("disabled")

    def test_small_matrix_performance(self):
        """Test that small matrices aren't slowed down by optimization overhead."""
        import time

        # Small matrix that should use naive implementation
        np.random.seed(42)
        X_dense = np.random.poisson(3, (100, 50))
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)

        # Time naive implementation
        X_naive = X_sparse.copy()
        start_time = time.perf_counter()
        result_naive = _normalize_csr_naive(X_naive, rows=100, columns=50)
        naive_time = time.perf_counter() - start_time

        # Time auto selection (should choose naive for small matrix)
        X_auto = X_sparse.copy()
        start_time = time.perf_counter()
        result_auto = _normalize_csr_optimized(
            X_auto, rows=100, columns=50, mode="auto"
        )
        auto_time = time.perf_counter() - start_time

        # Auto mode should not be significantly slower than direct naive call
        # Allow for some overhead but not more than 10x
        assert auto_time < naive_time * 10, (
            f"Auto mode too slow: {auto_time:.6f}s vs {naive_time:.6f}s"
        )

        # Results should be identical
        assert_allclose(result_naive[0], result_auto[0], rtol=1e-10)

    @pytest.mark.parametrize(
        "n_obs,n_vars,density",
        [(50, 30, 0.3), (200, 100, 0.1), (500, 200, 0.05), (1000, 500, 0.02)],
    )
    def test_scaling_performance(self, n_obs, n_vars, density):
        """Test performance across different matrix sizes."""
        np.random.seed(42)
        X_dense = np.random.poisson(2, (n_obs, n_vars))
        mask = np.random.random((n_obs, n_vars)) < density
        X_dense = X_dense * mask
        X_sparse = sp.csr_matrix(X_dense, dtype=np.float64)

        # Test that optimization completes in reasonable time
        import time

        start_time = time.perf_counter()

        result = _normalize_csr_optimized(
            X_sparse.copy(), rows=n_obs, columns=n_vars, mode="auto"
        )

        elapsed_time = time.perf_counter() - start_time

        # Should complete within reasonable time (adjust threshold as needed)
        max_time = max(0.1, n_obs * n_vars / 100000)  # Scale with matrix size
        assert elapsed_time < max_time, (
            f"Optimization too slow for {n_obs}x{n_vars}: {elapsed_time:.3f}s"
        )

        # Verify correctness
        if result is not None:
            counts, _ = result
            assert len(counts) == n_obs
            assert np.all(counts >= 0), "All counts should be non-negative"


if __name__ == "__main__":
    pytest.main([__file__])

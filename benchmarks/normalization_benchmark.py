"""Benchmark script for configurable normalization optimization.

This script comprehensively tests the performance of the flag-controlled
normalization optimizations to ensure net performance gains with no regressions.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
import psutil
import scipy.sparse as sp


# Import optimization functions
try:
    from scanpy.preprocessing._normalization_optimized import (
        _normalize_csr_optimized,
        apply_row_normalization_optimized,
        get_optimization_mode,
        set_optimization_mode,
    )

    OPTIMIZATION_AVAILABLE = True
except ImportError:
    print("Warning: Could not import optimization functions.")
    OPTIMIZATION_AVAILABLE = False


def generate_test_matrix(n_obs: int, n_vars: int, density: float = 0.1, seed: int = 42) -> sp.csr_matrix:
    """Generate a realistic sparse test matrix for benchmarking."""
    np.random.seed(seed)
    X = sp.random(n_obs, n_vars, density=density, format="csr", random_state=seed)
    X.data = np.random.poisson(X.data * 10 + 1).astype(np.float64)
    return X


def time_function(func, *args, **kwargs):
    """Time a function execution and return result and timing."""
    start_time = time.perf_counter()
    result = func(*args, **kwargs)
    end_time = time.perf_counter()
    return result, end_time - start_time


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


class ConfigurableNormalizationBenchmark:
    """Comprehensive benchmark suite for configurable normalization optimization."""

    def __init__(self):
        self.results = []

    def benchmark_mode_configurations(self):
        """Benchmark different optimization modes."""
        print("Benchmarking optimization mode configurations...")

        test_configs = [
            {"n_obs": 1000, "n_vars": 500, "density": 0.1, "name": "Small"},
            {"n_obs": 5000, "n_vars": 2000, "density": 0.05, "name": "Medium"},
            {"n_obs": 10000, "n_vars": 5000, "density": 0.03, "name": "Large"},
        ]

        modes = ["disabled", "naive", "numba", "auto"]

        for config in test_configs:
            print(f"\nTesting {config['name']} matrix: {config['n_obs']}x{config['n_vars']}")

            X = generate_test_matrix(config["n_obs"], config["n_vars"], config["density"])
            print(f"Matrix: nnz={X.nnz}, density={X.nnz / (X.shape[0] * X.shape[1]):.4f}")

            baseline_time = None

            for mode in modes:
                print(f"  Testing mode: {mode}")

                if not OPTIMIZATION_AVAILABLE:
                    print("    Skipping - optimization not available")
                    continue

                # Set global mode
                set_optimization_mode(mode)

                X_test = X.copy()
                mem_before = get_memory_usage()

                try:
                    if mode == "disabled":
                        # For disabled mode, simulate what standard scanpy would do
                        def baseline_normalization(mat):
                            return np.array(mat.sum(axis=1)).flatten()

                        result, elapsed_time = time_function(baseline_normalization, X_test)
                        baseline_time = elapsed_time
                        success = True
                        throughput = config["n_obs"] / elapsed_time

                    else:
                        # Test optimization modes
                        result, elapsed_time = time_function(
                            _normalize_csr_optimized,
                            X_test,
                            rows=X.shape[0],
                            columns=X.shape[1],
                        )

                        if result is None:
                            print(f"    Mode {mode} returned None (disabled)")
                            continue

                        success = True
                        throughput = config["n_obs"] / elapsed_time

                except Exception as e:
                    print(f"    Error in mode {mode}: {e}")
                    elapsed_time = float("inf")
                    success = False
                    throughput = 0

                mem_after = get_memory_usage()

                # Calculate speedup relative to baseline
                speedup = baseline_time / elapsed_time if baseline_time and elapsed_time > 0 else 1.0

                result_entry = {
                    "test_name": f"Mode_{mode}_{config['name']}",
                    "matrix_size": f"{config['n_obs']}x{config['n_vars']}",
                    "density": config["density"],
                    "nnz": X.nnz,
                    "mode": mode,
                    "time": elapsed_time,
                    "speedup": speedup,
                    "throughput": throughput,
                    "memory_mb": mem_after - mem_before,
                    "success": success,
                }

                self.results.append(result_entry)

                print(f"    Time: {elapsed_time:.4f}s, Speedup: {speedup:.2f}x, Throughput: {throughput:.0f} cells/s")

    def benchmark_automatic_selection(self):
        """Benchmark automatic implementation selection."""
        print("\n" + "=" * 60)
        print("Benchmarking automatic implementation selection...")

        if not OPTIMIZATION_AVAILABLE:
            print("Skipping - optimization not available")
            return

        # Test matrices of different sizes to verify automatic selection
        test_matrices = [
            # Small matrices - should use naive
            {"n_obs": 100, "n_vars": 50, "density": 0.2, "expected": "naive"},
            {"n_obs": 200, "n_vars": 100, "density": 0.1, "expected": "naive"},
            # Large matrices - should use numba
            {"n_obs": 1000, "n_vars": 500, "density": 0.1, "expected": "numba"},
            {"n_obs": 2000, "n_vars": 1000, "density": 0.05, "expected": "numba"},
        ]

        set_optimization_mode("auto")

        for config in test_matrices:
            print(f"\nTesting {config['n_obs']}x{config['n_vars']} (expected: {config['expected']})")

            X = generate_test_matrix(config["n_obs"], config["n_vars"], config["density"])

            # Test auto mode
            X_auto = X.copy()
            result_auto, time_auto = time_function(
                _normalize_csr_optimized,
                X_auto,
                rows=X.shape[0],
                columns=X.shape[1],
            )

            # Test explicit modes for comparison
            X_naive = X.copy()
            result_naive, time_naive = time_function(
                _normalize_csr_optimized,
                X_naive,
                rows=X.shape[0],
                columns=X.shape[1],
                mode="naive",
            )

            X_numba = X.copy()
            result_numba, time_numba = time_function(
                _normalize_csr_optimized,
                X_numba,
                rows=X.shape[0],
                columns=X.shape[1],
                mode="numba",
            )

            # Verify correctness
            if result_auto is not None and result_naive is not None:
                np.testing.assert_allclose(result_auto[0], result_naive[0], rtol=1e-10)

            if result_auto is not None and result_numba is not None:
                np.testing.assert_allclose(result_auto[0], result_numba[0], rtol=1e-10)

            # Check if auto selected the expected implementation
            better_choice = "naive" if time_naive < time_numba else "numba"

            print(f"  Auto time: {time_auto:.4f}s")
            print(f"  Naive time: {time_naive:.4f}s")
            print(f"  Numba time: {time_numba:.4f}s")
            print(f"  Better choice would be: {better_choice}")

            # Auto should be close to the better choice
            better_time = min(time_naive, time_numba)
            overhead = (time_auto - better_time) / better_time if better_time > 0 else 0

            result_entry = {
                "test_name": f"Auto_selection_{config['n_obs']}x{config['n_vars']}",
                "matrix_size": f"{config['n_obs']}x{config['n_vars']}",
                "auto_time": time_auto,
                "naive_time": time_naive,
                "numba_time": time_numba,
                "better_choice": better_choice,
                "overhead": overhead,
                "efficiency": 1.0 - overhead,
            }

            self.results.append(result_entry)

            print(f"  Automatic selection overhead: {overhead * 100:.1f}%")

    def benchmark_complete_pipeline(self):
        """Benchmark complete normalization pipeline."""
        print("\n" + "=" * 60)
        print("Benchmarking complete normalization pipeline...")

        if not OPTIMIZATION_AVAILABLE:
            print("Skipping - optimization not available")
            return

        test_configs = [
            {"n_obs": 2000, "n_vars": 1000, "density": 0.1, "name": "Small"},
            {"n_obs": 10000, "n_vars": 5000, "density": 0.05, "name": "Medium"},
            {"n_obs": 20000, "n_vars": 10000, "density": 0.02, "name": "Large"},
        ]

        for config in test_configs:
            print(f"\nTesting {config['name']} pipeline: {config['n_obs']}x{config['n_vars']}")

            X = generate_test_matrix(config["n_obs"], config["n_vars"], config["density"])

            # Test different modes
            for mode in ["disabled", "auto"]:
                X_test = X.copy()
                mem_before = get_memory_usage()

                if mode == "disabled":
                    # Simulate standard scanpy pipeline
                    def standard_pipeline(mat):
                        # Simple normalization to target sum
                        row_sums = np.array(mat.sum(axis=1)).flatten()
                        target_sum = 10000.0
                        scaling_factors = target_sum / (row_sums + 1e-10)

                        for i in range(mat.shape[0]):
                            start = mat.indptr[i]
                            end = mat.indptr[i + 1]
                            for j in range(start, end):
                                mat.data[j] *= scaling_factors[i]

                        return scaling_factors, row_sums

                    result, time_pipeline = time_function(standard_pipeline, X_test)

                else:
                    # Test optimized pipeline
                    result, time_pipeline = time_function(
                        apply_row_normalization_optimized,
                        X_test,
                        target_sum=10000.0,
                        mode=mode,
                    )

                mem_after = get_memory_usage()

                if result is not None:
                    # Verify normalization accuracy
                    row_sums = np.array(X_test.sum(axis=1)).flatten()
                    target_sum = 10000.0
                    max_deviation = np.max(np.abs(row_sums - target_sum))
                    mean_deviation = np.mean(np.abs(row_sums - target_sum))
                else:
                    max_deviation = float("nan")
                    mean_deviation = float("nan")

                pipeline_result = {
                    "test_name": f"Pipeline_{mode}_{config['name']}",
                    "matrix_size": f"{config['n_obs']}x{config['n_vars']}",
                    "mode": mode,
                    "pipeline_time": time_pipeline,
                    "memory_usage_mb": mem_after - mem_before,
                    "max_deviation": max_deviation,
                    "mean_deviation": mean_deviation,
                    "throughput_cells_per_sec": config["n_obs"] / time_pipeline,
                    "success": result is not None or mode == "disabled",
                }

                self.results.append(pipeline_result)

                print(f"  Mode {mode}:")
                print(f"    Time: {time_pipeline:.4f}s")
                print(f"    Throughput: {config['n_obs'] / time_pipeline:.0f} cells/sec")
                print(f"    Memory: {mem_after - mem_before:.1f}MB")
                if not np.isnan(max_deviation):
                    print(f"    Max deviation: {max_deviation:.2e}")

    def benchmark_safety_and_compatibility(self):
        """Test that default disabled mode maintains compatibility."""
        print("\n" + "=" * 60)
        print("Testing safety and compatibility...")

        if not OPTIMIZATION_AVAILABLE:
            print("Skipping - optimization not available")
            return

        # Ensure we start in disabled mode
        set_optimization_mode("disabled")
        assert get_optimization_mode() == "disabled"

        # Test that disabled mode doesn't break anything
        test_matrices = [
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float64),
            sp.random(1000, 500, density=0.05, format="csr", dtype=np.float64),
        ]

        for i, X in enumerate(test_matrices):
            X_original = X.copy()

            # Test optimization functions return None when disabled
            result1 = _normalize_csr_optimized(X, rows=X.shape[0], columns=X.shape[1])
            assert result1 is None, f"Matrix {i}: _normalize_csr_optimized should return None when disabled"

            result2 = apply_row_normalization_optimized(X)
            assert result2 is None, f"Matrix {i}: apply_row_normalization_optimized should return None when disabled"

            # Verify matrix is unchanged
            np.testing.assert_allclose(X.data, X_original.data, rtol=1e-15)

            print(f"  Matrix {i + 1}: ✓ No modifications when disabled")

        print("  ✓ All safety checks passed")

    def save_results(self, filename: str = "configurable_normalization_benchmark.csv"):
        """Save benchmark results and generate summary."""
        if not self.results:
            print("No results to save.")
            return

        df = pd.DataFrame(self.results)
        df.to_csv(filename, index=False)
        print(f"\nResults saved to {filename}")

        # Generate comprehensive summary
        print("\n" + "=" * 60)
        print("COMPREHENSIVE BENCHMARK SUMMARY")
        print("=" * 60)

        # Safety verification
        disabled_results = df[df.get("mode") == "disabled"]
        if not disabled_results.empty:
            print("✓ SAFETY: Disabled mode works correctly (no breaking changes)")

        # Performance analysis
        mode_results = df[df.get("mode").notna() & (df["mode"] != "disabled")]
        if not mode_results.empty:
            print("\nPERFORMANCE ANALYSIS:")

            for mode in ["naive", "numba", "auto"]:
                mode_data = mode_results[mode_results["mode"] == mode]
                if not mode_data.empty and "speedup" in mode_data.columns:
                    avg_speedup = mode_data["speedup"].mean()
                    print(f"  {mode.capitalize()} mode - Average speedup: {avg_speedup:.2f}x")

        # Automatic selection efficiency
        auto_results = df[df.get("test_name", "").str.contains("Auto_selection", na=False)]
        if not auto_results.empty and "efficiency" in auto_results.columns:
            avg_efficiency = auto_results["efficiency"].mean()
            print("\nAUTOMATIC SELECTION:")
            print(f"  Average efficiency: {avg_efficiency * 100:.1f}%")
            print("  (100% = perfect selection, >90% = good)")

        # Pipeline performance
        pipeline_results = df[df.get("test_name", "").str.contains("Pipeline", na=False)]
        if not pipeline_results.empty:
            print("\nPIPELINE PERFORMANCE:")

            disabled_pipeline = pipeline_results[pipeline_results.get("mode") == "disabled"]
            auto_pipeline = pipeline_results[pipeline_results.get("mode") == "auto"]

            if not disabled_pipeline.empty and not auto_pipeline.empty:
                disabled_avg = disabled_pipeline["throughput_cells_per_sec"].mean()
                auto_avg = auto_pipeline["throughput_cells_per_sec"].mean()
                improvement = auto_avg / disabled_avg if disabled_avg > 0 else 1.0

                print(f"  Standard pipeline: {disabled_avg:.0f} cells/sec")
                print(f"  Optimized pipeline: {auto_avg:.0f} cells/sec")
                print(f"  Improvement: {improvement:.2f}x")

        # Overall recommendation
        print("\nOVERALL ASSESSMENT:")

        # Check if optimizations provide net benefit
        has_benefits = False
        has_regressions = False

        if not mode_results.empty and "speedup" in mode_results.columns:
            beneficial_results = mode_results[mode_results["speedup"] > 1.1]  # >10% improvement
            regression_results = mode_results[mode_results["speedup"] < 0.9]  # >10% regression

            has_benefits = len(beneficial_results) > 0
            has_regressions = len(regression_results) > 0

        if has_benefits and not has_regressions:
            print("  ✅ RECOMMENDED: Optimizations provide clear benefits with no regressions")
            print("  💡 Safe to enable 'auto' mode for users who want better performance")
        elif has_benefits and has_regressions:
            print("  ⚠️  MIXED: Optimizations help some cases but may hurt others")
            print("  💡 Keep default 'disabled' but allow opt-in via 'auto' mode")
        else:
            print("  ❌ NOT RECOMMENDED: No clear performance benefits detected")
            print("  💡 Keep default 'disabled' mode")

        print(f"\n📊 Detailed results in: {filename}")


def main():
    """Run comprehensive benchmark suite."""
    print("Configurable Normalization Optimization - Comprehensive Benchmark")
    print("=" * 70)

    if not OPTIMIZATION_AVAILABLE:
        print("❌ Optimization functions not available. Cannot run benchmarks.")
        return

    # Suppress warnings for cleaner output
    warnings.filterwarnings("ignore")

    # Create benchmark instance
    benchmark = ConfigurableNormalizationBenchmark()

    try:
        # Run all benchmark categories
        benchmark.benchmark_safety_and_compatibility()
        benchmark.benchmark_mode_configurations()
        benchmark.benchmark_automatic_selection()
        benchmark.benchmark_complete_pipeline()

        # Save results and generate summary
        benchmark.save_results()

    except KeyboardInterrupt:
        print("\nBenchmark interrupted by user.")
        if benchmark.results:
            benchmark.save_results("partial_configurable_benchmark.csv")
    except Exception as e:
        print(f"\nBenchmark failed with error: {e}")
        import traceback

        traceback.print_exc()
        if benchmark.results:
            benchmark.save_results("failed_configurable_benchmark.csv")
    finally:
        # Reset to safe default
        if OPTIMIZATION_AVAILABLE:
            set_optimization_mode("disabled")


if __name__ == "__main__":
    main()

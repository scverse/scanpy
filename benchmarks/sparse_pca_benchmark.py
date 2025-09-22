#!/usr/bin/env python3
"""
Benchmark script for Sparse PCA Memory & Performance Optimization

This script compares the performance of standard PCA implementations
against the optimized sparse PCA methods for large single-cell datasets.
"""

import gc
import psutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# Add scanpy to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import scanpy as sc
from scanpy.preprocessing._pca._sparse_optimized import (
    SparsePCAConfig,
    estimate_memory_usage,
    optimized_pca,
)


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


def create_test_matrix(n_obs, n_vars, density=0.1, seed=42):
    """Create a test sparse matrix with realistic single-cell characteristics."""
    np.random.seed(seed)

    # Create base count matrix with negative binomial distribution
    # This mimics real single-cell count data
    X_dense = np.random.negative_binomial(n=5, p=0.3, size=(n_obs, n_vars))

    # Add some gene-specific scaling to create variance structure
    gene_scales = np.random.lognormal(mean=0, sigma=1, size=n_vars)
    X_dense = X_dense * gene_scales[np.newaxis, :]

    # Create sparsity by zeroing values below threshold
    # This creates realistic sparsity patterns
    threshold = np.percentile(X_dense, (1 - density) * 100)
    X_dense[X_dense < threshold] = 0

    # Convert to sparse matrix
    X_sparse = sparse.csr_matrix(X_dense, dtype=np.float32)

    return X_sparse


class SparsePCABenchmark:
    """Comprehensive benchmark suite for sparse PCA optimization."""

    def __init__(self):
        self.results = []

    def benchmark_method_comparison(self):
        """Benchmark different PCA methods on various matrix sizes."""
        print("Benchmarking PCA method comparison...")

        test_configs = [
            {"n_obs": 5000, "n_vars": 2000, "density": 0.1, "name": "Medium"},
            {"n_obs": 20000, "n_vars": 5000, "density": 0.05, "name": "Large"},
            {"n_obs": 50000, "n_vars": 10000, "density": 0.03, "name": "Very Large"},
        ]

        methods = [
            {"name": "Standard PCA (Dense)", "method": "standard_dense"},
            {"name": "Standard PCA (Sparse)", "method": "standard_sparse"},
            {"name": "Optimized Sparse PCA", "method": "optimized_sparse"},
        ]

        for config in test_configs:
            print(f"\n  Testing {config['name']} matrix: {config['n_obs']}x{config['n_vars']}")

            # Create test matrix
            X_sparse = create_test_matrix(
                config['n_obs'],
                config['n_vars'],
                config['density']
            )

            actual_density = X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1])
            print(f"    Matrix: nnz={X_sparse.nnz:,}, density={actual_density:.3f}")

            # Memory estimates
            memory_est = estimate_memory_usage(
                config['n_obs'],
                config['n_vars'],
                50,  # n_comps
                actual_density
            )

            print(f"    Estimated memory: sparse={memory_est['sparse_matrix']:.2f}GB, "
                  f"dense={memory_est['dense_matrix']:.2f}GB")

            for method_config in methods:
                try:
                    result = self._benchmark_single_method(
                        X_sparse, method_config, config, n_comps=50
                    )
                    self.results.append(result)

                    print(f"    {method_config['name']}: "
                          f"{result['time']:.2f}s, "
                          f"{result['peak_memory']:.1f}MB, "
                          f"{result['throughput']:.0f} cells/s")

                except Exception as e:
                    print(f"    {method_config['name']}: FAILED ({e})")
                    self.results.append({
                        'matrix_name': config['name'],
                        'method': method_config['name'],
                        'status': 'failed',
                        'error': str(e),
                        **config
                    })

                # Force garbage collection between methods
                gc.collect()

    def _benchmark_single_method(self, X_sparse, method_config, matrix_config, n_comps=50):
        """Benchmark a single PCA method."""
        method = method_config['method']

        # Record initial memory
        initial_memory = get_memory_usage()

        # Time the computation
        start_time = time.time()

        if method == "standard_dense":
            # Standard PCA on dense matrix (memory intensive)
            X_dense = X_sparse.toarray()
            from sklearn.decomposition import PCA
            pca = PCA(n_components=n_comps)
            pca.fit(X_dense)
            components = pca.components_
            explained_var = pca.explained_variance_

        elif method == "standard_sparse":
            # Standard sparse PCA using TruncatedSVD
            from sklearn.decomposition import TruncatedSVD
            svd = TruncatedSVD(n_components=n_comps, random_state=42)
            svd.fit(X_sparse)
            components = svd.components_
            explained_var = svd.explained_variance_

        elif method == "optimized_sparse":
            # Our optimized sparse PCA
            config = SparsePCAConfig(
                memory_limit_gb=16.0,
                use_randomized_svd=True
            )
            components, explained_var, _ = optimized_pca(
                X_sparse, n_comps, zero_center=False, config=config
            )

        end_time = time.time()

        # Record peak memory
        peak_memory = get_memory_usage()

        # Calculate metrics
        computation_time = end_time - start_time
        memory_used = peak_memory - initial_memory
        throughput = X_sparse.shape[0] / computation_time

        return {
            'matrix_name': matrix_config['name'],
            'method': method_config['name'],
            'n_obs': X_sparse.shape[0],
            'n_vars': X_sparse.shape[1],
            'nnz': X_sparse.nnz,
            'density': X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1]),
            'n_comps': n_comps,
            'time': computation_time,
            'peak_memory': peak_memory,
            'memory_used': memory_used,
            'throughput': throughput,
            'components_shape': components.shape,
            'explained_var_sum': np.sum(explained_var),
            'status': 'success'
        }

    def benchmark_memory_scaling(self):
        """Benchmark memory usage scaling with matrix size."""
        print("\nBenchmarking memory scaling...")

        base_size = 1000
        scaling_factors = [1, 2, 4, 8]

        for factor in scaling_factors:
            n_obs = base_size * factor
            n_vars = base_size * factor

            print(f"  Testing {n_obs}x{n_vars} matrix...")

            try:
                X_sparse = create_test_matrix(n_obs, n_vars, density=0.1)

                # Test memory usage
                initial_memory = get_memory_usage()

                config = SparsePCAConfig(memory_limit_gb=32.0)
                components, _, _ = optimized_pca(
                    X_sparse, 50, zero_center=False, config=config
                )

                peak_memory = get_memory_usage()
                memory_used = peak_memory - initial_memory

                # Theoretical memory estimate
                actual_density = X_sparse.nnz / (X_sparse.shape[0] * X_sparse.shape[1])
                memory_est = estimate_memory_usage(n_obs, n_vars, 50, actual_density)

                result = {
                    'scaling_factor': factor,
                    'n_obs': n_obs,
                    'n_vars': n_vars,
                    'actual_memory_mb': memory_used,
                    'estimated_memory_gb': memory_est['truncated_svd'],
                    'efficiency_ratio': memory_used / (memory_est['truncated_svd'] * 1024)
                }

                print(f"    Actual: {memory_used:.1f}MB, "
                      f"Estimated: {memory_est['truncated_svd']:.2f}GB, "
                      f"Efficiency: {result['efficiency_ratio']:.2f}")

                self.results.append(result)

            except Exception as e:
                print(f"    FAILED: {e}")

            gc.collect()

    def benchmark_optimization_methods(self):
        """Benchmark different optimization strategies."""
        print("\nBenchmarking optimization strategies...")

        # Create a large test matrix
        X_sparse = create_test_matrix(20000, 5000, density=0.05)
        print(f"  Test matrix: {X_sparse.shape[0]}x{X_sparse.shape[1]}, "
              f"nnz={X_sparse.nnz:,}")

        optimization_configs = [
            {
                "name": "Default",
                "config": SparsePCAConfig()
            },
            {
                "name": "Randomized SVD",
                "config": SparsePCAConfig(
                    use_randomized_svd=True,
                    n_oversamples=20,
                    n_iter=4
                )
            },
            {
                "name": "Memory Constrained",
                "config": SparsePCAConfig(
                    memory_limit_gb=4.0,
                    chunk_size=5000
                )
            },
            {
                "name": "High Performance",
                "config": SparsePCAConfig(
                    memory_limit_gb=32.0,
                    use_randomized_svd=True,
                    n_oversamples=15
                )
            }
        ]

        for opt_config in optimization_configs:
            print(f"  Testing {opt_config['name']} configuration...")

            try:
                start_time = time.time()
                initial_memory = get_memory_usage()

                components, explained_var, _, info = optimized_pca(
                    X_sparse, 50, zero_center=False,
                    config=opt_config['config'], return_info=True
                )

                end_time = time.time()
                peak_memory = get_memory_usage()

                result = {
                    'optimization': opt_config['name'],
                    'method_used': info['method'],
                    'time': end_time - start_time,
                    'memory_used': peak_memory - initial_memory,
                    'throughput': X_sparse.shape[0] / (end_time - start_time),
                    'explained_var_sum': np.sum(explained_var)
                }

                print(f"    Method: {info['method']}, "
                      f"Time: {result['time']:.2f}s, "
                      f"Memory: {result['memory_used']:.1f}MB")

                self.results.append(result)

            except Exception as e:
                print(f"    FAILED: {e}")

            gc.collect()

    def benchmark_centering_vs_no_centering(self):
        """Compare performance with and without centering."""
        print("\nBenchmarking centering vs no centering...")

        test_sizes = [
            (5000, 2000, "Medium"),
            (15000, 5000, "Large")
        ]

        for n_obs, n_vars, size_name in test_sizes:
            print(f"  Testing {size_name} matrix: {n_obs}x{n_vars}")

            X_sparse = create_test_matrix(n_obs, n_vars, density=0.08)

            for zero_center in [False, True]:
                center_name = "With Centering" if zero_center else "No Centering"
                print(f"    {center_name}...")

                try:
                    start_time = time.time()
                    initial_memory = get_memory_usage()

                    components, explained_var, _, info = optimized_pca(
                        X_sparse, 50, zero_center=zero_center, return_info=True
                    )

                    end_time = time.time()
                    peak_memory = get_memory_usage()

                    result = {
                        'matrix_size': size_name,
                        'n_obs': n_obs,
                        'n_vars': n_vars,
                        'zero_center': zero_center,
                        'method_used': info['method'],
                        'time': end_time - start_time,
                        'memory_used': peak_memory - initial_memory,
                        'throughput': n_obs / (end_time - start_time)
                    }

                    print(f"      Method: {info['method']}, "
                          f"Time: {result['time']:.2f}s, "
                          f"Throughput: {result['throughput']:.0f} cells/s")

                    self.results.append(result)

                except Exception as e:
                    print(f"      FAILED: {e}")

                gc.collect()

    def save_results(self, filename="sparse_pca_benchmark_results.csv"):
        """Save benchmark results to CSV file."""
        if self.results:
            df = pd.DataFrame(self.results)
            df.to_csv(filename, index=False)
            print(f"\n📊 Results saved to: {filename}")

            # Print summary statistics
            if 'time' in df.columns:
                print(f"📈 Performance Summary:")
                print(f"   Average computation time: {df['time'].mean():.2f}s")
                print(f"   Fastest method: {df.loc[df['time'].idxmin(), 'method']} "
                      f"({df['time'].min():.2f}s)")
                if 'throughput' in df.columns:
                    print(f"   Max throughput: {df['throughput'].max():.0f} cells/s")

    def print_summary(self):
        """Print a summary of benchmark results."""
        if not self.results:
            print("No results to summarize")
            return

        print("\n" + "=" * 80)
        print("SPARSE PCA OPTIMIZATION BENCHMARK SUMMARY")
        print("=" * 80)

        # Group results by benchmark type
        method_results = [r for r in self.results if 'method' in r and r.get('status') == 'success']

        if method_results:
            print(f"\n🚀 METHOD COMPARISON RESULTS:")
            df_methods = pd.DataFrame(method_results)

            for matrix_name in df_methods['matrix_name'].unique():
                matrix_results = df_methods[df_methods['matrix_name'] == matrix_name]
                print(f"\n  {matrix_name} Matrix:")

                for _, row in matrix_results.iterrows():
                    speedup = "N/A"
                    if matrix_name in df_methods['matrix_name'].values:
                        baseline = df_methods[
                            (df_methods['matrix_name'] == matrix_name) &
                            (df_methods['method'] == 'Standard PCA (Dense)')
                        ]
                        if not baseline.empty:
                            speedup = f"{baseline.iloc[0]['time'] / row['time']:.1f}x"

                    print(f"    {row['method']}: {row['time']:.2f}s "
                          f"({row['throughput']:.0f} cells/s, {speedup} speedup)")

        print(f"\n✅ Total benchmarks completed: {len([r for r in self.results if r.get('status') == 'success'])}")
        print(f"❌ Failed benchmarks: {len([r for r in self.results if r.get('status') == 'failed'])}")


def main():
    """Run the complete benchmark suite."""
    print("🧬 Sparse PCA Memory & Performance Optimization - Comprehensive Benchmark")
    print("=" * 80)

    benchmark = SparsePCABenchmark()

    try:
        # Run all benchmark suites
        benchmark.benchmark_method_comparison()
        benchmark.benchmark_memory_scaling()
        benchmark.benchmark_optimization_methods()
        benchmark.benchmark_centering_vs_no_centering()

        # Save and summarize results
        benchmark.save_results()
        benchmark.print_summary()

        print("\n🎉 Benchmark completed successfully!")

    except KeyboardInterrupt:
        print("\n⏹️  Benchmark interrupted by user")
        benchmark.save_results("sparse_pca_benchmark_partial.csv")
    except Exception as e:
        print(f"\n❌ Benchmark failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

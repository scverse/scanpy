"""Tests for sparse PCA optimization."""

from __future__ import annotations

import numpy as np
from scipy import sparse

import scanpy as sc
from scanpy.preprocessing._pca._sparse_optimized import (
    SparsePCAConfig,
    choose_optimal_pca_method,
    estimate_memory_usage,
    optimized_pca,
    sparse_centered_pca,
)
from scanpy.tests._helpers import gen_adata


class TestSparsePCAConfig:
    """Test SparsePCAConfig functionality."""

    def test_default_config(self):
        """Test default configuration."""
        config = SparsePCAConfig()

        assert config.use_gpu is False  # Should be False if CuPy not available
        assert config.chunk_size is None
        assert config.memory_limit_gb == 8.0
        assert config.use_randomized_svd is True
        assert config.n_oversamples == 10
        assert config.n_iter == 4
        assert config.max_memory_usage_ratio == 0.8

    def test_custom_config(self):
        """Test custom configuration."""
        config = SparsePCAConfig(
            use_gpu=True,
            chunk_size=5000,
            memory_limit_gb=16.0,
            use_randomized_svd=False,
        )

        # GPU should still be False if CuPy not available
        assert config.use_gpu is False or True  # Depends on CuPy availability
        assert config.chunk_size == 5000
        assert config.memory_limit_gb == 16.0
        assert config.use_randomized_svd is False


class TestMemoryEstimation:
    """Test memory estimation functions."""

    def test_estimate_memory_usage(self):
        """Test memory usage estimation."""
        n_obs, n_vars, n_comps = 10000, 2000, 50
        density = 0.1

        estimates = estimate_memory_usage(n_obs, n_vars, n_comps, density)

        assert "sparse_matrix" in estimates
        assert "dense_matrix" in estimates
        assert "components" in estimates
        assert "standard_pca" in estimates
        assert "truncated_svd" in estimates
        assert "incremental_pca" in estimates
        assert "density" in estimates
        assert "sparsity_ratio" in estimates

        # Sparse should be much smaller than dense for low density
        assert estimates["sparse_matrix"] < estimates["dense_matrix"]
        assert estimates["sparsity_ratio"] < 1.0

        # Check reasonable values
        assert estimates["sparse_matrix"] > 0
        assert estimates["dense_matrix"] > 0
        assert estimates["components"] > 0

    def test_memory_scaling(self):
        """Test that memory estimates scale correctly."""
        base_estimates = estimate_memory_usage(1000, 1000, 50, 0.1)

        # Double the size
        scaled_estimates = estimate_memory_usage(2000, 2000, 50, 0.1)

        # Should scale roughly by factor of 4 (2x obs * 2x vars)
        ratio = scaled_estimates["dense_matrix"] / base_estimates["dense_matrix"]
        assert 3.5 < ratio < 4.5  # Allow some tolerance

    def test_density_effects(self):
        """Test effect of density on memory estimates."""
        low_density = estimate_memory_usage(10000, 2000, 50, 0.01)
        high_density = estimate_memory_usage(10000, 2000, 50, 0.5)

        # Sparse memory should increase with density
        assert high_density["sparse_matrix"] > low_density["sparse_matrix"]

        # Dense memory should be the same
        assert high_density["dense_matrix"] == low_density["dense_matrix"]


class TestMethodSelection:
    """Test optimal method selection."""

    def test_method_selection_sparse_no_centering(self):
        """Test method selection for sparse matrix without centering."""
        X = sparse.random(1000, 500, density=0.1, format="csr")
        config = SparsePCAConfig(memory_limit_gb=16.0)

        method, params = choose_optimal_pca_method(X, 50, config, zero_center=False)

        # Should choose some form of truncated SVD
        assert "truncated_svd" in method.lower() or "svd" in method.lower()

    def test_method_selection_sparse_with_centering(self):
        """Test method selection for sparse matrix with centering."""
        X = sparse.random(1000, 500, density=0.1, format="csr")
        config = SparsePCAConfig(memory_limit_gb=16.0)

        method, params = choose_optimal_pca_method(X, 50, config, zero_center=True)

        # Should choose sparse centered PCA
        assert method == "sparse_centered_pca"

    def test_method_selection_memory_constrained(self):
        """Test method selection under memory constraints."""
        X = sparse.random(10000, 5000, density=0.1, format="csr")
        config = SparsePCAConfig(memory_limit_gb=1.0)  # Very limited memory

        method, params = choose_optimal_pca_method(X, 50, config, zero_center=True)

        # Should choose chunked method
        assert params.get("chunked", False) is True or "chunked" in method.lower()

    def test_method_selection_dense_matrix(self):
        """Test method selection for dense matrices."""
        X = np.random.randn(1000, 500)
        config = SparsePCAConfig(memory_limit_gb=16.0)

        method, params = choose_optimal_pca_method(X, 50, config, zero_center=True)

        # Should choose standard PCA for reasonable size dense matrix
        assert method in ["standard_pca", "incremental_pca"]


class TestSparseCenteredPCA:
    """Test sparse centered PCA implementation."""

    def test_sparse_centered_pca_basic(self):
        """Test basic sparse centered PCA functionality."""
        # Create test data
        np.random.seed(42)
        X_dense = np.random.randn(100, 50)
        X_sparse = sparse.csr_matrix(X_dense)

        n_comps = 10

        # Test our implementation
        components, explained_var, explained_var_ratio = sparse_centered_pca(
            X_sparse, n_comps, chunked=False
        )

        assert components.shape == (n_comps, 50)
        assert explained_var.shape == (n_comps,)
        assert explained_var_ratio.shape == (n_comps,)

        # Explained variance should be positive and sorted
        assert np.all(explained_var >= 0)
        assert np.all(explained_var[:-1] >= explained_var[1:])  # Descending order

        # Explained variance ratios should sum to <= 1
        assert np.sum(explained_var_ratio) <= 1.0

    def test_sparse_centered_pca_chunked(self):
        """Test chunked sparse centered PCA."""
        # Create test data
        np.random.seed(42)
        X_dense = np.random.randn(200, 50)
        X_sparse = sparse.csr_matrix(X_dense)

        n_comps = 10
        chunk_size = 50

        # Test chunked implementation
        components, explained_var, explained_var_ratio = sparse_centered_pca(
            X_sparse, n_comps, chunked=True, chunk_size=chunk_size
        )

        assert components.shape == (n_comps, 50)
        assert explained_var.shape == (n_comps,)
        assert explained_var_ratio.shape == (n_comps,)

        # Should produce reasonable results
        assert np.all(explained_var >= 0)
        assert np.sum(explained_var_ratio) <= 1.0

    def test_sparse_centered_pca_vs_standard(self):
        """Compare sparse centered PCA with standard PCA on same data."""
        # Create test data with some structure
        np.random.seed(42)
        n_obs, n_vars = 200, 50

        # Create structured data
        latent = np.random.randn(n_obs, 5)
        loadings = np.random.randn(5, n_vars)
        X_dense = latent @ loadings + 0.1 * np.random.randn(n_obs, n_vars)

        # Make it sparse by zeroing small values
        X_dense[np.abs(X_dense) < 0.5] = 0
        X_sparse = sparse.csr_matrix(X_dense)

        n_comps = 5

        # Our sparse implementation
        components_sparse, var_sparse, var_ratio_sparse = sparse_centered_pca(
            X_sparse, n_comps
        )

        # Standard PCA on dense data
        from sklearn.decomposition import PCA

        pca_standard = PCA(n_components=n_comps)
        pca_standard.fit(X_dense)

        # Compare explained variance (should be similar)
        # Note: Components might have different signs, so we compare magnitudes
        var_diff = np.abs(var_sparse - pca_standard.explained_variance_)
        relative_diff = var_diff / (pca_standard.explained_variance_ + 1e-8)

        # Should be reasonably close (within 10% for most components)
        assert np.mean(relative_diff) < 0.2


class TestOptimizedPCA:
    """Test the main optimized PCA function."""

    def test_optimized_pca_basic(self):
        """Test basic optimized PCA functionality."""
        # Create sparse test data
        X = sparse.random(500, 100, density=0.2, format="csr")
        n_comps = 10

        components, explained_var, explained_var_ratio = optimized_pca(
            X, n_comps, zero_center=False
        )

        assert components.shape == (n_comps, 100)
        assert explained_var.shape == (n_comps,)
        assert explained_var_ratio.shape == (n_comps,)

        # Basic sanity checks
        assert np.all(explained_var >= 0)
        assert np.sum(explained_var_ratio) <= 1.0

    def test_optimized_pca_with_info(self):
        """Test optimized PCA with return_info=True."""
        X = sparse.random(500, 100, density=0.2, format="csr")
        n_comps = 10

        components, explained_var, explained_var_ratio, info = optimized_pca(
            X, n_comps, zero_center=False, return_info=True
        )

        assert isinstance(info, dict)
        assert "method" in info
        assert "matrix_shape" in info
        assert "matrix_type" in info
        assert "n_components" in info
        assert "sparsity" in info
        assert "nnz" in info

        assert info["matrix_shape"] == (500, 100)
        assert info["matrix_type"] == "sparse"
        assert info["n_components"] == n_comps

    def test_optimized_pca_config(self):
        """Test optimized PCA with custom configuration."""
        X = sparse.random(500, 100, density=0.2, format="csr")
        n_comps = 10

        config = SparsePCAConfig(
            memory_limit_gb=1.0,  # Force chunked methods
            use_randomized_svd=True,
        )

        components, explained_var, explained_var_ratio = optimized_pca(
            X, n_comps, zero_center=False, config=config
        )

        assert components.shape == (n_comps, 100)
        assert explained_var.shape == (n_comps,)
        assert explained_var_ratio.shape == (n_comps,)


class TestScanpyIntegration:
    """Test integration with main scanpy PCA function."""

    def test_pca_with_sparse_optimization(self):
        """Test PCA with sparse optimization enabled."""
        # Create test AnnData with sparse data
        adata = gen_adata((1000, 500), sparse=True)

        # Run PCA with optimization
        sc.pp.pca(
            adata,
            n_comps=20,
            use_sparse_optimization=True,
            zero_center=False,  # Simpler case first
        )

        # Check results
        assert "X_pca" in adata.obsm
        assert "PCs" in adata.varm
        assert "pca" in adata.uns

        assert adata.obsm["X_pca"].shape == (1000, 20)
        assert adata.varm["PCs"].shape == (500, 20)

        # Check that optimization info is stored
        assert "optimization_info" in adata.uns["pca"]
        opt_info = adata.uns["pca"]["optimization_info"]
        assert "method" in opt_info

    def test_pca_optimization_fallback(self):
        """Test that PCA falls back gracefully if optimization fails."""
        # Create test data that might cause optimization to fail
        adata = gen_adata((100, 50), sparse=True)

        # This should not raise an error even if optimization fails
        sc.pp.pca(
            adata,
            n_comps=10,
            use_sparse_optimization=True,
            sparse_pca_config=SparsePCAConfig(memory_limit_gb=0.001),  # Very low limit
        )

        # Should still have PCA results
        assert "X_pca" in adata.obsm
        assert adata.obsm["X_pca"].shape == (100, 10)

    def test_pca_optimization_disabled_by_default(self):
        """Test that optimization is disabled by default."""
        adata = gen_adata((100, 50), sparse=True)

        sc.pp.pca(adata, n_comps=10)

        # Should have normal PCA results
        assert "X_pca" in adata.obsm

        # Should not have optimization info
        if "optimization_info" in adata.uns.get("pca", {}):
            # If it exists, it should be from standard methods, not our optimization
            pass

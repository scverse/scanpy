"""Optimized normalization functions with configurable implementation selection.

This module provides enhanced normalization functions with improved performance
and numerical stability. It includes multiple implementation strategies that can
be controlled via configuration flags to ensure compatibility and optimal performance.

Implementation Modes:
- 'disabled': Use only standard scanpy functions (no optimizations)
- 'naive': Use simple optimized functions without Numba
- 'numba': Use Numba-optimized functions (best for large matrices)
- 'auto': Automatically select best implementation based on matrix size (default)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numba
import numpy as np

from .._compat import CSRBase


if TYPE_CHECKING:
    from numpy.typing import NDArray

# Type for implementation selection
ImplementationMode = Literal["disabled", "naive", "numba", "auto"]

# Global configuration for optimization mode
_OPTIMIZATION_MODE: ImplementationMode = "auto"  # Intelligent default - automatic selection

# Performance thresholds for automatic implementation selection
_NUMBA_THRESHOLD_ELEMENTS = 100_000  # Use Numba for matrices with >100k elements
_NUMBA_THRESHOLD_NNZ = 50_000  # Use Numba for matrices with >50k non-zeros


def set_optimization_mode(mode: ImplementationMode) -> None:
    """Set the global optimization mode for normalization functions.

    Parameters
    ----------
    mode : {'disabled', 'naive', 'numba', 'auto'}
        Optimization mode to use:
        - 'disabled': Use standard scanpy functions only (default, safe)
        - 'naive': Use simple optimized functions without Numba
        - 'numba': Use Numba-optimized functions (best for large matrices)
        - 'auto': Automatically select best implementation based on matrix size

    Examples
    --------
    >>> import scanpy as sc
    >>> # Enable automatic optimization (recommended)
    >>> sc.pp.set_optimization_mode("auto")
    >>>
    >>> # Force Numba optimization for all matrices
    >>> sc.pp.set_optimization_mode("numba")
    >>>
    >>> # Disable optimizations (safe default)
    >>> sc.pp.set_optimization_mode("disabled")
    """
    global _OPTIMIZATION_MODE
    if mode not in ["disabled", "naive", "numba", "auto"]:
        raise ValueError(
            f"Invalid optimization mode: {mode}. Must be one of: 'disabled', 'naive', 'numba', 'auto'"
        )
    _OPTIMIZATION_MODE = mode


def get_optimization_mode() -> ImplementationMode:
    """Get the current optimization mode.

    Returns
    -------
    str
        Current optimization mode
    """
    return _OPTIMIZATION_MODE


def _should_use_numba_optimization(mat: CSRBase) -> bool:
    """Determine whether to use Numba optimization based on matrix characteristics.

    Numba has compilation overhead that makes it slower for small matrices,
    but provides significant speedups for larger matrices. This function
    automatically selects the best implementation.

    Parameters
    ----------
    mat : CSRBase
        Sparse CSR matrix to analyze

    Returns
    -------
    bool
        True if Numba optimization should be used, False for naive implementation
    """
    total_elements = mat.shape[0] * mat.shape[1]
    nnz = mat.nnz

    # Use Numba if matrix is large enough to overcome compilation overhead
    return total_elements >= _NUMBA_THRESHOLD_ELEMENTS or nnz >= _NUMBA_THRESHOLD_NNZ


@numba.njit(cache=True, parallel=False)
def _csr_sum_and_squared_sum_optimized(
    data: NDArray[np.floating],
    indices: NDArray[np.integer],
    indptr: NDArray[np.integer],
    n_rows: int,
    n_cols: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Numba-accelerated sum and squared sum computation for CSR matrix.

    This optimized version provides significant speedup over the standard
    implementation by using efficient memory access patterns and compiled code.
    Note: Currently uses sequential processing to avoid race conditions with
    parallel writes to the same array elements.

    Parameters
    ----------
    data : NDArray[np.floating]
        CSR matrix data array
    indices : NDArray[np.integer]
        CSR matrix indices array
    indptr : NDArray[np.integer]
        CSR matrix index pointer array
    n_rows : int
        Number of rows in the matrix
    n_cols : int
        Number of columns in the matrix

    Returns
    -------
    tuple[NDArray[np.float64], NDArray[np.float64]]
        Column sums and squared sums
    """
    sums = np.zeros(n_cols, dtype=np.float64)
    sq_sums = np.zeros(n_cols, dtype=np.float64)

    for row in range(n_rows):  # Sequential to avoid race conditions
        start = indptr[row]
        end = indptr[row + 1]
        for idx in range(start, end):
            col = indices[idx]
            val = data[idx]
            sums[col] += val
            sq_sums[col] += val * val

    return sums, sq_sums


@numba.njit(cache=True, parallel=True)
def _csr_row_scaling_optimized(
    data: NDArray[np.floating],
    indptr: NDArray[np.integer],
    scaling_factors: NDArray[np.floating],
    n_rows: int,
) -> None:
    """Numba-accelerated in-place row scaling for CSR matrix.

    This function efficiently scales each row of a CSR matrix by the
    corresponding scaling factor, using optimized memory access patterns.

    Parameters
    ----------
    data : NDArray[np.floating]
        CSR matrix data array (modified in-place)
    indptr : NDArray[np.integer]
        CSR matrix index pointer array
    scaling_factors : NDArray[np.floating]
        Scaling factor for each row
    n_rows : int
        Number of rows in the matrix
    """
    for row in numba.prange(n_rows):
        start_idx = indptr[row]
        end_idx = indptr[row + 1]
        scale = scaling_factors[row]
        for idx in range(start_idx, end_idx):
            data[idx] *= scale


@numba.njit(cache=True, parallel=True)
def _compute_row_sums_optimized(
    data: NDArray[np.floating], indptr: NDArray[np.integer], n_rows: int
) -> NDArray[np.float64]:
    """Optimized computation of row sums for CSR matrix.

    Parameters
    ----------
    data : NDArray[np.floating]
        CSR matrix data array
    indptr : NDArray[np.integer]
        CSR matrix index pointer array
    n_rows : int
        Number of rows in the matrix

    Returns
    -------
    NDArray[np.float64]
        Sum for each row
    """
    row_sums = np.zeros(n_rows, dtype=np.float64)

    for row in numba.prange(n_rows):
        start = indptr[row]
        end = indptr[row + 1]
        row_sum = 0.0
        for idx in range(start, end):
            row_sum += data[idx]
        row_sums[row] = row_sum

    return row_sums


@numba.njit(cache=True, parallel=False)
def _count_highly_expressed_optimized(
    data: NDArray[np.floating],
    indices: NDArray[np.integer],
    indptr: NDArray[np.integer],
    counts_per_cell: NDArray[np.floating],
    max_fraction: float,
    n_rows: int,
    n_cols: int,
) -> NDArray[np.int32]:
    """Optimized counting of highly expressed genes per cell.

    Parameters
    ----------
    data : NDArray[np.floating]
        CSR matrix data array
    indices : NDArray[np.integer]
        CSR matrix indices array
    indptr : NDArray[np.integer]
        CSR matrix index pointer array
    counts_per_cell : NDArray[np.floating]
        Total counts per cell
    max_fraction : float
        Maximum fraction threshold for highly expressed genes
    n_rows : int
        Number of rows (cells)
    n_cols : int
        Number of columns (genes)

    Returns
    -------
    NDArray[np.int32]
        Count of cells where each gene is highly expressed
    """
    counts_per_cols = np.zeros(n_cols, dtype=np.int32)

    for row in range(n_rows):  # Sequential to avoid race conditions
        threshold = max_fraction * counts_per_cell[row]
        start = indptr[row]
        end = indptr[row + 1]

        for idx in range(start, end):
            if data[idx] > threshold:
                col = indices[idx]
                counts_per_cols[col] += 1

    return counts_per_cols


def _normalize_csr_naive(
    mat: CSRBase,
    *,
    rows: int,
    columns: int,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
) -> tuple[NDArray[np.floating], NDArray[np.int32] | None]:
    """Naive CSR matrix normalization using standard NumPy operations.

    This implementation is faster for small matrices where Numba overhead
    dominates performance.

    Parameters
    ----------
    mat : CSRBase
        Sparse CSR matrix to normalize
    rows : int
        Number of rows in the matrix
    columns : int
        Number of columns in the matrix
    exclude_highly_expressed : bool, default False
        Whether to exclude highly expressed genes from normalization
    max_fraction : float, default 0.05
        Maximum fraction for highly expressed gene threshold
    n_threads : int, default 10
        Number of threads for parallel processing (unused in naive implementation)

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.int32] | None]
        Counts per cell and optionally counts per column for highly expressed genes
    """
    # Simple row sum computation using scipy
    counts_per_cell = np.array(mat.sum(axis=1)).flatten().astype(np.float64)

    counts_per_cols = None
    if exclude_highly_expressed:
        # Find highly expressed genes
        counts_per_cols = np.zeros(columns, dtype=np.int32)

        # Iterate through matrix to count highly expressed genes
        for i in range(rows):
            threshold = max_fraction * counts_per_cell[i]
            start = mat.indptr[i]
            end = mat.indptr[i + 1]

            for j in range(start, end):
                if mat.data[j] > threshold:
                    gene_idx = mat.indices[j]
                    counts_per_cols[gene_idx] += 1

        # Recompute counts excluding highly expressed genes
        counts_per_cell_filtered = np.zeros(rows, dtype=np.float64)

        for i in range(rows):
            count = 0.0
            start = mat.indptr[i]
            end = mat.indptr[i + 1]

            for j in range(start, end):
                gene_idx = mat.indices[j]
                if counts_per_cols[gene_idx] == 0:  # Not highly expressed
                    count += mat.data[j]
            counts_per_cell_filtered[i] = count

        counts_per_cell = counts_per_cell_filtered

    return counts_per_cell, counts_per_cols


def _normalize_csr_numba(
    mat: CSRBase,
    *,
    rows: int,
    columns: int,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
) -> tuple[NDArray[np.floating], NDArray[np.int32] | None]:
    """Numba-optimized CSR matrix normalization for large matrices.

    This function provides significant performance improvements for large
    matrices through Numba compilation and optimized algorithms.

    Parameters
    ----------
    mat : CSRBase
        Sparse CSR matrix to normalize
    rows : int
        Number of rows in the matrix
    columns : int
        Number of columns in the matrix
    exclude_highly_expressed : bool, default False
        Whether to exclude highly expressed genes from normalization
    max_fraction : float, default 0.05
        Maximum fraction for highly expressed gene threshold
    n_threads : int, default 10
        Number of threads for parallel processing (currently unused)

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.int32] | None]
        Counts per cell and optionally counts per column for highly expressed genes
    """
    # Compute row sums efficiently using Numba
    counts_per_cell = _compute_row_sums_optimized(mat.data, mat.indptr, rows)

    counts_per_cols = None
    if exclude_highly_expressed:
        # Count highly expressed genes using Numba
        counts_per_cols = _count_highly_expressed_optimized(
            mat.data,
            mat.indices,
            mat.indptr,
            counts_per_cell,
            max_fraction,
            rows,
            columns,
        )

        # Recompute row sums excluding highly expressed genes
        counts_per_cell_filtered = np.zeros(rows, dtype=np.float64)

        for i in range(rows):
            count = 0.0
            for j in range(mat.indptr[i], mat.indptr[i + 1]):
                gene_idx = mat.indices[j]
                if counts_per_cols[gene_idx] == 0:  # Not highly expressed
                    count += mat.data[j]
            counts_per_cell_filtered[i] = count

        counts_per_cell = counts_per_cell_filtered

    return counts_per_cell, counts_per_cols


def _normalize_csr_optimized(
    mat: CSRBase,
    *,
    rows: int,
    columns: int,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
    mode: ImplementationMode | None = None,
) -> tuple[NDArray[np.floating], NDArray[np.int32] | None]:
    """CSR matrix normalization with configurable implementation selection.

    This function chooses between different implementations based on the
    global optimization mode or an explicit mode parameter. This ensures
    compatibility and optimal performance across all use cases.

    Parameters
    ----------
    mat : CSRBase
        Sparse CSR matrix to normalize
    rows : int
        Number of rows in the matrix
    columns : int
        Number of columns in the matrix
    exclude_highly_expressed : bool, default False
        Whether to exclude highly expressed genes from normalization
    max_fraction : float, default 0.05
        Maximum fraction for highly expressed gene threshold
    n_threads : int, default 10
        Number of threads for parallel processing
    mode : {'disabled', 'naive', 'numba', 'auto'} or None
        Override global optimization mode for this call

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.int32] | None]
        Counts per cell and optionally counts per column for highly expressed genes

    Notes
    -----
    Implementation selection modes:
    - 'disabled': Returns None (caller should use standard scanpy functions)
    - 'naive': Uses simple optimized functions without Numba
    - 'numba': Uses Numba-optimized functions
    - 'auto': Automatically selects based on matrix size

    Automatic selection criteria:
    - Total matrix elements (>100k triggers Numba)
    - Number of non-zero elements (>50k triggers Numba)

    Examples
    --------
    >>> import scipy.sparse as sp
    >>> import numpy as np
    >>> # Use global mode (default: disabled)
    >>> X = sp.random(1000, 500, density=0.1, format="csr")
    >>> result = _normalize_csr_optimized(X, rows=1000, columns=500)
    >>> if result is None:
    ...     # Use standard scanpy function
    ...     pass
    >>>
    >>> # Override mode for this call
    >>> result = _normalize_csr_optimized(X, rows=1000, columns=500, mode="auto")
    """
    # Use provided mode or global mode
    current_mode = mode if mode is not None else _OPTIMIZATION_MODE

    # Handle disabled mode - return None to signal caller should use standard functions
    if current_mode == "disabled":
        return None

    # Select implementation based on mode
    if current_mode == "naive":
        return _normalize_csr_naive(
            mat,
            rows=rows,
            columns=columns,
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            n_threads=n_threads,
        )
    elif current_mode == "numba":
        return _normalize_csr_numba(
            mat,
            rows=rows,
            columns=columns,
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            n_threads=n_threads,
        )
    elif current_mode == "auto":
        # Automatic selection based on matrix size
        if _should_use_numba_optimization(mat):
            return _normalize_csr_numba(
                mat,
                rows=rows,
                columns=columns,
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                n_threads=n_threads,
            )
        else:
            return _normalize_csr_naive(
                mat,
                rows=rows,
                columns=columns,
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                n_threads=n_threads,
            )
    else:
        raise ValueError(f"Invalid optimization mode: {current_mode}")


def apply_row_normalization_optimized(
    mat: CSRBase,
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    mode: ImplementationMode | None = None,
) -> tuple[NDArray[np.floating], NDArray[np.floating]] | None:
    """Apply optimized row normalization to CSR matrix with configurable implementation selection.

    This function normalizes each row (cell) to have the specified target sum,
    with optional exclusion of highly expressed genes from the normalization
    calculation. Implementation is controlled by the global optimization mode.

    Parameters
    ----------
    mat : CSRBase
        Sparse CSR matrix to normalize (modified in-place)
    target_sum : float, default 1e4
        Target sum for each row after normalization
    exclude_highly_expressed : bool, default False
        Whether to exclude highly expressed genes from normalization
    max_fraction : float, default 0.05
        Maximum fraction for highly expressed gene threshold
    mode : {'disabled', 'naive', 'numba', 'auto'} or None
        Override global optimization mode for this call

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.floating]] or None
        Normalization factors and counts per cell, or None if disabled

    Notes
    -----
    Implementation selection modes:
    - 'disabled': Returns None (caller should use standard scanpy functions)
    - 'naive': Uses simple optimized functions without Numba
    - 'numba': Uses Numba-optimized functions
    - 'auto': Automatically selects based on matrix size

    For 'auto' mode:
    - Small matrices (<100k elements or <50k non-zeros): Uses naive NumPy
    - Large matrices: Uses Numba-optimized implementation

    This ensures optimal performance across all matrix sizes while maintaining
    identical results regardless of implementation used.

    Examples
    --------
    >>> import scipy.sparse as sp
    >>> import numpy as np
    >>> # Default mode (disabled) - returns None
    >>> X = sp.random(1000, 500, density=0.1, format="csr", dtype=np.float64)
    >>> result = apply_row_normalization_optimized(X)
    >>> if result is None:
    ...     # Use standard scanpy normalization
    ...     pass
    >>>
    >>> # Enable optimization for this call
    >>> result = apply_row_normalization_optimized(X, mode="auto")
    >>> if result is not None:
    ...     scaling_factors, counts = result
    """
    rows, columns = mat.shape

    # Get counts per cell using configurable selection
    result = _normalize_csr_optimized(
        mat,
        rows=rows,
        columns=columns,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        mode=mode,
    )

    # Handle disabled mode
    if result is None:
        return None

    counts_per_cell, _ = result

    # Compute scaling factors
    # Avoid division by zero by using a small epsilon
    eps = 1e-10
    scaling_factors = target_sum / (counts_per_cell + eps)

    # Apply scaling in-place using appropriate method
    current_mode = mode if mode is not None else _OPTIMIZATION_MODE

    if current_mode == "numba" or (
        current_mode == "auto" and _should_use_numba_optimization(mat)
    ):
        # Use Numba-optimized scaling for large matrices or forced Numba mode
        _csr_row_scaling_optimized(mat.data, mat.indptr, scaling_factors, rows)
    else:
        # Use simple scaling for small matrices or naive mode
        for i in range(rows):
            start_idx = mat.indptr[i]
            end_idx = mat.indptr[i + 1]
            scale = scaling_factors[i]
            for idx in range(start_idx, end_idx):
                mat.data[idx] *= scale

    return scaling_factors, counts_per_cell

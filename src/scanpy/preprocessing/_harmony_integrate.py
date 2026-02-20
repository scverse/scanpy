from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import DTypeLike


def harmony_integrate(  # noqa: PLR0913
    adata: AnnData,
    key: str | list[str],
    *,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    dtype: DTypeLike = np.float64,
    theta: float | list[float] | None = None,
    sigma: float = 0.1,
    n_clusters: int | None = None,
    max_iter_harmony: int = 10,
    max_iter_clustering: int = 200,
    tol_harmony: float = 1e-4,
    tol_clustering: float = 1e-5,
    ridge_lambda: float = 1.0,
    correction_method: Literal["fast", "original"] = "original",
    block_proportion: float = 0.05,
    random_state: int | None = 0,
    sparse: bool = False,
) -> None:
    """
    Integrate different experiments using the Harmony algorithm.

    This CPU implementation is based on the harmony-pytorch & rapids_singlecell version,
    using NumPy for efficient computation.

    Parameters
    ----------
    adata
        The annotated data matrix.
    key
        The key(s) of the column(s) in ``adata.obs`` that differentiates
        among experiments/batches.
    basis
        The name of the field in ``adata.obsm`` where the PCA table is
        stored. Defaults to ``'X_pca'``.
    adjusted_basis
        The name of the field in ``adata.obsm`` where the adjusted PCA
        table will be stored. Defaults to ``X_pca_harmony``.
    dtype
        The data type to use for Harmony computation.
    theta
        Diversity penalty weight(s). Default is 2 for each batch variable.
    sigma
        Width of soft clustering kernel. Default 0.1.
    n_clusters
        Number of clusters. Default is min(100, n_cells/30).
    max_iter_harmony
        Maximum Harmony iterations. Default 10.
    max_iter_clustering
        Maximum clustering iterations per Harmony round. Default 200.
    tol_harmony
        Convergence tolerance for Harmony. Default 1e-4.
    tol_clustering
        Convergence tolerance for clustering. Default 1e-5.
    ridge_lambda
        Ridge regression regularization. Default 1.0.
    correction_method
        Choose which method for the correction step: ``original`` for
        original method, ``fast`` for improved method.
    block_proportion
        Fraction of cells processed per clustering iteration. Default 0.05.
    random_state
        Random seed for reproducibility.
    sparse
        Use sparse matrices for batch encoding. Reduces memory for large datasets.

    Returns
    -------
    Updates adata with the field ``adata.obsm[adjusted_basis]``,
    containing principal components adjusted by Harmony.
    """
    from ._harmony import harmonize

    # Ensure the basis exists in adata.obsm
    if basis not in adata.obsm:
        msg = (
            f"The specified basis {basis!r} is not available in `adata.obsm`. "
            f"Available bases: {list(adata.obsm.keys())}"
        )
        raise ValueError(msg)

    # Get the input data
    input_data = adata.obsm[basis]

    # Convert to numpy array with specified dtype
    try:
        x = np.ascontiguousarray(input_data, dtype=dtype)
    except Exception as e:
        msg = (
            f"Could not convert input of type {type(input_data).__name__} "
            "to NumPy array."
        )
        raise TypeError(msg) from e

    # Check for NaN values
    if np.isnan(x).any():
        msg = (
            "Input data contains NaN values. Please handle these before "
            "running harmony_integrate."
        )
        raise ValueError(msg)

    # Run Harmony
    harmony_out = harmonize(
        x,
        adata.obs,
        key,
        theta=theta,
        sigma=sigma,
        n_clusters=n_clusters,
        max_iter_harmony=max_iter_harmony,
        max_iter_clustering=max_iter_clustering,
        tol_harmony=tol_harmony,
        tol_clustering=tol_clustering,
        ridge_lambda=ridge_lambda,
        correction_method=correction_method,
        block_proportion=block_proportion,
        random_state=random_state,
        sparse=sparse,
    )

    # Store result
    adata.obsm[adjusted_basis] = harmony_out

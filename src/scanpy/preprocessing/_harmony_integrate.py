from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


def harmony_integrate(
    adata: AnnData,
    key: str | list[str],
    *,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    dtype: type = np.float64,
    correction_method: Literal["fast", "original"] = "original",
    sparse: bool = False,
    **kwargs,
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
    correction_method
        Choose which method for the correction step: ``original`` for
        original method, ``fast`` for improved method.
    sparse
        Use sparse matrices for batch encoding. Reduces memory for large datasets.
    **kwargs
        Additional arguments passed to ``harmonize()``.

    Returns
    -------
    Updates adata with the field ``adata.obsm[adjusted_basis]``,
    containing principal components adjusted by Harmony.
    """
    from ._harmony import harmonize

    # Ensure the basis exists in adata.obsm
    if basis not in adata.obsm:
        msg = (
            f"The specified basis '{basis}' is not available in adata.obsm. "
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
        correction_method=correction_method,
        sparse=sparse,
        **kwargs,
    )

    # Store result
    adata.obsm[adjusted_basis] = harmony_out

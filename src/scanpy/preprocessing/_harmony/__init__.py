from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..._compat import warn

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import DTypeLike

    from ..._utils.random import RNGLike, SeedLike


def harmony_integrate(  # noqa: PLR0913
    adata: AnnData,
    key: str | Sequence[str],
    *,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    dtype: DTypeLike = np.float64,
    flavor: Literal["harmony2", "harmony1"] = "harmony2",
    n_clusters: int | None = None,
    max_iter_harmony: int = 10,
    max_iter_clustering: int = 200,
    tol_harmony: float = 1e-4,
    tol_clustering: float = 1e-5,
    sigma: float = 0.1,
    theta: float | Sequence[float] = 2.0,
    tau: int = 0,
    ridge_lambda: float = 1.0,
    alpha: float = 0.2,
    batch_prune_threshold: float | None = 1e-5,
    correction_method: Literal["fast"] = "fast",
    block_proportion: float = 0.05,
    rng: SeedLike | RNGLike | None = None,
) -> None:
    """Integrate different experiments using the Harmony algorithm :cite:p:`Korsunsky2019,Patikas2026`.

    This CPU implementation was originally based on the harmony-pytorch and
    rapids-singlecell implementations, using NumPy for efficient computation.
    Multiple batch variables follow the per-covariate formulation in the Harmony
    papers: each key is modeled separately instead of combining keys into one
    joint category.
    As Harmony works by adjusting the principal components,
    this function should be run after performing PCA but before computing the neighbor graph.

    By default, the Harmony2 algorithm is used,
    which includes a stabilized diversity penalty,
    dynamic per-cluster-per-batch ridge regularization,
    and automatic batch pruning.
    To revert to the original Harmony behavior::

        sc.pp.harmony_integrate(adata, key, flavor="harmony1")

    .. array-support:: pp.harmony_integrate

    Parameters
    ----------
    adata
        The annotated data matrix.
    key
        The key(s) of the column(s) in ``adata.obs`` that differentiate(s) among experiments/batches.
        Multiple keys are modeled as separate batch variables, with one active
        categorical level per variable and cell. To retain the joint-combination
        behavior of earlier releases, combine the desired columns into one
        categorical column and pass that single key.
    basis
        The name of the field in ``adata.obsm`` where the PCA table is stored.
    adjusted_basis
        The name of the field in ``adata.obsm`` where the adjusted PCA table will be stored.
    dtype
        The data type to use for Harmony computation.
        If you use 32-bit you may experience numerical instability.
    flavor
        Which version of the Harmony algorithm to use.
        ``"harmony2"`` (default) enables the stabilized diversity penalty,
        dynamic per-cluster-per-batch ridge regularization,
        and automatic batch pruning from :cite:p:`Patikas2026`.
        ``"harmony1"`` uses the original algorithm from :cite:p:`Korsunsky2019`.
    n_clusters
        Number of clusters used for soft k-means in the Harmony algorithm.
        If ``None``, uses ``min(100, N / 30)``.
        More clusters capture finer-grained structure but increase computation time.
    max_iter_harmony
        Maximum number of outer Harmony iterations
        (each consisting of a clustering step followed by a correction step).
    max_iter_clustering
        Maximum iterations for the clustering step within each Harmony iteration.
    tol_harmony
        Convergence tolerance for the Harmony objective function.
        The algorithm stops when the relative change in objective falls below this value.
    tol_clustering
        Convergence tolerance for the clustering step within each
        Harmony iteration.
    sigma
        Width of the soft-clustering kernel.
        Controls the entropy of cluster assignments:
        smaller values produce harder assignments (cells assigned to fewer clusters),
        while larger values produce softer assignments (cells spread across more clusters).
    theta
        Diversity penalty weight per batch variable.
        Controls how strongly Harmony encourages each cluster
        to contain a balanced representation of all batches.
        Higher values (e.g. ``4``) produce more aggressive mixing;
        lower values (e.g. ``0.5``) allow more batch-specific clusters.
        Set to ``0`` to disable the diversity penalty for a batch variable.
        A scalar is applied to every key. A sequence may contain one value per
        key, expanded over that key's categorical levels, or one value per
        categorical level across all keys.
    tau
        Discounting factor on ``theta``.
        When ``tau > 0``,
        the diversity penalty is down-weighted for batches with fewer cells,
        preventing over-correction of small batches.
        By default (``0``), there is no discounting.
    ridge_lambda
        Ridge regression regularization for the correction step.
        Larger values produce more conservative (smaller) corrections,
        preventing over-fitting.
        Only used with ``flavor="harmony1"``.
    alpha
        Scaling factor for the dynamic per-cluster-per-batch ridge regularization.
        The effective regularization for each cluster-batch pair is ``alpha * E_kb``
        where ``E_kb`` is the expected number of cells.
        Larger values produce more conservative corrections.
        Only used with ``flavor="harmony2"``.
    batch_prune_threshold
        Fraction threshold below which a batch-cluster pair is pruned (correction suppressed).
        When the fraction of a batch’s cells assigned to a cluster (``O_kb / N_b``) falls below this threshold,
        that batch-cluster pair receives no correction, preventing spurious adjustments.
        Only used with ``flavor="harmony2"``.
        Set to ``None`` to disable pruning.
    correction_method
        Method for the correction step.
        ``"fast"`` uses a precomputed factorization that avoids the full inversion,
        which is efficient for datasets with one batch variable. Multiple keys
        automatically use the exact general-design solve because this optimization
        only applies to a single batch variable.
    block_proportion
        Proportion of cells updated per clustering sub-iteration.
        Smaller values produce more stochastic updates.
        Larger values are faster but may converge to different solutions.
    rng
        Random number generator or seed for deterministic behavior.

    Returns
    -------
    Updates adata with the field ``adata.obsm[adjusted_basis]``,
    containing principal components adjusted by Harmony
    such that different experiments are integrated.
    """
    from .core import Harmony

    # Resolve flavor into internal flags
    if flavor not in {"harmony1", "harmony2"}:
        msg = f"flavor must be 'harmony1' or 'harmony2', got {flavor!r}."
        raise ValueError(msg)
    if correction_method != "fast":
        msg = f"correction_method must be 'fast', got {correction_method!r}."
        raise ValueError(msg)
    stabilized_penalty = flavor == "harmony2"
    dynamic_lambda = flavor == "harmony2"

    # Warn when flavor-incompatible parameters are explicitly set
    if flavor == "harmony2" and ridge_lambda != 1.0:
        warn(
            "ridge_lambda is ignored when flavor='harmony2'; "
            "use alpha to control regularization strength.",
            UserWarning,
        )
    if flavor == "harmony1":
        if alpha != 0.2:
            warn(
                "alpha is ignored when flavor='harmony1'; use ridge_lambda instead.",
                UserWarning,
            )
        if batch_prune_threshold != 1e-5:
            warn(
                "batch_prune_threshold is ignored when flavor='harmony1'.",
                UserWarning,
            )

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
        msg = f"Could not convert input of type {type(input_data).__name__} to NumPy array."
        raise TypeError(msg) from e

    # Check for NaN values
    if np.isnan(x).any():
        msg = (
            "Input data contains NaN values. Please handle these before "
            "running harmony_integrate."
        )
        raise ValueError(msg)

    # Run Harmony
    harmony = Harmony(
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
        block_proportion=block_proportion,
        tau=tau,
        rng=rng,
        stabilized_penalty=stabilized_penalty,
        dynamic_lambda=dynamic_lambda,
        alpha=alpha,
        batch_prune_threshold=batch_prune_threshold,
    )
    harmony_out = harmony.fit(x)

    # Store result
    adata.obsm[adjusted_basis] = harmony_out

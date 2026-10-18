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
    max_iter_clustering: int | None = None,
    tol_harmony: float | None = None,
    tol_clustering: float | None = None,
    sigma: float = 0.1,
    theta: float | Sequence[float] = 2.0,
    tau: int = 0,
    ridge_lambda: float = 1.0,
    alpha: float = 0.2,
    batch_prune_threshold: float | None = 1e-5,
    correction_method: Literal["fast"] = "fast",
    block_proportion: float = 0.05,
    stop_rule: Literal[
        "reference_record", "reference", "paper", "paper_window3"
    ] = "reference_record",
    seeding: Literal["afkmc2", "kmeans", "kmeans++", "paper_nearest"] | None = None,
    seeding_lloyd_steps: int = 0,
    seeding_cells_per_cluster: int | None = None,
    penalty_combine: Literal["product", "sum"] = "product",
    recompute_counts: bool = True,
    low_memory: bool = False,
    rng: SeedLike | RNGLike | None = None,
) -> None:
    """Integrate different experiments using the Harmony algorithm :cite:p:`Korsunsky2019,Patikas2026`.

    Integration adjusts a PCA embedding to reduce experiment-specific technical variation,
    allowing cells from different experiments to be compared.
    This CPU implementation runs in parallel numba kernels.
    Multiple batch variables follow the per-covariate formulation described in the Harmony papers:
    each key is modeled separately instead of combining all keys into one joint category.
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
        Multiple keys are modeled as separate batch variables,
        with one active categorical level per variable and cell.
    basis
        The name of the field in ``adata.obsm`` where the PCA table is stored.
    adjusted_basis
        The name of the field in ``adata.obsm`` where the adjusted PCA table will be stored.
    dtype
        The data type to use for Harmony computation, ``float32`` or ``float64``.
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
        If :data:`None`, uses the value of the reference implementation for the
        chosen ``flavor``: ``4`` for ``"harmony2"``, ``200`` for ``"harmony1"``
        (``200`` for both with the ``"paper"`` stop rules).
    tol_harmony
        Convergence tolerance for the Harmony objective function.
        The algorithm stops when the relative change in objective falls below this value.
        If :data:`None`, ``1e-2`` for ``"harmony2"`` and ``1e-4`` for ``"harmony1"``
        (``1e-4`` for both with the ``"paper"`` stop rules).
    tol_clustering
        Convergence tolerance for the clustering step within each
        Harmony iteration.
        If :data:`None`, ``1e-3`` for ``"harmony2"`` and ``1e-5`` for ``"harmony1"``
        (``1e-5`` for both with the ``"paper"`` stop rules).
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
        A scalar is applied to every key. A sequence may contain one value per key,
        expanded over that key’s categorical levels,
        or one value per categorical level across all keys.
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
        Must be finite and greater than zero.
        Only used with ``flavor="harmony1"``.
    alpha
        Scaling factor for the dynamic per-cluster-per-batch ridge regularization.
        The effective regularization for each cluster-batch pair is ``alpha * E_kb``
        where ``E_kb`` is the expected number of cells.
        Larger values produce more conservative corrections.
        Only used with ``flavor="harmony2"``.
    batch_prune_threshold
        Fraction threshold below which a batch-cluster pair is pruned (correction suppressed).
        When the fraction of a batch's cells assigned to a cluster (``O_kb / N_b``) falls below this threshold,
        that batch-cluster pair receives no correction, preventing spurious adjustments.
        Only used with ``flavor="harmony2"``.
        Set to ``None`` to disable pruning.
    correction_method
        Method for the correction step.
        Only ``"fast"`` is supported: the ridge regression of every cluster is solved exactly,
        eliminating the batch variable with the most levels in closed form,
        for one or several batch keys.
    block_proportion
        Proportion of cells updated per clustering sub-iteration.
        Smaller values produce more stochastic updates.
        Larger values are faster but may converge to different solutions.
    stop_rule
        When the clustering and the outer Harmony loop stop.
        ``"reference_record"`` (default) and ``"reference"`` use the iteration limits and tolerances
        listed under ``max_iter_clustering``, ``tol_clustering`` and ``tol_harmony``,
        and test clustering convergence over a window of 3 iterations.
        ``"reference_record"`` records every round's objective for the outer test;
        ``"reference"`` records it only when the clustering converged.
        ``"paper"`` and ``"paper_window3"`` use the tolerances of :cite:t:`Korsunsky2019`
        (200 clustering iterations, ``1e-5`` and ``1e-4``) for both flavors,
        testing clustering convergence on the last iteration or a window of 3.
    seeding
        How the initial cluster centroids are chosen.
        If :data:`None`, ``"afkmc2"`` for ``"harmony2"`` and ``"kmeans"`` for ``"harmony1"``.
        ``"afkmc2"`` approximates k-means++ seeding with Markov chains :cite:p:`Bachem2016`
        and needs a single pass over the cells.
        ``"kmeans"`` runs :class:`~sklearn.cluster.KMeans` (k-means++ and up to 25 Lloyd iterations).
        ``"kmeans++"`` is exact k-means++ seeding.
        ``"paper_nearest"`` follows the sampler of :cite:t:`Patikas2026`,
        taking the cell with the *highest* ``log(u) / D`` (the paper prints "lowest"),
        where ``D`` is the distance to the nearest centroid chosen so far:
        this draws cells with probability proportional to that distance.
    seeding_lloyd_steps
        Hard spherical k-means iterations applied to the seeds before Harmony starts.
    seeding_cells_per_cluster
        If set, ``seeding="kmeans"`` and the ``seeding_lloyd_steps`` run on a random subsample of
        ``seeding_cells_per_cluster * n_clusters`` cells.
    penalty_combine
        How the diversity penalties of several batch keys are combined in the soft assignment.
        ``"product"`` (default) multiplies the per-key penalty factors, as in :cite:t:`Korsunsky2019`;
        ``"sum"`` adds them, as the assignment rule is printed in :cite:t:`Patikas2026`,
        leaving out levels with ``theta == 0``.
        Identical for a single key.
    recompute_counts
        Rebuild the per-batch cluster counts from the assignments after every clustering iteration.
        If ``False``, the counts are only updated incrementally, which accumulates rounding error
        in float32.
    low_memory
        Drop one cell x component buffer during the iterations
        at the cost of one extra pass at the end; the result is identical.
    rng
        Random number generator or seed for deterministic behavior.

    Returns
    -------
    Updates adata with the field ``adata.obsm[adjusted_basis]``,
    containing principal components adjusted by Harmony
    such that different experiments are integrated.
    """
    from .core import harmonize

    # Resolve flavor into internal flags
    if flavor not in {"harmony1", "harmony2"}:
        msg = f"flavor must be 'harmony1' or 'harmony2', got {flavor!r}."
        raise ValueError(msg)
    if correction_method != "fast":
        msg = f"correction_method must be 'fast', got {correction_method!r}."
        raise ValueError(msg)

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

    adata.obsm[adjusted_basis], _ = harmonize(
        x,
        adata.obs,
        key,
        flavor=flavor,
        n_clusters=n_clusters,
        sigma=sigma,
        theta=theta,
        tau=tau,
        ridge_lambda=ridge_lambda,
        alpha=alpha,
        batch_prune_threshold=batch_prune_threshold,
        block_proportion=block_proportion,
        max_iter_harmony=max_iter_harmony,
        max_iter_clustering=max_iter_clustering,
        tol_harmony=tol_harmony,
        tol_clustering=tol_clustering,
        stop_rule=stop_rule,
        seeding=seeding,
        penalty_combine=penalty_combine,
        seeding_lloyd_steps=seeding_lloyd_steps,
        seeding_cells_per_cluster=seeding_cells_per_cluster,
        recompute_counts=recompute_counts,
        low_memory=low_memory,
        rng=rng,
    )

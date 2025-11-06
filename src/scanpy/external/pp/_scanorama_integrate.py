"""Use Scanorama to integrate cells from different experiments."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..._compat import old_positionals
from ..._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from anndata import AnnData


@old_positionals(
    "basis", "adjusted_basis", "knn", "sigma", "approx", "alpha", "batch_size"
)
@doctest_needs("scanorama")
def scanorama_integrate(
    adata: AnnData,
    key: str,
    *,
    basis: str = "X_pca",
    adjusted_basis: str = "X_scanorama",
    knn: int = 20,
    sigma: float = 15,
    approx: bool = True,
    alpha: float = 0.10,
    batch_size: int = 5000,
    **kwargs,
) -> None:
    """Use Scanorama :cite:p:`Hie2019` to integrate different experiments.

    Scanorama :cite:p:`Hie2019` is an algorithm for integrating single-cell
    data from multiple experiments stored in an AnnData object. This
    function should be run after performing PCA but before computing
    the neighbor graph, as illustrated in the example below.

    This uses the implementation of scanorama_ :cite:p:`Hie2019`.

    .. _scanorama: https://github.com/brianhie/scanorama

    Parameters
    ----------
    adata
        The annotated data matrix.
    key
        The name of the column in ``adata.obs`` that differentiates
        among experiments/batches. Cells from the same batch must be
        contiguously stored in ``adata``.
    basis
        The name of the field in ``adata.obsm`` where the PCA table is
        stored. Defaults to ``'X_pca'``, which is the default for
        ``sc.pp.pca()``.
    adjusted_basis
        The name of the field in ``adata.obsm`` where the integrated
        embeddings will be stored after running this function. Defaults
        to ``X_scanorama``.
    knn
        Number of nearest neighbors to use for matching.
    sigma
        Correction smoothing parameter on Gaussian kernel.
    approx
        Use approximate nearest neighbors with Python ``annoy``;
        greatly speeds up matching runtime.
    alpha
        Alignment score minimum cutoff.
    batch_size
        The batch size used in the alignment vector computation. Useful
        when integrating very large (>100k samples) datasets. Set to
        large value that runs within available memory.
    kwargs
        Any additional arguments will be passed to
        ``scanorama.assemble()``.

    Returns
    -------
    Updates adata with the field ``adata.obsm[adjusted_basis]``,
    containing Scanorama embeddings such that different experiments
    are integrated.

    Example
    -------
    First, load libraries and example dataset, and preprocess.

    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.recipe_zheng17(adata)
    >>> sc.pp.pca(adata)

    We now arbitrarily assign a batch metadata variable to each cell
    for the sake of example, but during real usage there would already
    be a column in ``adata.obs`` giving the experiment each cell came
    from.

    >>> adata.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]

    Finally, run Scanorama. Afterwards, there will be a new table in
    ``adata.obsm`` containing the Scanorama embeddings.

    >>> sce.pp.scanorama_integrate(adata, "batch", verbose=1)
    Processing datasets a <=> b
    >>> "X_scanorama" in adata.obsm
    True

    """
    try:
        import scanorama
    except ImportError as e:
        msg = "\nplease install Scanorama:\n\n\tpip install scanorama"
        raise ImportError(msg) from e

    # Get batch indices in linear time.
    curr_batch = None
    batch_names = []
    name2idx = {}
    for idx in range(adata.X.shape[0]):
        batch_name = adata.obs[key].iat[idx]
        if batch_name != curr_batch:
            curr_batch = batch_name
            if batch_name in batch_names:
                # Contiguous batches important for preserving cell order.
                msg = "Detected non-contiguous batches."
                raise ValueError(msg)
            batch_names.append(batch_name)  # Preserve name order.
            name2idx[batch_name] = []
        name2idx[batch_name].append(idx)

    # Separate batches.
    datasets_dimred = [
        adata.obsm[basis][name2idx[batch_name]] for batch_name in batch_names
    ]

    # Integrate.
    integrated = scanorama.assemble(
        datasets_dimred,  # Assemble in low dimensional space.
        knn=knn,
        sigma=sigma,
        approx=approx,
        alpha=alpha,
        ds_names=batch_names,
        batch_size=batch_size,
        **kwargs,
    )

    adata.obsm[adjusted_basis] = np.concatenate(integrated)

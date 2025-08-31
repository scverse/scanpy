"""Use harmony to integrate cells from different experiments."""

from __future__ import annotations

from collections.abc import Sequence  # noqa: TC003
from typing import TYPE_CHECKING

import numpy as np

from ..._compat import old_positionals
from ..._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from anndata import AnnData


@old_positionals("basis", "adjusted_basis")
@doctest_needs("harmonypy")
def harmony_integrate(
    adata: AnnData,
    key: str | Sequence[str],
    *,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    **kwargs,
):
    """Use harmonypy :cite:p:`Korsunsky2019` to integrate different experiments.

    Harmony :cite:p:`Korsunsky2019` is an algorithm for integrating single-cell
    data from multiple experiments. This function uses the python
    port of Harmony, ``harmonypy``, to integrate single-cell data
    stored in an AnnData object. As Harmony works by adjusting the
    principal components, this function should be run after performing
    PCA but before computing the neighbor graph, as illustrated in the
    example below.

    Parameters
    ----------
    adata
        The annotated data matrix.
    key
        The name of the column in ``adata.obs`` that differentiates
        among experiments/batches. To integrate over two or more covariates,
        you can pass multiple column names as a list. See ``vars_use``
        parameter of the ``harmonypy`` package for more details.
    basis
        The name of the field in ``adata.obsm`` where the PCA table is
        stored. Defaults to ``'X_pca'``, which is the default for
        ``sc.pp.pca()``.
    adjusted_basis
        The name of the field in ``adata.obsm`` where the adjusted PCA
        table will be stored after running this function. Defaults to
        ``X_pca_harmony``.
    kwargs
        Any additional arguments will be passed to
        ``harmonypy.run_harmony()``.

    Returns
    -------
    Updates adata with the field ``adata.obsm[obsm_out_field]``,
    containing principal components adjusted by Harmony such that
    different experiments are integrated.

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

    Finally, run harmony. Afterwards, there will be a new table in
    ``adata.obsm`` containing the adjusted PC's.

    >>> sce.pp.harmony_integrate(adata, "batch")
    >>> "X_pca_harmony" in adata.obsm
    True

    """
    try:
        import harmonypy
    except ImportError as e:
        msg = "\nplease install harmonypy:\n\n\tpip install harmonypy"
        raise ImportError(msg) from e

    X = adata.obsm[basis].astype(np.float64)

    harmony_out = harmonypy.run_harmony(X, adata.obs, key, **kwargs)

    adata.obsm[adjusted_basis] = harmony_out.Z_corr.T

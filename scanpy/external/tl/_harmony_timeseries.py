"""\
Harmony time series for data visualization with augmented affinity matrix at
discrete time points
"""

from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ... import logging as logg


def harmony_timeseries(
    adata: AnnData,
    tp: str,
    n_neighbors: int = 30,
    n_components: Optional[int] = 1000,
    n_jobs: int = -2,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Harmony time series for data visualization with augmented affinity matrix
    at discrete time points [Nowotschin18i]_.

    Harmony time series is a framework for data visualization, trajectory
    detection and interpretation for scRNA-seq data measured at discrete
    time points. Harmony constructs an augmented affinity matrix by augmenting
    the kNN graph affinity matrix with mutually nearest neighbors between
    successive time points. This augmented affinity matrix forms the basis for
    generated a force directed layout for visualization and also serves as input
    for computing the diffusion operator which can be used for trajectory
    detection using Palantir_.

    .. _Palantir: https://github.com/dpeerlab/Palantir

    .. note::
       More information and bug reports `here
       <https://github.com/dpeerlab/Harmony>`__.

    Parameters
    ----------
    adata
        Annotated data matrix of shape n_obs `×` n_vars. Rows correspond to
        cells and columns to genes. Rows represent two or more time points,
        where replicates of the same time point are consecutive in order.
    tp
        key name of observation annotation `.obs` representing time points. Time
        points should be categorical of `dtype=category`. The unique categories for
        the categorical will be used as the time points to construct the timepoint
        connections.
    n_neighbors
        Number of nearest neighbors for graph construction.
    n_components
        Minimum number of principal components to use. Specify `None` to use
        pre-computed components. The higher the value the better to capture 85% of the
        variance.
    n_jobs
        Nearest Neighbors will be computed in parallel using n_jobs.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `.obsm`, `.obsp` and `.uns` with the following:

    **X_harmony** - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.obsm`, dtype `float`)
        force directed layout
    **harmony_aff** - :class:`~scipy.sparse.spmatrix` (:attr:`~anndata.AnnData.obsp`, dtype `float`)
        affinity matrix
    **harmony_aff_aug** - :class:`~scipy.sparse.spmatrix` (:attr:`~anndata.AnnData.obsp`, dtype `float`)
        augmented affinity matrix
    **harmony_timepoint_var** - `str` (:attr:`~anndata.AnnData.uns`)
        The name of the variable passed as `tp`
    **harmony_timepoint_connections** - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.uns`, dtype `str`)
        The links between time points

    Example
    -------

    >>> from itertools import product
    >>> import pandas as pd
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> import scanpy.external as sce

    **Load** `AnnData`

    A sample with real data is available here_.

    .. _here: https://github.com/dpeerlab/Harmony/tree/master/data

    Random data sets of three time points with two replicates each:

    >>> adata_ref = sc.datasets.pbmc3k()
    >>> start = [596, 615, 1682, 1663, 1409, 1432]
    >>> adata = AnnData.concatenate(
    ...     *(adata_ref[i : i + 1000] for i in start),
    ...     join="outer",
    ...     batch_key="sample",
    ...     batch_categories=[f"sa{i}_Rep{j}" for i, j in product((1, 2, 3), (1, 2))],
    ... )
    >>> time_points = adata.obs["sample"].str.split("_", expand=True)[0]
    >>> adata.obs["time_points"] = pd.Categorical(
    ....    time_points, categories=['sa1', 'sa2', 'sa3']
    ... )

    Normalize and filter for highly expressed genes

    >>> sc.pp.normalize_total(adata, target_sum=10000)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    Run harmony_timeseries

    >>> sce.tl.harmony_timeseries(adata, tp="time_points", n_components=500)

    Plot time points:

    >>> sce.pl.harmony_timeseries(adata)

    For further demonstration of Harmony visualizations please follow the notebook
    `Harmony_sample_notebook.ipynb
    <https://github.com/dpeerlab/Harmony/blob/master/notebooks/
    Harmony_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*,
    amongst other things.
    """

    try:
        import harmony
    except ImportError:
        raise ImportError("\nplease install harmony:\n\n\tpip install harmonyTS")

    adata = adata.copy() if copy else adata
    logg.info("Harmony augmented affinity matrix")

    if adata.obs[tp].dtype.name != 'category':
        raise ValueError(f'{tp!r} column does not contain Categorical data')
    timepoints = adata.obs[tp].cat.categories.tolist()
    timepoint_connections = pd.DataFrame(np.array([timepoints[:-1], timepoints[1:]]).T)

    # compute the augmented and non-augmented affinity matrices
    aug_aff, aff = harmony.core.augmented_affinity_matrix(
        data_df=adata.to_df(),
        timepoints=adata.obs[tp],
        timepoint_connections=timepoint_connections,
        n_neighbors=n_neighbors,
        n_jobs=n_jobs,
        pc_components=n_components,
    )

    # Force directed layouts
    layout = harmony.plot.force_directed_layout(aug_aff, adata.obs.index)

    adata.obsm["X_harmony"] = np.asarray(layout)
    adata.obsp["harmony_aff"] = aff
    adata.obsp["harmony_aff_aug"] = aug_aff
    adata.uns["harmony_timepoint_var"] = tp
    adata.uns["harmony_timepoint_connections"] = np.asarray(timepoint_connections)

    return adata if copy else None

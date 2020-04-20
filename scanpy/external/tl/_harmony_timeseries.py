"""\
Harmony time series for data visualization with augmented affinity matrix at
discrete time points
"""

from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ... import logging as logg


def harmony_timeseries(adata: AnnData, tp: str, n_components: Optional[int] = 1000):
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
        key name of observation annotation `.obs` representing time points.
    n_components
        Minimum number of principal components to use. Specify `None` to use
        pre-computed components.

    Returns
    -------
    Updates `.obsm`, `.obsp` and `.uns` with the following:

    `.obsm['X_harmony']`
        force directed layout
    `.obsp['harmony_aff']`
        affinity matrix
    `.obsp['harmony_aff_aug']`
        augmented affinity matrix
    `.uns['harmony_timepoint_var']`
        The name of the variable passed as `tp`
    `.uns['harmony_timepoint_connections']`
        The links between time points

    Example
    -------

    >>> from itertools import product
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
    >>> adata.obs["time_points"] = adata.obs["sample"].str.split("_", expand=True)[0]

    Normalize and filter for highly expressed genes

    >>> sc.pp.normalize_total(adata, target_sum=10000)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    Run harmony_timeseries

    >>> sce.tl.harmony_timeseries(adata, tp="time_points", n_components=None)

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
        raise ImportError(
            "\nplease install harmony:\n\n"
            "\tpip install git+https://github.com/dpeerlab/Harmony.git"
        )

    logg.info("Harmony augmented affinity matrix")

    timepoints = adata.obs[tp].unique().tolist()
    timepoint_connections = pd.DataFrame(np.array([timepoints[:-1], timepoints[1:]]).T)

    # compute the augmented and non-augmented affinity matrices
    aug_aff, aff = harmony.core.augmented_affinity_matrix(
        adata.to_df(), adata.obs[tp], timepoint_connections, pc_components=n_components,
    )

    # Force directed layouts
    layout = harmony.plot.force_directed_layout(aug_aff, adata.obs.index)

    adata.obsm["X_harmony"] = np.asarray(layout)
    adata.obsp["harmony_aff"] = aff
    adata.obsp["harmony_aff_aug"] = aug_aff
    adata.uns["harmony_timepoint_var"] = tp
    adata.uns["harmony_timepoint_connections"] = np.asarray(timepoint_connections)

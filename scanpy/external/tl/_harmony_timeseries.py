"""Harmony time series for data visualization with augmented affinity matrix at
    discrete time points
"""

from typing import Optional

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
    Updates `.obsm` and `.uns` with the following:

    **aff_harmony**
        affinity matrix
    **aff_aug_harmony**
        augmented affinity matrix
    **X_harmony**
        force directed layout
    **timepoint_connections**
        links between time points

    Example
    -------

    >>> import scanpy as sc
    >>> import scanpy.external as sce

    **Load** `AnnData`

    A sample with real data is available here_.

    .. _here: https://github.com/dpeerlab/Harmony/tree/master/data

    Random data sets of three time points with two replicates each:

    >>> adata_ref = sc.datasets.pbmc3k()
    >>> start = [596, 615, 1682, 1663, 1409, 1432]
    >>> adatas = [adata_ref[i : i + 1000] for i in start]
    >>> sample_names = [
            "sa1_Rep1",
            "sa1_Rep2",
            "sa2_Rep1",
            "sa2_Rep2",
            "sa3_Rep1",
            "sa3_Rep2",
        ]
    >>> timepoints = [i.split("_")[0] for i in sample_names]
    >>> for ad, sn, tp in zip(adatas, sample_names, timepoints):
            ad.obs["time_points"] = tp
            ad.obs["sample_name"] = sn
    >>> adata = adatas[0].concatenate(*adatas[1:], join="outer")

    Normalize and filter for highly expressed genes

    >>> sc.pp.normalize_total(adata, target_sum=10000)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    Run harmony_timeseries

    >>> d = sce.tl.harmony_timeseries(
            adata=adata, tp="time_points", n_components=None
        )

    Plot time points:

    >>> from matplotlib import pyplot as plt
    >>> fig, ax = plt.subplots(1, 3, figsize=(12, 5))
    >>> tps = ['sa1', 'sa2', 'sa3']
    >>> for i, tp in enumerate(tps):
            p = sc.pl.scatter(
                adata=adata,
                x='x', y='y',
                color = 'time_points',
                groups=tp, title=tp,
                show=False, ax=ax[i],
                legend_loc='none',
            )
            p.set_axis_off()

    For further demonstration of Harmony visualizations please follow this
    notebook
    `Harmony_sample_notebook.ipynb
    <https://github.com/dpeerlab/Harmony/blob/master/notebooks/
    Harmony_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst
    other things.
    """

    try:
        import harmony
    except ImportError:
        raise ImportError(
            "\nplease install harmony: \n\n\t"
            "git clone https://github.com/dpeerlab/Harmony.git\n\t"
            "cd Harmony\n\t"
            "pip install .\n"
        )

    logg.info("Harmony augmented affinity matrix")

    timepoint_connections = pd.DataFrame(columns=[0, 1])
    index = 0
    timepoints = adata.obs[tp].unique().tolist()
    for i in range(len(timepoints) - 1):
        timepoint_connections.loc[index, :] = timepoints[i : i + 2]
        index += 1
    adata.uns['timepoint_connections'] = timepoint_connections

    # compute the augmented and non-augmented affinity matrices
    (aug_aff, aff,) = harmony.core.augmented_affinity_matrix(
        adata.to_df(),
        adata.obs[tp],
        adata.uns['timepoint_connections'],
        pc_components=n_components,
    )

    # Force directed layouts
    layout = harmony.plot.force_directed_layout(aug_aff, adata.obs.index)

    adata.obsm["X_harmony"] = layout
    adata.uns['aff_harmony'] = aff
    adata.uns['aff_aug_harmony'] = aug_aff

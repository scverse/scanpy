from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ... import logging
from ..._compat import old_positionals
from ..._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from collections.abc import Iterable

    from anndata import AnnData


@old_positionals("branch", "k", "components", "num_waypoints")
@doctest_needs("wishbone")
def wishbone(
    adata: AnnData,
    start_cell: str,
    *,
    branch: bool = True,
    k: int = 15,
    components: Iterable[int] = (1, 2, 3),
    num_waypoints: int | Collection = 250,
):
    """Identify bifurcating developmental trajectories from single-cell data :cite:p:`Setty2016`.

    Wishbone is an algorithm for positioning single cells along bifurcating
    developmental trajectories with high resolution. Wishbone uses multi-dimensional
    single-cell data, such as mass cytometry or RNA-Seq data, as input and orders cells
    according to their developmental progression, and it pinpoints bifurcation points
    by labeling each cell as pre-bifurcation or as one of two post-bifurcation cell
    fates.

    .. note::
       More information and bug reports `here
       <https://github.com/dpeerlab/wishbone>`__.

    Parameters
    ----------
    adata
        Annotated data matrix.
    start_cell
        Desired start cell from `obs_names`.
    branch
        Use True for Wishbone and False for Wanderlust.
    k
        Number of nearest neighbors for graph construction.
    components
        Components to use for running Wishbone.
    num_waypoints
        Number of waypoints to sample.

    Returns
    -------
    Updates `adata` with the following fields:

    `trajectory_wishbone` : (`adata.obs`, dtype `float64`)
        Computed trajectory positions.
    `branch_wishbone` : (`adata.obs`, dtype `int64`)
        Assigned branches.

    Example
    -------

    >>> import scanpy.external as sce
    >>> import scanpy as sc

    **Loading Data and Pre-processing**

    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.normalize_per_cell(adata)
    >>> sc.pp.pca(adata)
    >>> sc.tl.tsne(adata=adata, n_pcs=5, perplexity=30)
    >>> sc.pp.neighbors(adata, n_pcs=15, n_neighbors=10)
    >>> sc.tl.diffmap(adata, n_comps=10)

    **Running Wishbone Core Function**

    Usually, the start cell for a dataset should be chosen based on high expression of
    the gene of interest:

    >>> sce.tl.wishbone(
    ...     adata=adata, start_cell='ACAAGAGACTTATC-1',
    ...     components=[2, 3], num_waypoints=150,
    ... )

    **Visualizing Wishbone results**

    >>> sc.pl.tsne(adata, color=['trajectory_wishbone', 'branch_wishbone'])
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ', 'MALAT1']
    >>> sce.pl.wishbone_marker_trajectory(adata, markers, show=True)

    For further demonstration of Wishbone methods and visualization please follow the
    notebooks in the package `Wishbone_for_single_cell_RNAseq.ipynb
    <https://github.com/dpeerlab/wishbone/tree/master/notebooks>`_.\

    """
    try:
        from wishbone.core import wishbone as c_wishbone
    except ImportError as e:
        msg = "\nplease install wishbone:\n\n\thttps://github.com/dpeerlab/wishbone"
        raise ImportError(msg) from e

    # Start cell index
    s = np.where(adata.obs_names == start_cell)[0]
    if len(s) == 0:
        msg = (
            f"Start cell {start_cell} not found in data. "
            "Please rerun with correct start cell."
        )
        raise RuntimeError(msg)
    if isinstance(num_waypoints, Collection):
        diff = np.setdiff1d(num_waypoints, adata.obs.index)
        if diff.size > 0:
            logging.warning(
                "Some of the specified waypoints are not in the data. "
                "These will be removed"
            )
            num_waypoints = diff.tolist()
    elif num_waypoints > adata.shape[0]:
        msg = (
            "num_waypoints parameter is higher than the number of cells in the "
            "dataset. Please select a smaller number"
        )
        raise RuntimeError(msg)
    s = s[0]

    # Run the algorithm
    components = list(components)
    res = c_wishbone(
        adata.obsm["X_diffmap"][:, components],
        s=s,
        k=k,
        l=k,
        num_waypoints=num_waypoints,
        branch=branch,
    )

    # Assign results
    trajectory = res["Trajectory"]
    trajectory = (trajectory - np.min(trajectory)) / (
        np.max(trajectory) - np.min(trajectory)
    )
    adata.obs["trajectory_wishbone"] = np.asarray(trajectory)

    # branch_ = None
    if branch:
        branches = res["Branches"].astype(int)
        adata.obs["branch_wishbone"] = np.asarray(branches)


def _anndata_to_wishbone(adata: AnnData):
    from wishbone.wb import SCData, Wishbone

    scdata = SCData(adata.to_df())
    scdata.diffusion_eigenvectors = pd.DataFrame(
        adata.obsm["X_diffmap"], index=adata.obs_names
    )
    wb = Wishbone(scdata)
    wb.trajectory = adata.obs["trajectory_wishbone"]
    wb.branch = adata.obs["branch_wishbone"]
    return wb

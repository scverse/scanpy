import warnings

import numpy as np
from scipy import sparse
from typing import Optional

from anndata import AnnData

from ..get import _get_obs_rep


def visium_connectivity(
    adata: AnnData,
    *,
    obsm: str = "spatial",
    key_added: str = "spatial_connectivity",
    n_steps: int = 1,
):
    """
    Creates an adjacency matrix of wells for visium data.

    Params
    ------
    adata
        The AnnData object.
    obsm
        Key in obsm which stores coordinates of wells.
    key_added
        Key in obsp to add spatial graph as.
    n_steps
        How many steps from the starting node should be taken?

    Usage
    -----
    >>> adata = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    >>> adata
    AnnData object with n_obs × n_vars = 2698 × 31053
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial'
    obsm: 'spatial'
    >>> sc.pp.spatial_connectivity(adata)
    >>> adata
    AnnData object with n_obs × n_vars = 2698 × 31053
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial'
    obsm: 'spatial'
    obsp: 'spatial_connectivity'
    """
    adj = _hex_connectivity(_get_obs_rep(adata, obsm=obsm))
    if n_steps > 1:
        adj = walk_nsteps(adj, n_steps - 1)
    adata.obsp[key_added] = adj


def walk_nsteps(adj, n):
    """Expand adjacency matrix adj by walking out n steps from each node."""
    adj = adj.astype(bool)
    cur_step = adj
    result = adj.copy()
    for i in range(n):
        cur_step = adj @ cur_step
        result = result + cur_step
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
        result.setdiag(False)
        result.eliminate_zeros()
    return result


def _hex_connectivity(coords: np.ndarray):
    """
    Given the coordinates of hex based cells from a visium experiment, this
    returns an adjacency matrix for those cells.

    Usage
    -----
    >>> adata.obsp["spatial_connectivity"] = _hex_connectivity(adata.obsm["spatial"])
    """
    from sklearn.neighbors import NearestNeighbors

    N = coords.shape[0]
    dists, row_indices = (
        x.reshape(-1)
        for x in NearestNeighbors(n_neighbors=6, metric="euclidean")
        .fit(coords)
        .kneighbors()
    )
    col_indices = np.repeat(np.arange(N), 6)
    dist_cutoff = np.median(dists) * 1.3  # There's a small amount of sway
    mask = dists < dist_cutoff
    return sparse.csr_matrix(
        (np.ones(mask.sum(), dtype=bool), (row_indices[mask], col_indices[mask])),
        shape=(N, N),
        dtype=bool,
    )


# TODO: Test and document
def radius_neighbors(
    adata: AnnData,
    radius: float,
    *,
    mode: str = "distance",
    obsm: Optional[str] = None,
    layers: Optional[str] = None,
    key_added: Optional[str] = None,
    **kwargs,
):
    """Create a neighbor graph for each cell using chosen represenation.

    Params
    ------
    adata
    radius
        Radius to consider neighbors
    obsm
        Key for representation to find neighbors from.
    layers
        Key for representation to find neighbors from.
    key_added
        Key to add neighbors graph as in obsp.
    """
    from sklearn.neighbors import NearestNeighbors

    if key_added is None:
        key_added = f"{mode}_graph"

    X = _get_obs_rep(adata, obsm=obsm, layers=layers)
    nn = NearestNeighbors(radius=radius, **kwargs)
    nn.fit(X)
    adata.obsp[key_added] = nn.radius_neighbors_graph(mode=mode)

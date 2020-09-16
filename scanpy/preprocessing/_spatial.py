from ..get import _get_obs_rep

import numpy as np
from scipy import sparse
from typing import Optional


def spatial_connectivity(
    adata: "AnnData",
    obsm: str = "spatial",
    key_added: str = "spatial_connectivity",
    n_rings: int = 1,
    n_neigh: int = 6,
    radius: Optional[float] = None,
    coord_type: str = "visium",
):
    """
    Creates graph from spatial coordinates

    Params
    ------
    adata
        The AnnData object.
    obsm
        Key to spatial coordinates.
    key_added
        Key added to connectivity matrix in obsp.
    n_rings
        Number of rings of neighbors for Visium data
    n_neigh
        Number of neighborhoods to consider for non-Visium data
    radius
        Radius of neighbors for non-Visium data
    coord_type
        Type of coordinate system (Visium vs. general coordinates)

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
    coords = _get_obs_rep(adata, obsm=obsm)
    if coord_type == "visium":
        Adj = _build_connectivity(coords, 6, None, True)
        if n_rings > 1:
            # get up to n_rings order connections
            Adj += Adj ** n_rings
            Adj.setdiag(0)
            Adj.eliminate_zeros()
            Adj.data[:] = 1.0
        adata.obsp[key_added] = Adj
    else:
        adata.obsp[key_added] = _spatial_connectivity(coords, n_neigh, radius, False)


def _build_connectivity(
    coords: np.ndarray, n_neigh: int, radius: float, neigh_correct: bool
):
    """
    Build connectivity matrix from spatial coordinates
    """
    from sklearn.neighbors import NearestNeighbors

    N = coords.shape[0]

    tree = NearestNeighbors(
        n_neighbors=n_neigh or 6, radius=radius or 1, metric="euclidean"
    )
    tree.fit(coords)

    if radius is not None:
        results = tree.radius_neighbors()
        row_indices = np.concatenate(results[1])
        lengths = [len(x) for x in results[1]]
        col_indices = np.repeat(np.arange(N), lengths)
    else:
        results = tree.kneighbors()
        dists, row_indices = (result.reshape(-1) for result in results)
        col_indices = np.repeat(np.arange(N), n_neigh or 6)
        if neigh_correct:
            dist_cutoff = np.median(dists) * 1.3  # There's a small amount of sway
            mask = dists < dist_cutoff
            row_indices, col_indices = row_indices[mask], col_indices[mask]

    return sparse.csr_matrix(
        (np.ones(len(row_indices)), (row_indices, col_indices)), shape=(N, N)
    )

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix  # noqa: TID251

if TYPE_CHECKING:
    from anndata import AnnData

_PACKAGE_DATA = Path(__file__).parent / "_data"


def archetype(
    adata: AnnData,
    w=None,
    s=None,
):
    """Compute archetype scores and store them in adata.obsm["archetypes"]."""
    archs = archetype_computer(
        adata.X.transpose(),
        w=w,
        s=s,
    )

    if hasattr(archs, "toarray"):
        archs = archs.toarray()

    adata.obsm["archetypes"] = archs.T


def archetype_computer(
    data,
    w=None,
    s=None,
):
    """
    Compute archetype scores from expression data.

    Parameters
    ----------
    data
        Gene expression matrix.
    w
        Pretrained archetype weight matrix.
    s
        Gene scaling factors.

    Returns
    -------
    h
        Archetype score matrix.
    """
    # Load pretrained archetype weights if not supplied.
    if w is None:
        default_w_path = _PACKAGE_DATA / "Archetypes.csv"
        if not default_w_path.exists():
            msg = (
                "Archetype matrix `w` was not provided and default "
                f"file not found at {default_w_path}"
            )
            raise ValueError(msg)
        w = csr_matrix(np.loadtxt(default_w_path, delimiter=","))

    # Load pretrained gene scaling factors if not supplied.
    if s is None:
        default_s_path = _PACKAGE_DATA / "S.csv"
        if not default_s_path.exists():
            msg = (
                "Scaling vector `s` was not provided and default "
                f"file not found at {default_s_path}"
            )
            raise ValueError(msg)
        s = np.loadtxt(default_s_path, delimiter=",")

    # Reweight genes.
    dd = data / s[:, None]

    # Normalize each observation.
    dd = dd / dd.sum(axis=0)

    return nmf_new_weights(
        dd,
        w,
        6,
    )


def nmf_new_weights(
    v,
    w,
    k,
):
    """
    Compute archetype scores using pretrained weights.

    Parameters
    ----------
    v
        Input expression matrix.
    w
        Pretrained archetype weight matrix.
    k
        Number of archetypes.

    Returns
    -------
    h
        Archetype score matrix.
    """
    n_iter = 1000

    # Fixed random seed for reproducibility.
    rng = np.random.default_rng(0)

    eps = np.finfo(float).eps

    _, n = v.shape

    # Initialize dense archetype matrix.
    h = rng.random((k, n)) + eps

    # Pre-calculate W^T * V outside the loop for speed
    wv = w.T.dot(v)
    if hasattr(wv, "toarray"):
        wv = wv.toarray()

    # Pre-calculate W^T * W outside the loop
    w_tw = w.T.dot(w)
    if hasattr(w_tw, "toarray"):
        w_tw = w_tw.toarray()

    # Multiplicative NMF update while keeping W fixed.
    for _ in range(n_iter):
        wwh = w_tw.dot(h) + eps
        h *= wv / wwh

        # Normalize archetype contributions per cell.
        h_sum = h.sum(axis=0, keepdims=True)
        h_sum[h_sum == 0] = 1.0
        h /= h_sum

    return h


def add_archetypes_to_adata(
    adata: AnnData,
    w=None,
    s=None,
):
    """
    Add archetype scores to an AnnData object.

    Stores results in adata.obsm["archetypes"].
    """
    archs = archetype_computer(
        adata.X.transpose(),
        w=w,
        s=s,
    )

    if hasattr(archs, "toarray"):
        archs = archs.toarray()

    adata.obsm["archetypes"] = archs.T


def compute_hexagon_coordinates(
    archetypes,
):
    """Convert archetype scores into two-dimensional coordinates."""
    # Order archetypes according to the hexagon visualization layout.
    hh = archetypes[:, [1, 2, 0, 4, 3, 5]]

    directions = np.array([0, 1, 2, 3, 4, 5])

    direction_x = np.cos(np.pi / 2 - 2 * np.pi * directions / 6)

    direction_y = np.sin(np.pi / 2 - 2 * np.pi * directions / 6)

    coor1 = hh.dot(direction_x)
    coor2 = hh.dot(direction_y)

    return coor1, coor2


def add_coor_to_adata(
    adata: AnnData,
):
    """
    Add hexagon coordinates to an AnnData object.

    Stores coordinates in adata.obsm["Coor1"] and adata.obsm["Coor2"].
    """
    coor1, coor2 = compute_hexagon_coordinates(adata.obsm["archetypes"])

    adata.obsm["Coor1"] = coor1
    adata.obsm["Coor2"] = coor2

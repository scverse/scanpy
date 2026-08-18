"""Deprecated location of the HashSolo cell hashing demultiplexing method.

See :func:`scanpy.pp.hashsolo`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scverse_misc import Deprecation, deprecated

from ..._utils._doctests import doctest_skipif
from ...preprocessing._hashsolo import _hashsolo

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData


@deprecated(Deprecation("1.13.0", "Use `scanpy.pp.hashsolo` instead."))
@doctest_skipif(reason="Illustrative but not runnable doctest code")
def hashsolo(
    adata: AnnData,
    cell_hashing_columns: Sequence[str],
    *,
    priors: tuple[float, float, float] = (0.01, 0.8, 0.19),
    pre_existing_clusters: str | None = None,
    number_of_noise_barcodes: int | None = None,
    inplace: bool = True,
) -> AnnData | None:
    """Probabilistic demultiplexing of cell hashing data using HashSolo :cite:p:`Bernstein2020`.

    :func:`scanpy.pp.hashsolo` accepts :mod:`anndata.acc` references,
    so the hashing counts may live in :attr:`~anndata.AnnData.obs` *or* in the data matrix,
    and it writes namespaced result fields.

    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    cell_hashing_columns
        `.obs` columns that contain cell hashing counts.
    priors
        Prior probabilities of each hypothesis, in
        the order `[negative, singlet, doublet]`. The default is set to
        `[0.01, 0.8, 0.19]` assuming barcode counts are from cells that
        have passed QC in the transcriptome space, e.g. UMI counts, pct
        mito reads, etc.
    pre_existing_clusters
        The column in `.obs` containing pre-existing cluster assignments
        (e.g. Leiden clusters or cell types, but not batch assignments).
        If provided, demultiplexing will be performed separately for each
        cluster.
    number_of_noise_barcodes
        The number of barcodes used to create the noise distribution.
        Defaults to `len(cell_hashing_columns) - 2`.
    inplace
        Whether to update `adata` in-place or return a copy.

    Returns
    -------
    A copy of the input `adata` if `inplace=False`, otherwise the input
    `adata`. The following fields are added:

    `.obs["most_likely_hypothesis"]`
        Index of the most likely hypothesis, where `0` corresponds to negative,
        `1` to singlet, and `2` to doublet.
    `.obs["cluster_feature"]`
        The cluster assignments used for demultiplexing.
    `.obs["negative_hypothesis_probability"]`
        Probability of the negative hypothesis.
    `.obs["singlet_hypothesis_probability"]`
        Probability of the singlet hypothesis.
    `.obs["doublet_hypothesis_probability"]`
        Probability of the doublet hypothesis.
    `.obs["Classification"]`:
        Classification of the cell, one of the barcodes in `cell_hashing_columns`,
        `"Negative"`, or `"Doublet"`.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.io.read_h5ad("data.h5ad")
    >>> sc.external.pp.hashsolo(adata, ["Hash1", "Hash2", "Hash3"])
    >>> adata.obs.head()

    """
    print(
        "Please cite HashSolo paper:\nhttps://www.cell.com/cell-systems/fulltext/S2405-4712(20)30195-2"
    )
    adata = adata if inplace else adata.copy()
    cell_hashing_columns = list(cell_hashing_columns)
    data = adata.obs[cell_hashing_columns].to_numpy()
    clusters = (
        None
        if pre_existing_clusters is None
        else adata.obs[pre_existing_clusters].to_numpy()
    )
    probs = _hashsolo(
        data,
        priors=priors,
        clusters=clusters,
        number_of_noise_barcodes=number_of_noise_barcodes,
    )
    most_likely_hypothesis = np.argmax(probs, axis=1)

    adata.obs["most_likely_hypothesis"] = most_likely_hypothesis.astype(float)
    adata.obs["cluster_feature"] = 0.0 if clusters is None else clusters
    for i, hypothesis in enumerate(["negative", "singlet", "doublet"]):
        adata.obs[f"{hypothesis}_hypothesis_probability"] = probs[:, i]

    classification = np.asarray(
        np.array(cell_hashing_columns)[np.argmax(data, axis=1)], dtype=object
    )
    classification[most_likely_hypothesis == 0] = "Negative"
    classification[most_likely_hypothesis == 2] = "Doublet"
    adata.obs["Classification"] = classification

    return adata if not inplace else None

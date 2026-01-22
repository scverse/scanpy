"""Metrics which don't quite deserve their own file."""

from __future__ import annotations

from typing import TYPE_CHECKING, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from pandas.api.types import CategoricalDtype

from .._utils import NeighborsView

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    if TYPE_CHECKING:
        from numpy.typing import ArrayLike
    else:  # sphinx-autodoc-typehints will execute the outer block, but end up here:
        ArrayLike = type("ArrayLike", (), dict(__module__="numpy.typing"))

    from .._compat import SpBase


def confusion_matrix(
    orig: pd.Series | np.ndarray | Sequence,
    new: pd.Series | np.ndarray | Sequence,
    data: pd.DataFrame | None = None,
    *,
    normalize: bool = True,
) -> pd.DataFrame:
    """Given an original and new set of labels, create a labelled confusion matrix.

    Parameters `orig` and `new` can either be entries in data or categorical arrays
    of the same size.

    Params
    ------
    orig
        Original labels.
    new
        New labels.
    data
        Optional dataframe to fill entries from.
    normalize
        Should the confusion matrix be normalized?


    Examples
    --------

    .. plot::

        import scanpy as sc; import seaborn as sns
        pbmc = sc.datasets.pbmc68k_reduced()
        cmtx = sc.metrics.confusion_matrix("bulk_labels", "louvain", pbmc.obs)
        sns.heatmap(cmtx)

    """
    from sklearn.metrics import confusion_matrix as _confusion_matrix

    if data is not None:
        if isinstance(orig, str):
            orig = data[orig]
        if isinstance(new, str):
            new = data[new]

    # Coercing so I don't have to deal with it later
    orig, new = pd.Series(orig), pd.Series(new)
    assert len(orig) == len(new)

    unique_labels = pd.unique(np.concatenate((orig.to_numpy(), new.to_numpy())))

    # Compute
    mtx = _confusion_matrix(orig, new, labels=unique_labels)
    if normalize:
        sums = mtx.sum(axis=1)[:, np.newaxis]
        mtx = np.divide(mtx, sums, where=sums != 0)

    # Label
    orig_name = "Original labels" if orig.name is None else orig.name
    new_name = "New Labels" if new.name is None else new.name
    df = pd.DataFrame(
        mtx,
        index=pd.Index(unique_labels, name=orig_name),
        columns=pd.Index(unique_labels, name=new_name),
    )

    # Filter
    if isinstance(orig.dtype, CategoricalDtype):
        orig_idx = pd.Series(orig).cat.categories
    else:
        orig_idx = natsorted(pd.unique(orig))
    if isinstance(new.dtype, CategoricalDtype):
        new_idx = pd.Series(new).cat.categories
    else:
        new_idx = natsorted(pd.unique(new))
    df = df.loc[np.array(orig_idx), np.array(new_idx)]

    return df


@overload
def modularity(
    connectivities: ArrayLike | SpBase,
    /,
    labels: pd.Series | ArrayLike,
    *,
    is_directed: bool,
) -> float: ...


@overload
def modularity(
    adata: AnnData,
    /,
    labels: str | pd.Series | ArrayLike = "leiden",
    *,
    neighbors_key: str | None = None,
    mode: Literal["calculate", "update", "retrieve"] = "calculate",
) -> float: ...


def modularity(
    adata_or_connectivities: AnnData | ArrayLike | SpBase,
    /,
    labels: str | pd.Series | ArrayLike = "leiden",
    *,
    neighbors_key: str | None = None,
    is_directed: bool | None = None,
    mode: Literal["calculate", "update", "retrieve"] = "calculate",
) -> float:
    """Compute the modularity of a graph given its connectivities and labels.

    Parameters
    ----------
    adata_or_connectivities
        The AnnData object containing the data or a weighted adjacency matrix representing the graph.
    labels
        Cluster labels for each node in the graph.
        When `AnnData` is provided, this can be the key in `adata.obs` that contains the clustering labels and defaults to `"leiden"`.
    neighbors_key
        When `AnnData` is provided, the key in `adata.obsp` that contains the connectivities.
    is_directed
        Whether the connectivities are directed or undirected.
        Always `False` if `AnnData` is provided, as connectivities are derived from (symmetric) neighbors.
    mode
        When `AnnData` is provided,
        this controls if the stored modularity is retrieved,
        or if we should calculate it (and optionally update it in `adata.uns[labels]`).

    Returns
    -------
    The modularity of the graph based on the provided clustering.
    """
    if isinstance(adata_or_connectivities, AnnData):
        if is_directed:
            msg = f"Connectivities stored in `AnnData` are undirected, can’t specify `{is_directed=!r}`"
            raise ValueError(msg)
        return modularity_adata(
            adata_or_connectivities,
            labels=labels,
            neighbors_key=neighbors_key,
            mode=mode,
        )
    if isinstance(labels, str):
        msg = "`labels` must be provided as array when passing a connectivities array"
        raise TypeError(msg)
    if is_directed is None:
        msg = "`is_directed` must be provided when passing a connectivities array"
        raise TypeError(msg)
    return modularity_array(
        adata_or_connectivities, labels=labels, is_directed=is_directed
    )


def modularity_adata(
    adata: AnnData,
    /,
    *,
    labels: str | pd.Series | ArrayLike,
    neighbors_key: str | None,
    mode: Literal["calculate", "update", "retrieve"],
) -> float:
    if mode in {"retrieve", "update"} and not isinstance(labels, str):
        msg = "`labels` must be a string when `mode` is `'retrieve'` or `'update'`"
        raise ValueError(msg)
    if mode == "retrieve":
        return adata.uns[labels]["modularity"]

    labels_vec = adata.obs[labels] if isinstance(labels, str) else labels
    connectivities = NeighborsView(adata, neighbors_key)["connectivities"]

    # distances are treated as symmetric, so connectivities as well
    m = modularity(connectivities, labels_vec, is_directed=False)
    if mode == "update":
        adata.uns[labels]["modularity"] = m
    return m


def modularity_array(
    connectivities: ArrayLike | SpBase,
    /,
    *,
    labels: pd.Series | ArrayLike,
    is_directed: bool,
) -> float:
    try:
        import igraph as ig
    except ImportError as e:  # pragma: no cover
        msg = "igraph is require for computing modularity"
        raise ImportError(msg) from e
    igraph_mode = ig.ADJ_DIRECTED if is_directed else ig.ADJ_UNDIRECTED
    graph = ig.Graph.Weighted_Adjacency(connectivities, mode=igraph_mode)
    # cluster labels to integer codes required by igraph
    labels = pd.Categorical(np.asarray(labels)).codes
    return graph.modularity(labels)

"""Metrics which don't quite deserve their own file."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from natsort import natsorted
from pandas.api.types import CategoricalDtype
from scipy.sparse import coo_matrix

from .._compat import SpBase

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData
    from numpy.typing import ArrayLike


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


def modularity(
    connectivities: ArrayLike | SpBase,
    labels: pd.Series | ArrayLike,
    *,
    is_directed: bool,
) -> float:
    """Compute the modularity of a graph given its connectivities and labels.

    Parameters
    ----------
    connectivities:
        Weighted adjacency matrix representing the graph. Can be a dense NumPy array or a sparse CSR matrix.
    labels:
        Cluster labels for each node in the graph.
    is_directed:
        Whether the graph is directed or undirected. If True, the graph is treated as directed; otherwise, it is treated as undirected.

    Returns
    -------
    float
        The modularity of the graph based on the provided clustering.
    """
    try:
        # try to import igraph in case the user wants to calculate modularity
        # not in the main module to avoid import errors
        import igraph as ig
    except ImportError as e:
        msg = "igraph is require for computing modularity"
        raise ImportError(msg) from e
    if isinstance(connectivities, SpBase):
        # check if the connectivities is a sparse matrix
        coo = coo_matrix(connectivities)
        edges = list(zip(coo.row, coo.col, strict=False))
        # converting to the coo format to extract the edges and weights
        # storing only non-zero elements and their indices
        weights = coo.data.tolist()
        graph = ig.Graph(edges=edges, directed=is_directed)
        graph.es["weight"] = weights
    else:
        # if the graph is dense, creates it directly using igraph's adjacency matrix
        dense_array = np.asarray(connectivities)
        igraph_mode = ig.ADJ_DIRECTED if is_directed else ig.ADJ_UNDIRECTED
        graph = ig.Graph.Weighted_Adjacency(dense_array.tolist(), mode=igraph_mode)
    # cluster labels to integer codes required by igraph
    labels = pd.Categorical(np.asarray(labels)).codes

    return graph.modularity(labels)


def modularity_adata(
    adata: AnnData,
    *,
    label_key: str | ArrayLike = "leiden",
    key: str = "neighbors",
    is_directed: bool,
) -> float:
    # default to leiden labels and connectivities as it is more common
    """Compute modularity from an AnnData object using stored graph and clustering labels.

    Parameters
    ----------
    adata:
        The AnnData object containing the data.
    label_key:
        The key in `adata.obs` that contains the clustering labels. Defaults to "leiden".
    key:
        The key in `adata.obsp` that contains the connectivities. Defaults to "neighbors".
    is_directed:
        Whether the graph is directed or undirected. If True, the graph is treated as directed; otherwise, it is treated as undirected.

    Returns
    -------
    float
        The modularity of the graph based on the provided clustering.
    """
    label_array = adata.obs[label_key] if isinstance(label_key, str) else label_key
    connectivities_key = adata.uns[key]["connectivities_key"]
    params = adata.uns[key]["params"]
    connectivities = adata.obsp[connectivities_key]
    is_directed = params["is_directed"]

    if isinstance(connectivities, pd.DataFrame):
        connectivities = connectivities.values

    return modularity(connectivities, label_array, is_directed=is_directed)

"""Metrics which don't quite deserve their own file."""

from __future__ import annotations

from typing import TYPE_CHECKING

import igraph as ig
import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from pandas.api.types import CategoricalDtype

from scanpy._compat import CSRBase

if TYPE_CHECKING:
    from collections.abc import Sequence


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

    unique_labels = pd.unique(
        np.concatenate((np.asarray(orig.values), np.asarray(new.values)))
    )

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


def modularity(connectivities, labels, mode="UNDIRECTED") -> float:
    # default mode is undirected?? can be specified as directed or undirected
    """Compute the modularity of a graph given its connectivities and labels.

    Parameters
    ----------
    connectivities: array-like or sparse matrix
        Weighted adjacency matrix representing the graph. Can be a dense NumPy array or a sparse CSR matrix.
    labels: array-like or pandas.Series
        Cluster labels for each node in the graph.
    mode: str
        The mode of the graph. Can be "UNDIRECTED" or "DIRECTED". Default is "UNDIRECTED".

    Returns
    -------
    float
        The modularity of the graph based on the provided clustering.
    """
    if isinstance(connectivities, CSRBase):
        # Convert sparse matrix to dense format so that igraph can handle it
        # Weighted_Adjacency expects with nested lists or numpy arrays and not sparse matrices
        dense_connectivities = connectivities.toarray()
    else:
        dense_connectivities = connectivities
    # creating igraph graph from the dense connectivity matrix
    graph = ig.Graph.Weighted_Adjacency(dense_connectivities.tolist(), mode=mode)
    if isinstance(labels, pd.Series):
        labels = labels.values
    # making sure labels are in the right format, i.e., a list of integers
    labels = pd.Categorical(np.asarray(labels)).codes
    return graph.modularity(labels)


def modularity_adata(adata: AnnData, labels="leiden", obsp="connectivities") -> float:
    """Compute modularity from an AnnData object using stored graph and clustering labels.

    Parameters
    ----------
    adata: AnnData
        The AnnData object containing the data.
    labels: str or array-like
        The key in adata.obs that contains the cluster labels.
    obsp: str
        The key in adata.obsp that contains the connectivities.

    Returns
    -------
    float
        The modularity of the graph based on the provided clustering.
    """
    # if user passes leiden or louvain as a string, this will get the actual labels
    # from adata.obs
    label_array = adata.obs[labels] if isinstance(labels, str) else labels
    # extracting the connectivities from adata.obsp["connectivities"]
    connectivities = adata.obsp[obsp]
    return modularity(connectivities, label_array)

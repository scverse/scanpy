from numbers import Number
from typing import List, Iterable, Union
from itertools import repeat

import scanpy as sc
import igraph

import pandas as pd
import numpy as np
from scipy.sparse import spmatrix


def _validate_reps(reps, adata):
    reps_out = []
    correct_shape = (adata.n_obs, adata.n_obs)
    for i, rep in enumerate(reps):
        if isinstance(rep, str):
            rep = adata.obsp[rep]
        if isinstance(rep, spmatrix):
            reps_out.append(sc._utils.get_igraph_from_adjacency(rep, directed=True))
        elif isinstance(rep, igraph.Graph):
            if rep.shape != correct_shape:
                raise ValueError(
                    "Graph of invalid shape passed to leiden_multiplex.\n\n"
                    f"Rep {i}'s shape was {rep.shape}, but should be {correct_shape}."
                )
            elif "weight" not in rep.es.attribute_names():
                raise NotImplementedError(
                    "Graphs passed to leiden_multiplex must have edge weights under a 'weight' attribute."
                )
            reps_out.append(rep)
    return reps_out


def _validate_floats(vals, n_reps: int, arg_name: str):
    if isinstance(vals, Number):
        return list(repeat(vals, n_reps))
    vals_out = list(vals)
    if len(vals_out) != n_reps:
        raise ValueError(
            f"Incorrect number of {arg_name}. Expected 1 or {n_reps}, got {len(vals_out)}."
        )
    return vals_out


# TODO: Document better
# TODO: Test better
def leiden_multiplex(
    adata: "sc.AnnData",
    reps: Iterable[Union[str, igraph.Graph, spmatrix]],
    layer_weights: Union[float, Iterable[float]] = 1.0,
    resolutions: Union[float, Iterable[float]] = 1.0,
    key_added: str = "leiden_multiplex",
):
    """Perform a multiplexed clustering on multiple graphs representation of the same data.

    Overview: https://leidenalg.readthedocs.io/en/latest/multiplex.html.

    Params
    ------
    adata
    reps
        Connectivity graphs to use. Possible values are:
        * string key to .obsp
        * sparse matrix of shape (n_obs, n_obs)
        * igraph.Graph
    layer_weights
        Weights to use for each representation in the joint clustering.
    resolutions
        Resolution parameter to use for each representation.
    key_added
        Key in .obs to add this clustering in.

    Usage
    -----
    >>> sc.tl.leiden_multiplex(adata, ["rna_connectivities", "protein_connectivities"])
    >>> sc.pl.umap(adata, color="leiden_multiplex")
    """
    n_reps = len(reps)

    reps = _validate_reps(reps, adata)
    layer_weights = _validate_floats(layer_weights, n_reps, "layer_weights")
    resolutions = _validate_floats(resolutions, n_reps, "resolutions")

    clustering = cluster_joint(reps, layer_weights, resolutions)

    adata.obs[key_added] = pd.Categorical.from_codes(
        clustering, categories=list(map(str, np.unique(clustering))),
    )


def cluster_joint(
    reps: List[igraph.Graph], layer_weights: List[float], resolutions: List[float],
):
    """Actually do the clustering.
    """
    import leidenalg
    partitions = [
        leidenalg.RBConfigurationVertexPartition(
            rep, weights="weight", resolution_parameter=res
        )
        for rep, res in zip(reps, resolutions)
    ]
    optimiser = leidenalg.Optimiser()
    optimiser.optimise_partition_multiplex(partitions, layer_weights=layer_weights)
    return np.array(partitions[0].membership)

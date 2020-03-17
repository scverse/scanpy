import scanpy as sc
import leidenalg

import pandas as pd
import numpy as np

# TODO: Work for n reps
# TODO: Allow passing in the graphs
def leiden_multiplex(
    adata: "sc.AnnData",
    rep1: str,
    rep2: str,
    w1=1.0,
    w2=1.0,
    res1: float = 1.0,
    res2: float = None,
    key_added: str = "leiden_multiplex",
):
    """Perform a multiplexed clustering on multiple graphs representation of the same data.

    Overview: https://leidenalg.readthedocs.io/en/latest/multiplex.html.

    Params
    ------
    adata
    rep
        Adjacency matrix of shape n_obs, n_obs
    w
        Layer weight for graph rep
    res
        Resolution parameter for this graph rep
        If connectivity graphs calculated by umap are being used, you probably don't need to specify different values for this.
    key_added
        Key in .obs to add this clustering in.

    Usage
    -----
    >>> sc.tl.leiden_multiplex(adata, "rna_connectivities", "protein_connectivities")
    >>> sc.pl.umap(adata, color="leiden_multiplex")
    """
    if res2 is None:
        res2 = res1

    g1_rep = sc.get._get_obs_rep(adata, obsp=rep1)
    g2_rep = sc.get._get_obs_rep(adata, obsp=rep2)

    g1 = sc._utils.get_igraph_from_adjacency(g1_rep, directed=True)
    g2 = sc._utils.get_igraph_from_adjacency(g2_rep, directed=True)

    clustering = cluster_joint(g1, g2, res1, res2, w1=w1, w2=w2)

    adata.obs[key_added] = pd.Categorical.from_codes(
        clustering, categories=list(map(str, np.unique(clustering))),
    )


# Internal function, for after everything has been preprocessed and extracted
def cluster_joint(
    g1: "igraph.Graph",
    g2: "igraph.Graph",
    res1: float,
    res2: float,
    w1: float = 1.0,
    w2: float = 1.0,
):
    part1 = leidenalg.RBConfigurationVertexPartition(
        g1, weights="weight", resolution_parameter=res1
    )
    part2 = leidenalg.RBConfigurationVertexPartition(
        g2, weights="weight", resolution_parameter=res2
    )
    optimiser = leidenalg.Optimiser()
    optimiser.optimise_partition_multiplex([part1, part2], layer_weights=[w1, w2])
    return np.array(part1.membership)

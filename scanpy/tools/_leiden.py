from typing import Optional, Tuple, Sequence, Type

import numpy as np
import pandas as pd
from natsort import natsorted
from anndata import AnnData
from scipy import sparse

from .. import utils
from .. import logging as logg
from ..logging import _settings_verbosity_greater_or_equal_than

from ._utils_clustering import rename_groups, restrict_adjacency

try:
    from leidenalg.VertexPartition import MutableVertexPartition
except ImportError:
    class MutableVertexPartition: pass
    MutableVertexPartition.__module__ = 'leidenalg.VertexPartition'


def leiden(
    adata: AnnData,
    resolution: float = 1,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_state: int = 0,
    key_added: Optional[str] = 'leiden',
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = True,
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: Optional[Type[MutableVertexPartition]] = None,
    copy: bool = False,
    **partition_kwargs
):
    """Cluster cells into subgroups [Traag18]_.

    Cluster cells using the Leiden algorithm [Traag18]_, an improved version of the
    Louvain algorithm [Blondel08]_. The Louvain algorithm has been proposed for single-cell
    analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters. Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain `(obs_key, list_of_categories)`.
    key_added
        `adata.obs` key under which to add the cluster labels. (default: `'leiden'`)
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']`.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
    partition_type
        Type of partition to use. Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for :func:`~leidenalg.find_partition`.
    copy
        Whether to copy `adata` or modify it inplace.
    **partition_kwargs
        Any further arguments to pass to `~leidenalg.find_partition`
        (which in turn passes arguments to the `partition_type`).

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of samples) that stores the subgroup id (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['leiden']['params']`
        A dict with the values for the parameters `resolution`, `random_state`, and `n_iterations`.
    """
    try:
        import leidenalg
    except ImportError:
        raise ImportError('Please install the leiden algorithm: `pip3 install leidenalg`.')

    logg.info('running Leiden clustering', r=True)
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        if 'neighbors' not in adata.uns:
            raise ValueError('You need to run `pp.neighbors` first to compute a neighborhood graph.')
        adjacency = adata.uns['neighbors']['connectivities']
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories,
            adjacency
        )
    # convert it to igraph
    g = utils.get_igraph_from_adjacency(adjacency, directed=directed)
    # flip to the default partition type if not overriden by the user
    if partition_type is None:
        partition_type = leidenalg.RBConfigurationVertexPartition
    # prepare find_partition arguments as a dictionary, appending to whatever the user provided
    # it needs to be this way as this allows for the accounting of a None resolution
    # (in the case of a partition variant that doesn't take it on input)
    if partition_kwargs is None:
        partition_kwargs = {}
    if use_weights:
        partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
    partition_kwargs['n_iterations'] = n_iterations
    partition_kwargs['seed'] = random_state
    if resolution is not None:
        partition_kwargs['resolution_parameter'] = resolution
    # clustering proper
    part = leidenalg.find_partition(g, partition_type, **partition_kwargs)
    # store output into adata.obs
    groups = np.array(part.membership)
    if restrict_to is not None:
        groups = rename_groups(
            adata,
            key_added,
            restrict_key,
            restrict_categories,
            restrict_indices,
            groups
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(np.unique(groups).astype('U')),
    )
    # store information on the clustering parameters
    adata.uns['leiden'] = {}
    adata.uns['leiden']['params'] = dict(
        resolution=resolution,
        random_state=random_state,
        n_iterations=n_iterations,
    )
    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    logg.hint('found {} clusters and added\n'
              '    \'{}\', the cluster labels (adata.obs, categorical)'
              .format(len(np.unique(groups)), key_added))
    return adata if copy else None

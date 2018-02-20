from ..data_structs import data_graph


def neighbors(
        adata, n_neighbors=30, knn=True, weights={'probabilities'}, n_jobs=None,
        copy=False):
    """Compute neighborhood graph of observations.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    n_neighbors : `int`, optional (default: 30)
        Number of nearest neighbors in the knn graph. If `knn` is `False`, set
        the Gaussian kernel width to the distance of the `n_neighbors` neighbor.
    knn : `bool`, optional (default: `True`)
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    weights : subset of {`None`, `'probabilities'`, `'distances'`} (default: `'probabilities'`)
        Compute the neighborhood graph with different weights.

    Returns
    -------
    Depending on `weights`, returns the following:
    neighbors : sparse matrix (`adata.uns`, dtype `int8`)
        Unweighted adjacency matrix of the neighborhood graph of the
        `n_neighbors` nearest neighbors among data points.
    neighbors_probabilities : sparse matrix (`adata.uns`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as probabilities.
    neighbors_distances : sparse matrix (`adata.uns`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """
    Dsq, indices, distances = data_graph.get_distance_matrix_and_neighbors(
        adata.X, n_neighbors, sparse=knn, n_jobs=n_jobs)
    adata.uns['neighbors_distancs'] = Dsq

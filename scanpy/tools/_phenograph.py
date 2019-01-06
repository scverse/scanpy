"""Performm clustering using PhenoGraph
"""

from .. import logging as logg


def phenograph( data,
                k=30,
                directed=False,
                prune=False,
                min_cluster_size=10,
                jaccard=True,
                primary_metric='euclidean',
                n_jobs=-1,
                q_tol=1e-3,
                louvain_time_limit=2000,
                nn_method='kdtree'):
    """
        PhenoGraph clustering [Levine15]_.


        :param data:     Numpy ndarray of data to cluster,
                         or sparse matrix of k-nearest neighbor graph
                         If ndarray, n-by-d array of n cells in d dimensions
                         If sparse matrix, n-by-n adjacency matrix
        :param k:        Number of nearest neighbors to use in first
                         step of graph construction
        :param directed: Whether to use a symmetric (default) or
                         asymmetric ("directed") graph
                         The graph construction process produces a directed graph,
                         which is symmetrized by one of two methods (see below)
        :param prune:    Whether to symmetrize by taking the average (prune=False) or
                         product (prune=True) between the graph and its transpose
        :param min_cluster_size:
                         Cells that end up in a cluster smaller than min_cluster_size
                         are considered outliers and are assigned to -1 in the cluster labels
        :param jaccard:  If True, use Jaccard metric between k-neighborhoods to build graph
                         If False, use a Gaussian kernel
        :param primary_metric:
                         Distance metric to define nearest neighbors
                         Options include: {'euclidean','manhattan','correlation','cosine'}.
                         Note that performance will be slower for correlation and cosine
        :param n_jobs:   Nearest Neighbors and Jaccard coefficients will be computed
                         in parallel using n_jobs. If n_jobs=-1, the number of jobs is
                         determined automatically
        :param q_tol:    Tolerance (i.e., precision) for monitoring modularity optimization
        :param louvain_time_limit:
                         Maximum number of seconds to run modularity optimization.
                         If exceeded the best result so far is returned
        :param nn_method:
                         Whether to use brute force or kdtree for nearest neighbor
                         search. For very large high-dimensional data sets, brute force
                         (with parallel computation) performs faster than kdtree

        :return communities:
                        numpy integer array of community assignments for each
                        row in data
        :return graph:  numpy sparse array of the graph that was used for clustering
        :return Q:      the modularity score for communities on graph


        Example
        -------
        
        >>> import scanpy.api as sc
        >>> import numpy as np

        >>> # Cluster and cluster centrolds
        >>> df = np.random.rand(1000,40)
        >>> df.shape
        (1000, 40)
        >>> communities, graph, Q = sc.tl.phenograph(df, k=50)
        Finding 50 nearest neighbors using minkowski metric and 'auto' algorithm
        Neighbors computed in 0.16141605377197266 seconds
        Jaccard graph constructed in 0.7866239547729492 seconds
        Wrote graph to binary file in 0.42542195320129395 seconds
        Running Louvain modularity optimization
        After 1 runs, maximum modularity is Q = 0.223536
        After 2 runs, maximum modularity is Q = 0.235874
        Louvain completed 22 runs in 1.5609488487243652 seconds
        PhenoGraph complete in 2.9466471672058105 seconds

    """
    logg.info('PhenoGraph clustering', r=True)

    try:
        import phenograph
    except ImportError:
        raise ImportError(
            'please install phenograph: '
            'pip3 install git+https://github.com/jacoblevine/phenograph.git')

    communities, graph, Q = phenograph.cluster(
                                        data=data,
                                        k=k,
                                        directed=directed,
                                        prune=prune,
                                        min_cluster_size=min_cluster_size,
                                        jaccard=jaccard,
                                        primary_metric=primary_metric,
                                        n_jobs=n_jobs,
                                        q_tol=q_tol,
                                        louvain_time_limit=louvain_time_limit,
                                        nn_method=nn_method
                                        )

    logg.info('    finished', time=True)

    return communities, graph, Q

from ... import logging as logg
import numpy as np
import scipy.sparse as scp


def trimap(
    adata,
    n_dims=2,
    n_inliers=10,
    n_outliers=5,
    n_random=5,
    metric='euclidean',
    weight_adj=500.0,
    lr=1000.0,
    n_iters=400,
    init_pos=None,
    verbose=False,
    copy=False,
):
    """
    TriMap: Large-scale Dimensionality Reduction Using Triplets

    TriMap is a dimensionality reduction method that uses triplet constraints
    to form a low-dimensional embedding of a set of points. The triplet
    constraints are of the form "point i is closer to point j than point k".
    The triplets are sampled from the high-dimensional representation of the
    points and a weighting scheme is used to reflect the importance of each
    triplet.

    TriMap provides a significantly better global view of the data than the
    other dimensionality reduction methods such t-SNE, LargeVis, and UMAP.
    The global structure includes relative distances of the clusters, multiple
    scales in the data, and the existence of possible outliers. We define a
    global score to quantify the quality of an embedding in reflecting the
    global structure of the data.

    Reference
    ---------

    @article{2019TRIMAP,
     author = {{Amid}, E. and {Warmuth}, M. K.},
     title = "{TriMap: Large-scale Dimensionality Reduction Using Triplets}",
     journal = {ArXiv e-prints},
     archivePrefix = "arXiv",
     eprint = {1910.00204},
     year = 2019,
    }

    Input
    ------

    n_dims: Number of dimensions of the embedding (default = 2)

    n_inliers: Number of inlier points for triplet constraints (default = 10)

    n_outliers: Number of outlier points for triplet constraints (default = 5)

    n_random: Number of random triplet constraints per point (default = 5)

    metric: Distance measure ('euclidean' (default), 'manhattan', 'angular',
    'hamming')

    lr: Learning rate (default = 1000.0)

    n_iters: Number of iterations (default = 400)

    knn_tuple: Use the pre-computed nearest-neighbors information in form of a
    tuple (knn_nbrs, knn_distances) (default = None)

    apply_pca: Apply PCA to reduce the dimensions to 100 if necessary before the
    nearest-neighbor calculation (default = True)

    opt_method: Optimization method ('sd': steepest descent,  'momentum': GD
    with momentum, 'dbd': GD with momentum delta-bar-delta (default))

    verbose: Print the progress report (default = True)

    weight_adj: Adjusting the weights using a non-linear transformation
    (default = 500.0)

    return_seq: Return the sequence of maps recorded every 10 iterations
    (default = False)

    Example
    -------
    
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> pbmc = sce.tl.trimap(pbmc, copy=True)
    >>> sce.pl.trimap(pbmc, color=['bulk_labels'], s=10)

    """

    try:
        from trimap import TRIMAP
    except ImportError:
        raise ImportError(
            '\nplease install trimap: \n\n\tsudo pip install trimap'
        )
    adata = adata.copy() if copy else adata
    start = logg.info('computing TRIMAP')

    if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif isinstance(init_pos, str) and init_pos == 'paga':
        init_coords = get_init_pos_from_paga(adata, random_state=random_state)
    else:
        init_coords = init_pos  # Let TriMap handle it
    if hasattr(init_coords, "dtype"):
        init_coords = check_array(
            init_coords, dtype=np.float32, accept_sparse=False
        )

    if 'X_pca' in adata.obsm:
        n_dim_pca = adata.obsm['X_pca'].shape[1]
        X = adata.obsm['X_pca'][:, : min(n_dim_pca, 100)]
    else:
        X = adata.X
        if scp.issparse(X):
            raise ValueError(
                'trimap currently does not support sparse matrices. please'
                'use a dense matrix or apply pca first.'
            )
    if metric not in ['euclidean', 'manhattan', 'angular', 'hamming']:
        raise ValueError(
            'metric \'' + metric + '\' is currently not supported '
            'by trimap. acceptable values are \'euclidean\', \'manhattan\', '
            '\'angular\', and \'hamming\'.'
        )
    X_trimap = TRIMAP(
        n_dims=n_dims,
        n_inliers=n_inliers,
        n_outliers=n_outliers,
        n_random=n_random,
        lr=lr,
        distance=metric,
        weight_adj=weight_adj,
        n_iters=n_iters,
        verbose=verbose,
    ).fit_transform(X, init=init_coords)
    adata.obsm['X_trimap'] = X_trimap
    logg.info(
        '    finished',
        time=start,
        deep=('added\n' "    'X_trimap', TriMap coordinates (adata.obsm)"),
    )
    return adata if copy else None

"""Embed high-dimensional data using PHATE
"""


def PHATE(
        n_components=2,
        k=15,
        a=10,
        alpha_decay=None,
        n_landmark=2000,
        t='auto',
        potential_method='log',
        n_pca=100,
        knn_dist='euclidean',
        mds_dist='euclidean',
        mds='metric',
        n_jobs=1,
        random_state=None,
        verbose=True):
    """PHATE operator which performs dimensionality reduction.

    Potential of Heat-diffusion for Affinity-based Trajectory Embedding
    (PHATE). Embeds high dimensional single-cell data into two or three
    dimensions for visualization of biological progressions as described
    in [Moon17]_.

    PHATE operator is used with standard `sklearn` estimator methods:
    `fit(X)`, `transform(X)` and `fit_transform(X)`. Read more in the
    `PHATE Documentation <https://phate.readthedocs.io/>`_.

    Parameters
    ----------

    n_components : int, optional, default: 2
        number of dimensions in which the data will be embedded

    k : int, optional, default: 15
        number of nearest neighbors on which to build kernel

    a : int, optional, default: None
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used

    alpha_decay : boolean, default: None
        forces the use of alpha decaying kernel
        If None, alpha decaying kernel is used for small inputs
        (n_samples < n_landmark) and not used otherwise

    n_landmark : int, optional, default: 2000
        number of landmarks to use in fast PHATE

    t : int, optional, default: 'auto'
        power to which the diffusion operator is powered
        sets the level of diffusion

    potential_method : string, optional, default: 'log'
        choose from ['log', 'sqrt']
        which transformation of the diffusional operator is used
        to compute the diffusion potential

    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        log(n_samples) time.

    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph

    mds_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for MDS

    mds : string, optional, default: 'metric'
        choose from ['classic', 'metric', 'nonmetric']
        which MDS algorithm is used for dimensionality reduction

    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used

    random_state : integer or `numpy.RandomState`, optional
        The generator used to initialize SMACOF (metric, nonmetric) MDS
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator

    verbose : boolean, optional
        If true, print status messages

    Attributes
    ----------

    X : array-like, shape=[n_samples, n_dimensions]

    embedding : array-like, shape=[n_samples, n_components]
        Stores the position of the dataset in the embedding space

    diff_op : array-like, shape=[n_samples, n_samples] or [n_landmarks, n_landmarks]
        The diffusion operator fit on the input data

    diff_potential : array-like, shape=[n_samples, n_samples]
        Precomputed diffusion potential

    landmark_transitions : array-like, shape=[n_samples, n_landmarks]
        Transition matrix between input data and landmarks,
        if `n_landmark` is set, otherwise `None`

    Examples
    --------
    >>> import scanpy.api as sc
    >>> import matplotlib.pyplot as plt
    >>> from anndata import AnnData
    >>> import phate
    >>> tree_data, tree_clusters = phate.tree.gen_dla(n_dim=100,
                                                      n_branch=20,
                                                      branch_length=100)
    >>> tree_data.shape
    (2000, 100)
    >>> adata = AnnData(tree_data)
    >>> phate_operator = sc.tl.PHATE(k=5, a=20, t=100)
    >>> tree_phate = phate_operator.fit_transform(adata)
    >>> tree_phate.shape
    (2000, 2)
    >>> plt.scatter(tree_phate[:,0], tree_phate[:,1], c=tree_clusters)
    >>> plt.show()
    """
    try:
        import phate as ph
    except ImportError:
        raise ImportError('You need to install the package `phate`: please run'
                          ' `pip install - -user phate` in a terminal.')

    return ph.PHATE(
        n_components=n_components,
        k=k,
        a=a,
        alpha_decay=alpha_decay,
        n_landmark=n_landmark,
        t=t,
        potential_method=potential_method,
        n_pca=n_pca,
        knn_dist=knn_dist,
        mds_dist=mds_dist,
        mds=mds,
        n_jobs=n_jobs,
        random_state=random_state,
        verbose=verbose,
    )


def run_phate(
        adata,
        n_components=2,
        k=15,
        a=10,
        alpha_decay=None,
        n_landmark=2000,
        t='auto',
        potential_method='log',
        n_pca=100,
        knn_dist='euclidean',
        mds_dist='euclidean',
        mds='metric',
        n_jobs=1,
        random_state=None,
        verbose=True):
    """Runs PHATE dimensionality reduction.

    Embeds high dimensional single-cell data into two or three dimensions for
    visualization of biological progressions. Potential of Heat-diffusion for
    Affinity-based Trajectory Embedding (PHATE) [Moon17]_

    Parameters
    ----------
    data : AnnData or array-like [n_samples, n_dimensions]
        input data with `n_samples` samples and `n_dimensions`
        dimensions. Accepted data types: `numpy.ndarray`,
        `scipy.sparse.spmatrix`, `pd.DataFrame`, `anndata.AnnData`

    n_components : int, optional, default: 2
        number of dimensions in which the data will be embedded

    k : int, optional, default: 15
        number of nearest neighbors on which to build kernel

    a : int, optional, default: None
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used

    alpha_decay : boolean, default: None
        forces the use of alpha decaying kernel
        If None, alpha decaying kernel is used for small inputs
        (n_samples < n_landmark) and not used otherwise

    n_landmark : int, optional, default: 2000
        number of landmarks to use in fast PHATE

    t : int, optional, default: 'auto'
        power to which the diffusion operator is powered
        sets the level of diffusion

    potential_method : string, optional, default: 'log'
        choose from ['log', 'sqrt']
        which transformation of the diffusional operator is used
        to compute the diffusion potential

    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        log(n_samples) time.

    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph

    mds_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for MDS

    mds : string, optional, default: 'metric'
        choose from ['classic', 'metric', 'nonmetric']
        which MDS algorithm is used for dimensionality reduction

    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used

    random_state : integer or `numpy.RandomState`, optional
        The generator used to initialize SMACOF (metric, nonmetric) MDS
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator

    verbose : boolean, optional
        If true, print status messages

    Returns
    -------

    embedding : array-like, shape=[n_samples, n_components]
        Stores the position of the dataset in the embedding space

    Examples
    --------
    >>> import scanpy.api as sc
    >>> import matplotlib.pyplot as plt
    >>> from anndata import AnnData
    >>> import phate
    >>> tree_data, tree_clusters = phate.tree.gen_dla(n_dim=100,
                                                      n_branch=20,
                                                      branch_length=100)
    >>> tree_data.shape
    (2000, 100)
    >>> adata = AnnData(tree_data)
    >>> tree_phate = sc.tl.run_phate(adata, k=5, a=20, t=150)
    >>> tree_phate.shape
    (2000, 2)
    >>> plt.scatter(tree_phate[:,0], tree_phate[:,1], c=tree_clusters)
    >>> plt.show()
    """
    try:
        import phate as ph
    except ImportError:
        raise ImportError(
            'You need to install the package `phate`: please run `pip install --user phate` in a terminal.')

    return ph.PHATE(
        n_components=n_components,
        k=k,
        a=a,
        alpha_decay=alpha_decay,
        n_landmark=n_landmark,
        t=t,
        potential_method=potential_method,
        n_pca=n_pca,
        knn_dist=knn_dist,
        mds_dist=mds_dist,
        mds=mds,
        n_jobs=n_jobs,
        random_state=random_state,
        verbose=verbose,
    ).fit_transform(adata)

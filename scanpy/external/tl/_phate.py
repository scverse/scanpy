"""\
Embed high-dimensional data using PHATE
"""
from typing import Optional, Union

from anndata import AnnData

from ..._settings import settings
from ..._compat import Literal
from ... import logging as logg
from ..._utils import AnyRandom


def phate(
    adata: AnnData,
    n_components: int = 2,
    k: int = 5,
    a: int = 15,
    n_landmark: int = 2000,
    t: Union[int, str] = 'auto',
    gamma: float = 1.0,
    n_pca: int = 100,
    knn_dist: str = 'euclidean',
    mds_dist: str = 'euclidean',
    mds: Literal['classic', 'metric', 'nonmetric'] = 'metric',
    n_jobs: Optional[int] = None,
    random_state: AnyRandom = None,
    verbose: Union[bool, int, None] = None,
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """\
    PHATE [Moon17]_.

    Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
    embeds high dimensional single-cell data into two or three dimensions for
    visualization of biological progressions.

    For more information and access to the object-oriented interface, read the
    `PHATE documentation <https://phate.readthedocs.io/>`__.  For
    tutorials, bug reports, and R/MATLAB implementations, visit the `PHATE
    GitHub page <https://github.com/KrishnaswamyLab/PHATE/>`__. For help
    using PHATE, go `here <https://krishnaswamylab.org/get-help>`__.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_components
        number of dimensions in which the data will be embedded
    k
        number of nearest neighbors on which to build kernel
    a
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used
    n_landmark
        number of landmarks to use in fast PHATE
    t
        power to which the diffusion operator is powered
        sets the level of diffusion. If 'auto', t is selected
        according to the knee point in the Von Neumann Entropy of
        the diffusion operator
    gamma
        Informational distance constant between -1 and 1.
        `gamma=1` gives the PHATE log potential, `gamma=0` gives
        a square root potential.
    n_pca
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        log(n_samples) time.
    knn_dist
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph
    mds_dist
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for MDS
    mds
        Selects which MDS algorithm is used for dimensionality reduction.
    n_jobs
        The number of jobs to use for the computation.
        If `None`, `sc.settings.n_jobs` is used.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used
    random_state
        Random seed. Defaults to the global `numpy` random number generator
    verbose
        If `True` or an `int`/`Verbosity` ≥ 2/`hint`, print status messages.
        If `None`, `sc.settings.verbosity` is used.
    copy
        Return a copy instead of writing to `adata`.
    kwargs
        Additional arguments to `phate.PHATE`

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_phate** : `np.ndarray`, (`adata.obs`, shape=[n_samples, n_components], dtype `float`)
        PHATE coordinates of data.

    Examples
    --------
    >>> from anndata import AnnData
    >>> import scanpy.external as sce
    >>> import phate
    >>> tree_data, tree_clusters = phate.tree.gen_dla(
    ...     n_dim=100,
    ...     n_branch=20,
    ...     branch_length=100,
    ... )
    >>> tree_data.shape
    (2000, 100)
    >>> adata = AnnData(tree_data)
    >>> sce.tl.phate(adata, k=5, a=20, t=150)
    >>> adata.obsm['X_phate'].shape
    (2000, 2)
    >>> sce.pl.phate(adata)
    """
    start = logg.info('computing PHATE')
    adata = adata.copy() if copy else adata
    verbosity = settings.verbosity if verbose is None else verbose
    verbose = verbosity if isinstance(verbosity, bool) else verbosity >= 2
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    try:
        import phate
    except ImportError:
        raise ImportError(
            'You need to install the package `phate`: please run `pip install '
            '--user phate` in a terminal.'
        )
    X_phate = phate.PHATE(
        n_components=n_components,
        k=k,
        a=a,
        n_landmark=n_landmark,
        t=t,
        gamma=gamma,
        n_pca=n_pca,
        knn_dist=knn_dist,
        mds_dist=mds_dist,
        mds=mds,
        n_jobs=n_jobs,
        random_state=random_state,
        verbose=verbose,
        **kwargs,
    ).fit_transform(adata)
    # update AnnData instance
    adata.obsm['X_phate'] = X_phate  # annotate samples with PHATE coordinates
    logg.info(
        '    finished',
        time=start,
        deep=('added\n' "    'X_phate', PHATE coordinates (adata.obsm)"),
    )
    return adata if copy else None

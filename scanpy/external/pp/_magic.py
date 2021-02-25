"""\
Denoise high-dimensional data using MAGIC
"""
from typing import Union, Sequence, Optional

from anndata import AnnData
from legacy_api_wrap import legacy_api
from packaging import version

from ... import logging as logg
from ..._settings import settings
from ..._compat import Literal
from ..._utils import AnyRandom


MIN_VERSION = "2.0"


@legacy_api('k', 'a')
def magic(
    adata: AnnData,
    name_list: Union[Literal['all_genes', 'pca_only'], Sequence[str], None] = None,
    *,
    knn: int = 5,
    decay: Optional[float] = 1,
    knn_max: Optional[int] = None,
    t: Union[Literal['auto'], int] = 3,
    n_pca: Optional[int] = 100,
    solver: Literal['exact', 'approximate'] = 'exact',
    knn_dist: str = 'euclidean',
    random_state: AnyRandom = None,
    n_jobs: Optional[int] = None,
    verbose: bool = False,
    copy: Optional[bool] = None,
    **kwargs,
) -> Optional[AnnData]:
    """\
    Markov Affinity-based Graph Imputation of Cells (MAGIC) API [vanDijk18]_.

    MAGIC is an algorithm for denoising and transcript recover of single cells
    applied to single-cell sequencing data. MAGIC builds a graph from the data
    and uses diffusion to smooth out noise and recover the data manifold.

    The algorithm implemented here has changed primarily in two ways
    compared to the algorithm described in [vanDijk18]_. Firstly, we use
    the adaptive kernel described in Moon et al, 2019 [Moon17]_ for
    improved stability. Secondly, data diffusion is applied
    in the PCA space, rather than the data space, for speed and
    memory improvements.

    More information and bug reports
    `here <https://github.com/KrishnaswamyLab/MAGIC>`__. For help, visit
    <https://krishnaswamylab.org/get-help>.

    Parameters
    ----------
    adata
        An anndata file with `.raw` attribute representing raw counts.
    name_list
        Denoised genes to return. The default `'all_genes'`/`None`
        may require a large amount of memory if the input data is sparse.
        Another possibility is `'pca_only'`.
    knn
        number of nearest neighbors on which to build kernel.
    decay
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used.
    knn_max
        maximum number of nearest neighbors with nonzero connection.
        If `None`, will be set to 3 * `knn`.
    t
        power to which the diffusion operator is powered.
        This sets the level of diffusion. If 'auto', t is selected
        according to the Procrustes disparity of the diffused data.
    n_pca
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time. If `None`, no PCA is performed.
    solver
        Which solver to use. "exact" uses the implementation described
        in van Dijk et al. (2018) [vanDijk18]_. "approximate" uses a faster
        implementation that performs imputation in the PCA space and then
        projects back to the gene space. Note, the "approximate" solver may
        return negative values.
    knn_dist
        recommended values: 'euclidean', 'cosine', 'precomputed'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph. If 'precomputed',
        `data` should be an n_samples x n_samples distance or
        affinity matrix.
    random_state
        Random seed. Defaults to the global `numpy` random number generator.
    n_jobs
        Number of threads to use in training. All cores are used by default.
    verbose
        If `True` or an integer `>= 2`, print status messages.
        If `None`, `sc.settings.verbosity` is used.
    copy
        If true, a copy of anndata is returned. If `None`, `copy` is True if
        `genes` is not `'all_genes'` or `'pca_only'`. `copy` may only be False
        if `genes` is `'all_genes'` or `'pca_only'`, as the resultant data
        will otherwise have different column names from the input data.
    kwargs
        Additional arguments to `magic.MAGIC`.

    Returns
    -------
    If `copy` is True, AnnData object is returned.

    If `subset_genes` is not `all_genes`, PCA on MAGIC values of cells are
    stored in `adata.obsm['X_magic']` and `adata.X` is not modified.

    The raw counts are stored in `.raw` attribute of AnnData object.

    Examples
    --------
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> adata = sc.datasets.paul15()
    >>> sc.pp.normalize_per_cell(adata)
    >>> sc.pp.sqrt(adata)  # or sc.pp.log1p(adata)
    >>> adata_magic = sce.pp.magic(adata, name_list=['Mpo', 'Klf1', 'Ifitm1'], knn=5)
    >>> adata_magic.shape
    (2730, 3)
    >>> sce.pp.magic(adata, name_list='pca_only', knn=5)
    >>> adata.obsm['X_magic'].shape
    (2730, 100)
    >>> sce.pp.magic(adata, name_list='all_genes', knn=5)
    >>> adata.X.shape
    (2730, 3451)
    """

    try:
        from magic import MAGIC, __version__
    except ImportError:
        raise ImportError(
            'Please install magic package via `pip install --user '
            'git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python`'
        )
    else:
        if not version.parse(__version__) >= version.parse(MIN_VERSION):
            raise ImportError(
                'scanpy requires magic-impute >= '
                f'v{MIN_VERSION} (detected: v{__version__}). '
                'Please update magic package via `pip install --user '
                '--upgrade magic-impute`'
            )

    start = logg.info('computing MAGIC')
    all_or_pca = isinstance(name_list, (str, type(None)))
    if all_or_pca and name_list not in {"all_genes", "pca_only", None}:
        raise ValueError(
            "Invalid string value for `name_list`: "
            "Only `'all_genes'` and `'pca_only'` are allowed."
        )
    if copy is None:
        copy = not all_or_pca
    elif not all_or_pca and not copy:
        raise ValueError(
            "Can only perform MAGIC in-place with `name_list=='all_genes' or "
            f"`name_list=='pca_only'` (got {name_list}). Consider setting "
            "`copy=True`"
        )
    adata = adata.copy() if copy else adata
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs

    X_magic = MAGIC(
        knn=knn,
        decay=decay,
        knn_max=knn_max,
        t=t,
        n_pca=n_pca,
        solver=solver,
        knn_dist=knn_dist,
        random_state=random_state,
        n_jobs=n_jobs,
        verbose=verbose,
        **kwargs,
    ).fit_transform(adata, genes=name_list)
    logg.info(
        '    finished',
        time=start,
        deep=(
            "added\n    'X_magic', PCA on MAGIC coordinates (adata.obsm)"
            if name_list == "pca_only"
            else ''
        ),
    )
    # update AnnData instance
    if name_list == "pca_only":
        # special case – update adata.obsm with smoothed values
        adata.obsm["X_magic"] = X_magic.X
    elif copy:
        # just return X_magic
        X_magic.raw = adata
        adata = X_magic
    else:
        # replace data with smoothed data
        adata.raw = adata
        adata.X = X_magic.X

    if copy:
        return adata

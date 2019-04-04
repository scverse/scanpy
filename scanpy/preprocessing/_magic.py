"""Denoise high-dimensional data using MAGIC
"""

from .._settings import settings
from .. import logging as logg
from ..logging import _settings_verbosity_greater_or_equal_than


def magic(adata,
          name_list=None,
          k=10,
          a=15,
          t='auto',
          n_pca=100,
          knn_dist='euclidean',
          random_state=None,
          n_jobs=None,
          verbose=False,
          copy=None,
          **kwargs):
    """Markov Affinity-based Graph Imputation of Cells (MAGIC) API [vanDijk18]_.

    MAGIC is an algorithm for denoising and transcript recover of single cells
    applied to single-cell sequencing data. MAGIC builds a graph from the data
    and uses diffusion to smooth out noise and recover the data manifold.

    More information and bug reports
    `here <https://github.com/KrishnaswamyLab/MAGIC>`__. For help, visit
    <https://krishnaswamylab.org/get-help>.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        An anndata file with `.raw` attribute representing raw counts.
    name_list : `list`, `'all_genes'`, or `'pca_only'`, optional (default: `'all_genes'`)
        Denoised genes to return. Default is all genes, but this
        may require a large amount of memory if the input data is sparse.
    k : int, optional, default: 10
        number of nearest neighbors on which to build kernel
    a : int, optional, default: 15
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used
    t : int, optional, default: 'auto'
        power to which the diffusion operator is powered.
        This sets the level of diffusion. If 'auto', t is selected
        according to the Procrustes disparity of the diffused data
    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time.
    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean', 'cosine', 'precomputed'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph. If 'precomputed',
        `data` should be an n_samples x n_samples distance or
        affinity matrix
    random_state : `int`, `numpy.RandomState` or `None`, optional (default: `None`)
        Random seed. Defaults to the global `numpy` random number generator
    n_jobs : `int` or None, optional. Default: None
        Number of threads to use in training. All cores are used by default.
    verbose : `bool`, `int` or `None`, optional (default: `sc.settings.verbosity`)
        If `True` or an integer `>= 2`, print status messages.
        If `None`, `sc.settings.verbosity` is used.
    copy : `bool` or `None`, optional. Default: `None`.
        If true, a copy of anndata is returned. If `None`, `copy` is True if
        `genes` is not `'all_genes'` or `'pca_only'`. `copy` may only be False
        if `genes` is `'all_genes'` or `'pca_only'`, as the resultant data
        will otherwise have different column names from the input data.
    kwargs : additional arguments to `magic.MAGIC`

    Returns
    -------
    If `copy` is True, AnnData object is returned.

    If `subset_genes` is not `all_genes`, PCA on MAGIC values of cells are stored in
    `adata.obsm['X_magic']` and `adata.X` is not modified.

    The raw counts are stored in `.raw` attribute of AnnData object.

    Examples
    --------
    >>> import scanpy.api as sc
    >>> import magic
    >>> adata = sc.datasets.paul15()
    >>> sc.pp.normalize_per_cell(adata)
    >>> sc.pp.sqrt(adata)  # or sc.pp.log1p(adata)
    >>> adata_magic = sc.pp.magic(adata, name_list=['Mpo', 'Klf1', 'Ifitm1'], k=5)
    >>> adata_magic.shape
    (2730, 3)
    >>> sc.pp.magic(adata, name_list='pca_only', k=5)
    >>> adata.obsm['X_magic'].shape
    (2730, 100)
    >>> sc.pp.magic(adata, name_list='all_genes', k=5)
    >>> adata.X.shape
    (2730, 3451)
    """

    try:
        from magic import MAGIC
    except ImportError:
        raise ImportError(
            'Please install magic package via `pip install --user '
            'git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python`')

    logg.info('computing PHATE', r=True)
    needs_copy = not (name_list is None or
                      (isinstance(name_list, str) and
                       name_list in ["all_genes", "pca_only"]))
    if copy is None:
        copy = needs_copy
    elif needs_copy and not copy:
        raise ValueError(
            "Can only perform MAGIC in-place with `name_list=='all_genes' or "
            "`name_list=='pca_only'` (got {}). Consider setting "
            "`copy=True`".format(name_list))
    adata = adata.copy() if copy else adata
    verbose = settings.verbosity if verbose is None else verbose
    if isinstance(verbose, (str, int)):
        verbose = _settings_verbosity_greater_or_equal_than(2)
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs

    X_magic = MAGIC(k=k,
                    a=a,
                    t=t,
                    n_pca=n_pca,
                    knn_dist=knn_dist,
                    random_state=random_state,
                    n_jobs=n_jobs,
                    verbose=verbose,
                    **kwargs).fit_transform(adata,
                                            genes=name_list)
    logg.info('    finished', time=True,
              end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    # update AnnData instance
    if name_list == "pca_only":
        # special case - update adata.obsm with smoothed values
        adata.obsm["X_magic"] = X_magic.X
        logg.hint('added\n'
                  '    \'X_magic\', PCA on MAGIC coordinates (adata.obsm)')
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

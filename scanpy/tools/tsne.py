import numpy as np
from ..tools.pca import pca
from .. import settings
from .. import logging as logg


def tsne(adata, n_pcs=50, perplexity=30, early_exaggeration=12,
         learning_rate=1000, random_state=0, use_fast_tsne=True,
         recompute_pca=False, n_jobs=None, copy=False):
    """t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.

    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been
    proposed for visualizating single-cell data by [Amir13]_. Here, by default,
    we use the implementation of *scikit-learn* [Pedregosa11]_. You can achieve
    a huge speedup and better convergence if you install `*Multicore-tSNE*
    <https://github.com/DmitryUlyanov/Multicore-TSNE>`__ by [Ulyanov16]_, which
    will be automatically detected by Scanpy.

    Parameters
    ----------
    adata : `~scanpy.api.AnnData`
        Annotated data matrix.
    n_pcs : `int`, optional (default: 50)
        Number of principal components in preprocessing PCA. Set to 0 if you
        do not want preprocessing with PCA.
    perplexity : `float`, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    early_exaggeration : `float`, optional (default: 12.0)
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate : `float`, optional (default: 1000)
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state : `int` or `None`, optional (default: 0)
        Change this to use different intial states for the optimization. If `None`,
        the initial state is not reproducible.
    use_fast_tsne : `bool`, optional (default: `True`)
        Use the MulticoreTSNE package by D. Ulyanov if it is installed.
    n_jobs : `int` or `None` (default: `sc.settings.n_jobs`)
        Number of jobs.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_tsne : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    """
    logg.info('computing tSNE', r=True)
    adata = adata.copy() if copy else adata
    # preprocessing by PCA
    if (n_pcs > 0
        and 'X_pca' in adata.obsm_keys()
        and adata.obsm['X_pca'].shape[1] >= n_pcs
        and not recompute_pca):
        X = adata.obsm['X_pca'][:, :n_pcs]
        logg.info('    using \'X_pca\' with n_pcs = {} for tSNE'
                  .format(n_pcs))
    else:
        if n_pcs > 0 and adata.X.shape[1] > n_pcs:
            logg.info('    computing \'X_pca\' with n_pcs = {}'.format(n_pcs))
            logg.hint('avoid this by setting n_pcs = 0')
            X = pca(adata.X, random_state=random_state, n_comps=n_pcs)
            adata.obsm['X_pca'] = X
        else:
            X = adata.X
            logg.info('    using data matrix X directly (no PCA)')
    # params for sklearn
    params_sklearn = {'perplexity': perplexity,
                      'random_state': random_state,
                      'verbose': max(0, settings.verbosity-3),
                      'early_exaggeration': early_exaggeration,
                      'learning_rate': learning_rate,
                      }
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    # deal with different tSNE implementations
    multicore_failed = True
    if n_jobs >= 1 and use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE
            tsne = TSNE(n_jobs=n_jobs, **params_sklearn)
            logg.info('    using the "MulticoreTSNE" package by Ulyanov (2017)')
            X_tsne = tsne.fit_transform(X.astype(np.float64))
            multicore_failed = False
        except ImportError:
            logg.warn('Consider installing the package MulticoreTSNE '
                      '(https://github.com/DmitryUlyanov/Multicore-TSNE). '
                      'Even for n_jobs=1 this speeds up the computation considerably '
                      'and might yield better converged results.')
            pass
    if multicore_failed:
        from sklearn.manifold import TSNE
        from . import _tsne_fix   # fix by D. DeTomaso for sklearn < 0.19
        # unfortunately, sklearn does not allow to set a minimum number of iterations for barnes-hut tSNE
        tsne = TSNE(**params_sklearn)
        logg.info('    using sklearn.manifold.TSNE with a fix by D. DeTomaso')
        X_tsne = tsne.fit_transform(X)
    # update AnnData instance
    adata.obsm['X_tsne'] = X_tsne  # annotate samples with tSNE coordinates
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'X_tsne\', tSNE coordinates (adata.obs)')
    return adata if copy else None

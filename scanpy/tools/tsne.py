# Author: F. Alex Wolf (http://falexwolf.de)
"""tSNE

Notes
-----
This module automatically choose from three t-SNE versions from
- sklearn.manifold.TSNE
- Dmitry Ulyanov (multicore, fastest)
  https://github.com/DmitryUlyanov/Multicore-TSNE
  install via 'pip install psutil cffi', get code from github
"""

import numpy as np
from ..tools.pca import pca
from .. import settings as sett
from .. import logging as logg


def tsne(adata, random_state=0, n_pcs=50, perplexity=30, learning_rate=None,
         recompute_pca=False, use_fast_tsne=True, n_jobs=None, copy=False):
    """tSNE

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with adata.smp['X_pca'], which is
        written when running sc.pca(adata). Is directly used for tSNE.
    random_state : unsigned int or -1, optional (default: 0)
        Change to use different intial states for the optimization, if -1, use
        default behavior of implementation (sklearn uses np.random.seed,
        Multicore-TSNE produces a new plot at every call).
    n_pcs : int, optional (default: 50)
        Number of principal components in preprocessing PCA.
    perplexity : float, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    learning_rate : float, optional (default: 1000)
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    use_fast_tsne : bool, optional (default: True)
        Use the MulticoreTSNE package by D. Ulyanov if available.
    n_jobs : int or None (default: None)
        Use the multicore implementation, if it is installed. Defaults to
        sett.n_jobs.

    Returns
    -------
    Returns or updates adata depending on `copy` with
        "X_tsne", tSNE coordinates of data (adata.smp)

    References
    ----------
    L.J.P. van der Maaten and G.E. Hinton.
    Visualizing High-Dimensional Data Using t-SNE.
    Journal of Machine Learning Research 9(Nov):2579-2605, 2008.

    D. Ulyanov
    Multicore-TSNE
    GitHub (2017)
    """
    logg.info('computing tSNE', r=True)
    adata = adata.copy() if copy else adata
    # preprocessing by PCA
    if 'X_pca' in adata.smp and adata.smp['X_pca'].shape[1] >= n_pcs and not recompute_pca:
        X = adata.smp['X_pca'][:, :n_pcs]
        logg.info('    using X_pca for tSNE')
        logg.info('    using', n_pcs, 'principal components')
    else:
        if n_pcs > 0 and adata.X.shape[1] > n_pcs:
            logg.info('    preprocess using PCA with', n_pcs, 'PCs')
            logg.hint('avoid this by setting n_pcs = 0')
            X = pca(adata.X, random_state=random_state, n_comps=n_pcs)
            adata.smp['X_pca'] = X
            logg.info('    using', n_pcs, 'principal components')
        else:
            X = adata.X
            logg.info('    using data matrix X directly (no PCA)')
    # params for sklearn
    params_sklearn = {'perplexity': perplexity,
                      'random_state': None if random_state == -1 else random_state,
                      'verbose': max(0, sett.verbosity-3),
                      'early_exaggeration': 12,
                      }
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    # deal with different tSNE implementations
    multicore_failed = True
    if n_jobs >= 1 and use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE
            params_sklearn['learning_rate'] = 200 if learning_rate is None else learning_rate
            tsne = TSNE(n_jobs=n_jobs, **params_sklearn)
            logg.info('    using the "MulticoreTSNE" package by Ulyanov (2017)')
            X_tsne = tsne.fit_transform(X.astype(np.float64))
            multicore_failed = False
        except ImportError:
            pass
    if multicore_failed:
        from sklearn.manifold import TSNE
        # unfortunately, we cannot set a minimum number of iterations for barnes-hut
        params_sklearn['learning_rate'] = 1000 if learning_rate is None else learning_rate
        tsne = TSNE(**params_sklearn)
        logg.info('    using sklearn.manifold.TSNE')
        logg.warn('Consider installing the package MulticoreTSNE '
                  ' https://github.com/DmitryUlyanov/Multicore-TSNE '
                  ' Even for `n_jobs=1` this speeds up the computation considerably and will yield better converged results.')
        X_tsne = tsne.fit_transform(X)
    # update AnnData instance
    adata.smp['X_tsne'] = X_tsne
    logg.info('    finished', t=True, end=' ')
    logg.info('and added\n'
              '    "X_tsne", tSNE coordinates (adata.smp)')
    return adata if copy else None

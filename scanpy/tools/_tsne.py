from typing import Optional, Union

from anndata import AnnData

from .._utils import _doc_params, AnyRandom
from ..tools._utils import _choose_representation, doc_use_rep, doc_n_pcs
from .._settings import settings
from .. import logging as logg


@_doc_params(doc_n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def tsne(
    adata: AnnData,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    perplexity: Union[float, int] = 30,
    early_exaggeration: Union[float, int] = 12,
    learning_rate: Union[float, int] = 1000,
    random_state: AnyRandom = 0,
    use_fast_tsne: bool = True,
    n_jobs: Optional[int] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.

    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been
    proposed for visualizating single-cell data by [Amir13]_. Here, by default,
    we use the implementation of *scikit-learn* [Pedregosa11]_. You can achieve
    a huge speedup and better convergence if you install `Multicore-tSNE
    <https://github.com/DmitryUlyanov/Multicore-TSNE>`__ by [Ulyanov16]_, which
    will be automatically detected by Scanpy.

    Parameters
    ----------
    adata
        Annotated data matrix.
    {doc_n_pcs}
    {use_rep}
    perplexity
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    early_exaggeration
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state
        Change this to use different intial states for the optimization.
        If `None`, the initial state is not reproducible.
    use_fast_tsne
        Use the MulticoreTSNE package by D. Ulyanov if it is installed.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_tsne** : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    """
    start = logg.info('computing tSNE')
    adata = adata.copy() if copy else adata
    X = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
    # params for sklearn
    params_sklearn = dict(
        perplexity=perplexity,
        random_state=random_state,
        verbose=settings.verbosity > 3,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
    )
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    # deal with different tSNE implementations
    X_tsne = None
    if n_jobs >= 1 and use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE

            tsne = TSNE(n_jobs=n_jobs, **params_sklearn)
            logg.info("    using the 'MulticoreTSNE' package by Ulyanov (2017)")
            # need to transform to float64 for MulticoreTSNE...
            X_tsne = tsne.fit_transform(X.astype('float64'))
        except ImportError:
            logg.warning(
                'Consider installing the package MulticoreTSNE '
                '(https://github.com/DmitryUlyanov/Multicore-TSNE). '
                'Even for n_jobs=1 this speeds up the computation considerably '
                'and might yield better converged results.'
            )
    if X_tsne is None:
        from sklearn.manifold import TSNE
        from . import _tsne_fix  # fix by D. DeTomaso for sklearn < 0.19

        # unfortunately, sklearn does not allow to set a minimum number
        # of iterations for barnes-hut tSNE
        tsne = TSNE(**params_sklearn)
        logg.info('    using sklearn.manifold.TSNE with a fix by D. DeTomaso')
        X_tsne = tsne.fit_transform(X)
    # update AnnData instance
    adata.obsm['X_tsne'] = X_tsne  # annotate samples with tSNE coordinates
    logg.info(
        '    finished',
        time=start,
        deep="added\n    'X_tsne', tSNE coordinates (adata.obsm)",
    )
    return adata if copy else None

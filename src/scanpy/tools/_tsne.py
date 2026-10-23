from __future__ import annotations

from typing import TYPE_CHECKING

from scverse_misc import Deprecation, deprecated_arg

from .. import logging as logg
from .._compat import warn
from .._docs import DEPR_KEY_ADDED, doc_out_key, doc_rng
from .._keys import _embedding_keys
from .._settings import Default, settings
from .._utils import _doc_params, raise_not_implemented_error_if_backed_type
from .._utils.random import _accepts_legacy_random_state, _legacy_random_state
from ..get._out import _out_key, deprecated_copy
from ..get._uns import _use_to_uns
from ..get._use import _resolve_obs
from ..neighbors._doc import doc_n_pcs, doc_use_rep
from ._utils import _choose_representation_compat

if TYPE_CHECKING:
    import numpy as np
    from anndata import AnnData
    from anndata.acc import MultiAcc

    from .._utils.random import RNGLike, SeedLike
    from ..get._out import Prefix
    from ..get._types import RepAcc


@_accepts_legacy_random_state(0)
@_doc_params(
    doc_n_pcs=doc_n_pcs,
    use_rep=doc_use_rep,
    rng=doc_rng,
    out=doc_out_key(
        "the embedding",
        example="A.obsm['tsne']",
        prefixed="`.obsm['p_tsne']` and `.uns['p_tsne']`",
        key_added="`key_added='k'` is `out=A.obsm['k']`.",
        default=(
            "`.obsm['X_tsne']` and `.uns['tsne']` "
            "(`.obsm['tsne']` if :attr:`scanpy.settings.preset` is "
            ":attr:`~scanpy.Preset.ScanpyV2Preview`)"
        ),
    ),
)
@deprecated_arg("use_rep", Deprecation("1.13.0", "Use `use` instead."))
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_copy
def tsne(  # noqa: PLR0913
    adata: AnnData,
    n_pcs: int | None = None,
    *,
    n_components: int = 2,
    use: RepAcc | str | None = None,
    perplexity: float = 30,
    metric: str = "euclidean",
    early_exaggeration: float = 12,
    learning_rate: float = 1000,
    rng: SeedLike | RNGLike | None = None,
    use_fast_tsne: bool = False,
    n_jobs: int | None = None,
    out: MultiAcc | str | Prefix | bool = False,
    # deprecated
    key_added: str | None = None,
    copy: bool = False,
    use_rep: str | None = None,
) -> AnnData | np.ndarray | None:
    r"""t-SNE :cite:p:`vanDerMaaten2008,Amir2013,Pedregosa2011`.

    t-distributed stochastic neighborhood embedding (tSNE, :cite:t:`vanDerMaaten2008`) was
    proposed for visualizating single-cell data by :cite:t:`Amir2013`. Here, by default,
    we use the implementation of *scikit-learn* :cite:p:`Pedregosa2011`. You can achieve
    a huge speedup and better convergence if you install Multicore-tSNE_
    by :cite:t:`Ulyanov2016`, which will be automatically detected by Scanpy.

    .. _multicore-tsne: https://github.com/DmitryUlyanov/Multicore-TSNE

    .. array-support:: tl.tsne

    Parameters
    ----------
    adata
        Annotated data matrix.
    {doc_n_pcs}
    {use_rep}
    n_components
        The number of dimensions of the embedding.
    perplexity
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    metric
        Distance metric calculate neighbors on.
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
    {rng}
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy.settings.n_jobs`.
    {out}
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Returns the embedding if `out=True`, else `None`. Sets the following fields:

    `adata.obsm['X_tsne' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        tSNE coordinates of data.
    `adata.uns['tsne' | key_added]` : :class:`dict`
        tSNE parameters.

    """
    start = logg.info("computing tSNE")
    if use is not None:
        use = _resolve_obs(use, rep=True)
    key = _out_key(out, key_added, kind="obsm", base="tsne", default=Default())
    x = _choose_representation_compat(adata, use=use, use_rep=use_rep, n_pcs=n_pcs)
    raise_not_implemented_error_if_backed_type(x, "tsne")
    # params for sklearn
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    params_sklearn = dict(
        perplexity=perplexity,
        random_state=_legacy_random_state(rng),
        verbose=settings.verbosity > 3,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_jobs=n_jobs,
        metric=metric,
        n_components=n_components,
    )

    if use_fast_tsne:
        msg = (
            "Argument `use_fast_tsne` is deprecated, and support for MulticoreTSNE "
            "will be dropped in a future version of scanpy."
        )
        warn(msg, FutureWarning)

    # deal with different tSNE implementations
    if use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE  # noqa: N814
        except ImportError:
            use_fast_tsne = False
            msg = "Could not import 'MulticoreTSNE'. Falling back to scikit-learn."
            warn(msg, ImportWarning)
        else:
            tsne = TSNE(**params_sklearn)
            logg.info("    using the 'MulticoreTSNE' package by Ulyanov (2017)")
            # need to transform to float64 for MulticoreTSNE...
            x_tsne = tsne.fit_transform(x.astype("float64"))
    if use_fast_tsne is False:  # In case MultiCore failed to import
        from sklearn.manifold import TSNE

        # unfortunately, sklearn does not allow to set a minimum number
        # of iterations for barnes-hut tSNE
        tsne = TSNE(**params_sklearn)
        logg.info("    using sklearn.manifold.TSNE")
        x_tsne = tsne.fit_transform(x)

    # update AnnData instance
    params = dict(
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_jobs=n_jobs,
        metric=metric,
        **_use_to_uns(use, use_rep, rep=True),
        n_components=n_components,
    )
    if key is True:
        logg.info("    finished", time=start)
        return x_tsne
    keys = _embedding_keys("tsne", key)
    adata.obsm[keys.obsm] = x_tsne  # annotate samples with tSNE coordinates
    adata.uns[keys.uns] = dict(
        params={k: v for k, v in params.items() if v is not None}
    )

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added\n"
            f"    {keys.obsm!r}, tSNE coordinates (adata.obsm)\n"
            f"    {keys.uns!r}, tSNE parameters (adata.uns)"
        ),
    )

    return None

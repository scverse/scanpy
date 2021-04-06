"""\
Embed high-dimensional data via Shannon Components Analysis (SCA)
"""

import typing
from anndata import AnnData


def sca(
    adata: AnnData,
    n_comps: int = 50,
    iters: int = 1,
    nbhd_size: int = 15,
    nbhds=None,
    rep=None,
    n_pcs: int = 50,
    metric: str = 'cosine',
    max_bins: float = float('inf'),
    scaled_bins: bool = False,
    fast_version: bool = True,
    keep_scores: bool = False,
    keep_loadings: bool = True,
    keep_all_iters: bool = False,
    verbose: bool = False,
    multi_correct: bool = True,
    n_tests: int = 50,
    copy: bool = False,
    key_added: str = 'sca',
    layer: typing.Optional[str] = None,
    **kwargs,
):
    """\
    Reduce the data using Shannon Component Analysis (SCA) [DeMeo21]_.

    SCA compares local and global gene expression measurements to compute Shannon
    information scores for each measurement, and outputs a linear dimensionality
    reduction with high Shannon information. This function should be run after
    normalization but before computing the neighbor graph; it can be used in place of
    PCA. Since information scores are computed using binarized counts, any
    pre-normalizations must preserve zero counts, i.e. zero in adata.X represents
    zero recorded (or imputed) transcripts.

    This uses the implementation of `shannonca
    <https://github.com/bendemeo/shannonca>`__ [DeMeo21]_.

    Parameters
    ----------
    adata
        The annotated data matrix.
    n_comps
        Target dimensionality, i.e. number of top Shannon components.
    iters
        Number of SCA iterations to run. More iterations often increases
        the discriminative power of the representation. Typically stable after
        5 or so iterations.
    nbhd_size
        Size of neighborhoods used to compute local gene expression probabilities.
    nbhds
        (Optional) Array of size `(adata.shape[0], nbhd_size)` specifying the indices
        of the nearest neighbors of each cell. If not specified, computes
        neighborhoods at runtime.
    rep
        (Optional) Initial representation in which local neighborhoods are computed.
        Ignored if nbhds is supplied. If not specified, computes a PCA representation
    n_pcs
        (Optional) Number of principal components used when computing local
        neighborhoods in the first iteration. Ignored if `rep` or `nbhds` is supplied.
    metric
        Metric used to compute nearest neighbors. cosine distance performs well on
        most data.
    max_bins
        Maximum number of gene probability bins used in the computation. Bins
        are used to approximate global gene probabilities. Lower numbers reduce
        the memory footprint at the expense of accuracy. If infinite (default),
        gene probabilities are not binned and scores are computed exactly.
    scaled_bins
        Whether or not to scale bins so that approximations to gene probabilities
        are correct to within a constant factor. Defaults to
        evenly-spaced bins. Ignored if `max_bins` is infinite.
    fast_version
        If `True` (default), Shannon information scores are computed via a vectorized
        matrix multiplication; however, this can be memory-intensive. Setting
        `fast_version=False` instead iterates over cells and computes the scores
        individually, which takes longer but is less memory-intensive.
    keep_scores
        If True, add the information scores for each gene in cells as a layer
        to adata, keyed by `key_added+'_score'`. Scores are stored as a sparse matrix;
        however, score matrices are typically less sparse than `adata.X`, so this may
        increase the memory footprint of adata significantly.
    keep_loadings:
        If `True` (default), store the loadings of each gene in each shannon component
        under adata.varm[key_added+'_loadings']. Allows interpretation of components.
    keep_all_iters:
        If `True`, store intermediate representations after each iteration, under
        adata.obsm[key_added+'_i'] for i=1,2,...,iters.
    verbose:
        If `True`, print progress as SCA runs.
    multi_correct:
        Whether or not to correct information scores for multiple testing. If `True`,
        the Shannon information computations in each cell are corrected via a
        family-wise error rate correction.
    n_tests:
        Effective number of tests per cell used in FWER correction.
    copy:
        If `True`, returns a copy of `adata`. If `False`, updates `adata` in place.
    key_added:
        Key under which SCA results are stored in `adata`.
    layer:
        Layer on which SCA runs. If None, will run on `adata.X`.
    kwargs
        Any additional arguments will be passed to
        ``shannonca.dimred.reduce_scanpy``


    Returns
    -------
    If copy=True, returns a copy of adata with the SCA reduction in
    adata.obsm[key_added], or adata.obsm[key_added+'_'+iters] if
    keep_all_iters=True. Loadings, and scores are optionally added
    to adata.varm[key_added+'_loadings'] and adata.layers[key_added+'_scores')
    respectively.
    """

    try:
        from shannonca.dimred import reduce_scanpy
    except ImportError:
        raise ImportError('\n please install shannonca \n\n\tpip install shannonca')

    adata = adata.copy() if copy else adata

    reduce_scanpy(
        adata,
        n_comps=n_comps,
        iters=iters,
        nbhd_size=nbhd_size,
        nbhds=nbhds,
        rep=rep,
        n_pcs=n_pcs,
        metric=metric,
        max_bins=max_bins,
        scaled_bins=scaled_bins,
        fast_version=fast_version,
        keep_scores=keep_scores,
        keep_loadings=keep_loadings,
        keep_all_iters=keep_all_iters,
        verbose=verbose,
        multi_correct=multi_correct,
        n_tests=n_tests,
        key_added=key_added,
        layer=layer,
        **kwargs,
    )

    if copy:
        return adata.copy()

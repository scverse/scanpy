"""Calculate scores based on relative expression change of maker pairs
"""

def sandbag(
        adata,
        annotation,
        gene_names,
        sample_names,
        fraction=0.65,
        filter_genes=None,
        filter_samples=None):
    """Generate pairs of genes [Scialdone15]_ [Fechtner18]_.

    Calculates the pairs of genes serving as marker pairs for each phase,
    based on a matrix of gene counts and an annotation of known phases.

    This reproduces the approach of [Scialdone15]_ in the implementation of
    [Fechtner18]_.

    More information and bug reports `here
    <https://github.com/rfechtner/pypairs>`__.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    categories : `dict`
        Dictionary of lists, i.e. {phase: [sample, ...]},
        containing annotation of samples to their phase
    gene_names: `list`
        List of genes.
    sample_names: `list`
        List of samples.
    fraction : `float`, optional (default: 0.5)
        Fraction to be used as threshold.
    filter_genes : `list` or `None`, optional (default: `None`)
        Genes for sampling the reference set. Default is all genes.
    filter_samples : `list` or `None`, optional (default: `None`)
        Cells for sampling the reference set. Default is all samples.

    Returns
    -------
    `dict` of `list` of `tuple`, i.e.
    {phase: [(Gene1, Gene2), ...]},
    containing marker pairs per phase
    """
    try:
        from pypairs import __version__ as pypairsversion
        from distutils.version import LooseVersion

        if LooseVersion(pypairsversion) < LooseVersion("v3.0.9"):
            raise ImportError('Please only use `pypairs` >= v3.0.9 ')
    except ImportError:
        raise ImportError('You need to install the package `pypairs`.')


    from pypairs.pairs import sandbag
    from . import settings
    from pypairs import settings as pp_settings

    pp_settings.verbosity = settings.verbosity
    pp_settings.n_jobs = settings.n_jobs
    pp_settings.writedir = settings.writedir
    pp_settings.cachedir = settings.cachedir
    pp_settings.logfile = settings.logfile

    return sandbag(
        data = adata,
        annotation = annotation,
        gene_names = gene_names,
        sample_names = sample_names,
        fraction = fraction,
        filter_genes = filter_genes,
        filter_samples = filter_samples
    )


def cyclone(
        adata,
        marker_pairs,
        gene_names,
        sample_names,
        iterations=1000,
        min_iter=100,
        min_pairs=50):
    """Assigns scores and predicted class to observations [Scialdone15]_ [Fechtner18]_.

    Calculates scores for each observation and each phase and assigns prediction
    based on marker pairs indentified by sandbag.

    This reproduces the approach of [Scialdone15]_ in the implementation of
    [Fechtner18]_.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    marker_pairs : `dict`
        Dictionary of marker pairs. See :func:`~scanpy.api.sandbag` output.
    gene_names: `list`
        List of genes.
    sample_names: `list`
        List of samples.
    iterations : `int`, optional (default: 1000)
        An integer scalar specifying the number of
        iterations for random sampling to obtain a cycle score.
    min_iter : `int`, optional (default: 100)
        An integer scalar specifying the minimum number of iterations
        for score estimation
    min_pairs : `int`, optional (default: 50)
        An integer scalar specifying the minimum number of iterations
        for score estimation

    Returns
    -------
    A :class:`~pandas.DataFrame` with samples as index and categories as columns with scores for each category for each
    sample and a additional column with the name of the max scoring category for each sample.

    If marker pairs contain only the cell cycle categories G1, S and G2M an additional column
    ``pypairs_cc_prediction`` will be added. Where category S is assigned to samples where G1 and G2M score are
    below 0.5.
    """
    try:
        from pypairs import __version__ as pypairsversion
        from distutils.version import LooseVersion

        if LooseVersion(pypairsversion) < LooseVersion("v3.0.9"):
            raise ImportError('Please only use `pypairs` >= v3.0.9 ')
    except ImportError:
        raise ImportError('You need to install the package `pypairs`.')


    from pypairs.pairs import cyclone
    from . import settings
    from pypairs import settings as pp_settings

    pp_settings.verbosity = settings.verbosity
    pp_settings.n_jobs = settings.n_jobs
    pp_settings.writedir = settings.writedir
    pp_settings.cachedir = settings.cachedir
    pp_settings.logfile = settings.logfile

    return cyclone(
        data = adata,
        marker_pairs = marker_pairs,
        gene_names = gene_names,
        sample_names = sample_names,
        iterations = iterations,
        min_iter = min_iter,
        min_pairs = min_pairs
    )

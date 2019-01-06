"""Calculate scores based on relative expression change of maker pairs
"""

import pandas as pd


def sandbag(
        adata,
        categories,
        gene_names,
        sample_names,
        fraction=0.5,
        subset_genes=None,
        subset_samples=None,
        filter_genes_dispersion=False,
        always_keep_genes=None,
        min_diff=0,
        rm_zeros=True,
        n_jobs=1):
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
    subset_genes : `list` or `None`, optional (default: `None`)
        Genes for sampling the reference set. Default is all genes.
    subset_samples : `list` or `None`, optional (default: `None`)
        Cells for sampling the reference set. Default is all samples.
    filter_genes_dispersion: `bool`, optional (default: False)
        Filtered genes based on dispersion
    always_keep_genes: `list`, optional (default: None)
        List of genes to always keep
    min_diff : `int`, optional (default: 0)
        Minimum different
    rm_zeros: `bool`, optional (default: True)
        Remove zeros
    n_jobs : `int`, optional (default: 1)
        Number of concurrent n_jobs to be used. 0 = use all available cores.

    Returns
    -------
    `dict` of `list` of `tuple`, i.e.
    {phase: [(Gene1, Gene2), ...]},
    containing marker pairs per phase
    """
    try:
        from pypairs import sandbag
    except ImportError:
        raise ImportError('You need to install the package `pypairs`.')

    x = pd.DataFrame(adata.X)

    return sandbag.sandbag(
        matrix=x,
        categories=categories,
        gene_names=gene_names,
        sample_names=sample_names,
        fraction=fraction,
        subset_genes=subset_genes,
        subset_samples=subset_samples,
        rm_zeros=rm_zeros,
        filter_genes_dispersion=filter_genes_dispersion,
        always_keep_genes=always_keep_genes,
        min_diff=min_diff,
        processes=n_jobs
    )


def cyclone(
        adata,
        marker_pairs,
        gene_names,
        sample_names,
        subset_genes=None,
        subset_samples=None,
        iterations=1000,
        min_iter=100,
        min_pairs=50,
        rm_zeros=True,
        n_jobs=1):
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
    subset_genes : `list` or `None`, optional (default: `None`)
        Genes for sampling the reference set. Default is all genes.
    subset_samples : `list` or `None`, optional (default: `None`)
        Cells for sampling the reference set. Default is all samples.
    iterations : `int`, optional (default: 1000)
        An integer scalar specifying the number of
        iterations for random sampling to obtain a cycle score.
    min_iter : `int`, optional (default: 100)
        An integer scalar specifying the minimum number of iterations
        for score estimation
    min_pairs : `int`, optional (default: 50)
        An integer scalar specifying the minimum number of iterations
        for score estimation
    rm_zeros: `bool`, optional (default: True)
        Remove zeros
    n_jobs : `int`, optional (default: 1)
        Number of concurrent n_jobs to be used. 0 = use all available cores.

    Returns
    -------
    `dict` of `list`
    {
        "prediction": The predicted classes based on scores
        "prediction_normalized": The predicted classes based on normalized scores
        "scores": Prediction scores
        "normalized": Normalized prediction scores
    }
    """
    try:
        from pypairs import cyclone
    except ImportError:
        raise ImportError('You need to install the package `pypairs`.')

    x = pd.DataFrame(adata.X)

    return cyclone.cyclone(
        matrix=x,
        marker_pairs=marker_pairs,
        gene_names=gene_names,
        sample_names=sample_names,
        subset_genes=subset_genes,
        subset_samples=subset_samples,
        iterations=iterations,
        min_iter=min_iter,
        min_pairs=min_pairs,
        rm_zeros=rm_zeros,
        processes=n_jobs
    )

"""Calculate scores based on relative expression change of maker pairs
"""

import pypairs
import pandas as pd

def sandbag(
        adata,
        phases,
        fraction=0.5,
        subset_genes=None,
        subset_samples=None,
        processes=1):
    """Generate pairs of genes [Scialdone15]_.

    Calculates the pairs of genes serving as marker pairs for each phase,
    based on a matrix of gene counts and an annotation of known phases.

    This reproduces the approach in Pairs [Scialdone15]_ and has been
    implemented for Scanpy by Ron Fechtner.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    phases : `dict`
        Dictionary of lists, i.e. {phase: [sample, ...]},
        containing annotation of samples to their phase
    fraction : `float`, optional (default: 0.5)
        Fraction to be used as threshold.
    subset_genes : `list` or `None`, optional (default: `None`)
        Genes for sampling the reference set. Default is all genes.
    subset_samples : `list` or `None`, optional (default: `None`)
        Cells for sampling the reference set. Default is all samples.
    processes : `int`, optional (default: 1)
        Number of concurrent processes to be used. 0 = use all available cores.

    Returns
    -------
    `dict` of `list` of `tuple`, i.e.
    {phase: [(Gene1, Gene2), ...]},
    containing marker pairs per phase

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`_.
    """

    x = pd.DataFrame(adata.X)

    return pypairs.sandbag(
        x=x,
        phases=phases,
        subset_genes=subset_genes,
        subset_samples=subset_samples,
        processes=processes
    )


def cyclone(
        adata,
        marker_pairs,
        subset_genes=None,
        subset_samples=None,
        iterations=1000,
        min_iter=100,
        min_pairs=50,
        processes=1):
    """Assigns scores and predicted class to samples [Scialdone15]_.

    Calculates scores for each sample and each phase and assigns prediction
    based on marker pairs indentified by sandbag

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    marker_pairs : `dict`
        Dictionary of marker pairs. See :func:`~scanpy.api.sandbag` output.
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
    processes : `int`, optional (default: 1)
        Number of concurrent processes to be used. 0 = use all available cores.

    Returns
    -------
    Dictionary of List containing
    {
        "prediction": The predicted classes based on scores
        "prediction_normalized": The predicted classes based on normalized scores
        "scores": Prediction scores
        "normalized": Normalized prediction scores
    }

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`_.
    """
    x = pd.DataFrame(adata.X)

    return pypairs.cyclone(
        x=x,
        marker_pairs=marker_pairs,
        subset_genes=subset_genes,
        subset_samples=subset_samples,
        iterations=iterations,
        min_iter=min_iter,
        min_pairs=min_pairs,
        processes=processes
    )


def load_oscope_marker(fraction=0.6, cc_only=True, processes=1):
    """Load default marker pairs for :func:`~scanpy.api.cyclone`

    Generated marker pairs based on the `Oscope dataset <http://dx.doi.org/10.1038/nmeth.3549>`_

    Parameters
    ----------
    fraction : `float`, optional (default: 0.6)
        Fraction to be used as threshold.
    cc_only : `bool`, optional (default: True)
        Use cell cycle related genes only
    processes : `int`, optional (default: 1)
        Number of concurrent processes to be used. 0 = use all available cores.

    Returns
    -------
    Dictionary of List of Tuples with marker pairs for each phase G1, S, G2M.
    """
    gencounts_oscope = pd.read_csv(
        Path("../datasets/GSE64016_H1andFUCCI_normalized_EC_human.csv")
    )

    gencounts_oscope.set_index("Unnamed: 0", inplace=True)

    gencounts_oscope_sorted = gencounts_oscope.iloc[:,
                              [gencounts_oscope.columns.get_loc(c)
                               for c in gencounts_oscope.columns
                               if "G1_" in c or "G2_" in c or "S_" in c]]

    is_G1 = [
        gencounts_oscope_sorted.columns.get_loc(c)
        for c in gencounts_oscope_sorted.columns if "G1_" in c
    ]
    is_S = [
        gencounts_oscope_sorted.columns.get_loc(c)
        for c in gencounts_oscope_sorted.columns if "S_" in c]
    is_G2M = [
        gencounts_oscope_sorted.columns.get_loc(c)
        for c in gencounts_oscope_sorted.columns if "G2_" in c
    ]

    annotation = {
        "G1": list(is_G1),
        "S": list(is_S),
        "G2M": list(is_G2M)
    }

    if cc_only:
        go_0007049 = [line.replace("\n", "").replace("\r", "") for line in
                      open("../datasets/go_0007049_homoSapiens.csv", "r")]
        cycle_base = [
            line.split("\t")[0] for i, line
            in enumerate(open("./datasets/cyclebase_top1000_genes.tsv", "r"))
            if 0 < i
        ]
        cycle_genes = numpy.unique(numpy.concatenate((go_0007049, cycle_base), 0))

        return pairs.sandbag(
            gencounts_oscope_sorted, phases=annotation,
            subset_genes=list(cycle_genes), fraction=fraction,
            processes=processes
        )
    else:
        return pairs.sandbag(
            gencounts_oscope_sorted, phases=annotation,
            fraction=fraction, processes=processes
        )

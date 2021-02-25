"""\
Run Diffusion maps using the adaptive anisotropic kernel
"""
from typing import Optional, List

import pandas as pd
from anndata import AnnData

from ... import logging as logg


def palantir(
    adata: AnnData,
    n_components: int = 10,
    knn: int = 30,
    alpha: float = 0,
    use_adjacency_matrix: bool = False,
    distances_key: Optional[str] = None,
    n_eigs: int = None,
    impute_data: bool = True,
    n_steps: int = 3,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Run Diffusion maps using the adaptive anisotropic kernel [Setty18]_.

    Palantir is an algorithm to align cells along differentiation trajectories.
    Palantir models differentiation as a stochastic process where stem cells
    differentiate to terminally differentiated cells by a series of steps through
    a low dimensional phenotypic manifold. Palantir effectively captures the
    continuity in cell states and the stochasticity in cell fate determination.
    Palantir has been designed to work with multidimensional single cell data
    from diverse technologies such as Mass cytometry and single cell RNA-seq.

    .. note::
       More information and bug reports `here <https://github.com/dpeerlab/Palantir>`__.

    Parameters
    ----------
    adata
        An AnnData object.
    n_components
        Number of diffusion components.
    knn
        Number of nearest neighbors for graph construction.
    alpha
        Normalization parameter for the diffusion operator.
    use_adjacency_matrix
        Use adaptive anisotropic adjacency matrix, instead of PCA projections
        (default) to compute diffusion components.
    distances_key
        With `use_adjacency_matrix=True`, use the indicated distances key for `.obsp`.
        If `None`, `'distances'`.
    n_eigs
        Number of eigen vectors to use. If `None` specified, the number of eigen
        vectors will be determined using eigen gap. Passed to
        `palantir.utils.determine_multiscale_space`.
    impute_data
        Impute data using MAGIC.
    n_steps
        Number of steps in the diffusion operator. Passed to
        `palantir.utils.run_magic_imputation`.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields:

    **Diffusion maps**,
        used for magic imputation, and to generate multi-scale data matrix,

        - X_palantir_diff_comp - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.obsm`, dtype `float`)
            Array of Diffusion components.
        - palantir_EigenValues - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.uns`, dtype `float`)
            Array of corresponding eigen values.
        - palantir_diff_op - :class:`~scipy.sparse.spmatrix` (:attr:`~anndata.AnnData.obsp`, dtype `float`)
            The diffusion operator matrix.

    **Multi scale space results**,
        used to build tsne on diffusion components, and to compute branch probabilities
        and waypoints,

        - X_palantir_multiscale - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.obsm`, dtype `float`)
            Multi scale data matrix.

    **MAGIC imputation**,
        used for plotting gene expression on tsne, and gene expression trends,

        - palantir_imp - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.layers`, dtype `float`)
            Imputed data matrix (MAGIC imputation).

    Example
    -------
    >>> import scanpy.external as sce
    >>> import scanpy as sc

    A sample data is available `here <https://github.com/dpeerlab/Palantir/tree/master/data>`_.

    **Load sample data**

    >>> adata = sc.read_csv(filename="Palantir/data/marrow_sample_scseq_counts.csv.gz")

    *Cleanup and normalize*

    >>> sc.pp.filter_cells(adata, min_counts=1000)
    >>> sc.pp.filter_genes(adata, min_counts=10)
    >>> sc.pp.normalize_per_cell(adata)
    >>> sc.pp.log1p(adata)

    **Data preprocessing**

    Palantir builds diffusion maps using one of two optional inputs:

    *Principal component analysis*

    >>> sc.tl.pca(adata, n_comps=300)

    or,

    *Nearist neighbors graph*

    >>> sc.pp.neighbors(adata, knn=30)

    *Diffusion maps*

    Palantir determines the diffusion maps of the data as an estimate of the low
    dimensional phenotypic manifold of the data.

    >>> sce.tl.palantir(adata, n_components=5, knn=30)

    if pre-computed distances are to be used,

    >>> sce.tl.palantir(
    ...     adata,
    ...     n_components=5,
    ...     knn=30,
    ...     use_adjacency_matrix=True,
    ...     distances_key="distances",
    ... )

    **Visualizing Palantir results**

    *tSNE visualization*

    important for Palantir!

    Palantir constructs the tSNE map in the embedded space since these maps better
    represent the differentiation trajectories.

    >>> sc.tl.tsne(adata, n_pcs=2, use_rep='X_palantir_multiscale', perplexity=150)

    *tsne by cell size*

    >>> sc.pl.tsne(adata, color="n_counts")

    *Imputed gene expression visualized on tSNE maps*

    >>> sc.pl.tsne(
    ...     adata,
    ...     gene_symbols=['CD34', 'MPO', 'GATA1', 'IRF8'],
    ...     layer='palantir_imp',
    ...     color=['CD34', 'MPO', 'GATA1', 'IRF8']
    ... )

    **Running Palantir**

    Palantir can be run by specifying an approximate early cell. While Palantir
    automatically determines the terminal states, they can also be specified using the
    `termine_states` parameter.

    >>> start_cell = 'Run5_164698952452459'
    >>> pr_res = sce.tl.palantir_results(
    ...     adata,
    ...     early_cell=start_cell,
    ...     ms_data='X_palantir_multiscale',
    ...     num_waypoints=500,
    ... )

    .. note::
       A `start_cell` must be defined for every data set. The start cell for
       this dataset was chosen based on high expression of CD34.

    At this point the returned Palantir object `pr_res` can be used for all downstream
    analysis and plotting. Please consult this notebook
    `Palantir_sample_notebook.ipynb
    <https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst other
    things.
    """

    _check_import()
    from palantir.utils import (
        run_diffusion_maps,
        determine_multiscale_space,
        run_magic_imputation,
    )

    adata = adata.copy() if copy else adata

    logg.info('Palantir Diffusion Maps in progress ...')

    if use_adjacency_matrix:
        df = adata.obsp[distances_key] if distances_key else adata.obsp["distances"]
    else:
        df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)

    # Diffusion maps
    dm_res = run_diffusion_maps(
        data_df=df,
        n_components=n_components,
        knn=knn,
        alpha=alpha,
    )
    # Determine the multi scale space of the data
    ms_data = determine_multiscale_space(dm_res=dm_res, n_eigs=n_eigs)

    # MAGIC imputation
    if impute_data:
        imp_df = run_magic_imputation(
            data=adata.to_df(), dm_res=dm_res, n_steps=n_steps
        )
        adata.layers['palantir_imp'] = imp_df

    (
        adata.obsm['X_palantir_diff_comp'],
        adata.uns['palantir_EigenValues'],
        adata.obsp['palantir_diff_op'],
        adata.obsm['X_palantir_multiscale'],
    ) = (
        dm_res['EigenVectors'].to_numpy(),
        dm_res['EigenValues'].to_numpy(),
        dm_res['T'],
        ms_data.to_numpy(),
    )

    return adata if copy else None


def palantir_results(
    adata: AnnData,
    early_cell: str,
    ms_data: str = 'X_palantir_multiscale',
    terminal_states: List = None,
    knn: int = 30,
    num_waypoints: int = 1200,
    n_jobs: int = -1,
    scale_components: bool = True,
    use_early_cell_as_start: bool = False,
    max_iterations: int = 25,
) -> Optional[AnnData]:
    """\
    **Running Palantir**

    A convenience function that wraps `palantir.core.run_palantir` to compute branch
    probabilities and waypoints.

    Parameters
    ----------
    adata
        An AnnData object.
    early_cell
        Start cell for pseudotime construction.
    ms_data
        Palantir multi scale data matrix,
    terminal_states
        List of user defined terminal states
    knn
        Number of nearest neighbors for graph construction.
    num_waypoints
        Number of waypoints to sample.
    n_jobs
        Number of jobs for parallel processing.
    scale_components
        Transform features by scaling each feature to a given range. Consult the
        documentation for `sklearn.preprocessing.minmax_scale`.
    use_early_cell_as_start
        Use `early_cell` as `start_cell`, instead of determining it from the boundary
        cells closest to the defined `early_cell`.
    max_iterations
        Maximum number of iterations for pseudotime convergence.

    Returns
    -------
    PResults
        PResults object with pseudotime, entropy, branch probabilities and waypoints.
    """
    logg.info('Palantir computing waypoints..')

    _check_import()
    from palantir.core import run_palantir

    ms_data = pd.DataFrame(adata.obsm[ms_data], index=adata.obs_names)
    pr_res = run_palantir(
        ms_data=ms_data,
        early_cell=early_cell,
        terminal_states=terminal_states,
        knn=knn,
        num_waypoints=num_waypoints,
        n_jobs=n_jobs,
        scale_components=scale_components,
        use_early_cell_as_start=use_early_cell_as_start,
        max_iterations=max_iterations,
    )

    return pr_res


def _check_import():
    try:
        import palantir
    except ImportError:
        raise ImportError('\nplease install palantir:\n\tpip install palantir')

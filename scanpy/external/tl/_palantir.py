"""\
Run Diffusion maps using the adaptive anisotropic kernel
"""
from typing import Optional

from anndata import AnnData

from ... import logging as logg


def palantir(
    adata: AnnData,
    normalize: bool = False,
    log_transform: bool = False,
    filter_low: bool = False,
    inplace: bool = True,
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
        An AnnData object, or Dataframe of cells X genes.
    normalize
        property setter passed to palantir to normalize using palantir method
        `palantir.preprocess.normalize_counts`.
    log_transform
        property setter passed to palantir. Some datasets show better signal in the log
        scale. Applied using `palantir.preprocess.log_transform`
    filter_low
        property setter passed to palantir to remove low molecule count cells and low detection genes
    inplace
        Set fields in `adata` or return a copy?

    Returns
    -------
    `.uns['palantir_norm_data']`
        A `data_df` copy of adata if normalized
    `pca_results`
        PCA projections and explained variance ratio of adata:
        - `.uns['palantir_pca_results']['pca_projections']`
        - `.uns['palantir_pca_results']['variance_ratio']`
    `dm_res`
        Diffusion components, corresponding eigen values and diffusion operator:
        - `.uns['palantir_diff_maps']['EigenVectors']`
        - `.uns['palantir_diff_maps']['EigenValues']`
        - `.uns['palantir_diff_maps']['T']`
    `.uns['palantir_ms_data']`
        The `ms_data` – Multi scale data matrix
    `.uns['palantir_tsne']` : `tsne`
        tSNE on diffusion maps
    `.uns['palantir_imp_df']` : `imp_df`
        Imputed data matrix (MAGIC imputation)

    Example
    -------
    >>> import scanpy.external as sce
    >>> import scanpy as sc

    A sample data is available `here <https://github.com/dpeerlab/Palantir/tree/master/data>`_.

    To view the plots, it is recommended to run Jupyter notebook

    **Load sample data**

    >>> adata = sc.read_csv(filename="Palantir/data/marrow_sample_scseq_counts.csv.gz")

    **Pre-processing**

    The provided adata will be used as input to the embedded `palantir` methods:

    >>> d = sce.tl.palantir(adata)

    At this point, a new class object, `d`, will be instantiated. If the data
    needs pre-processing – filtering low genes/cells counts, or normalization,
    or log transformation, set the `filter_low`, `normalize`, or `log_transform`
    to `True`:

    >>> d.filter_low = True
    >>> d.normalize = True
    >>> d.log_transform = True

    The created object `d.palantir` can be used to override the default
    parameters used for pre-processing.

    Follow the next step to pass the data to palantir methods, to generate the
    return objects listed above.

    **Run Palantir**

    >>> d.process()

    By calling this method `palantir` will run and generate the various outputs.
    The generated objects will be pushed to `adata` and stored for further use.
    Once instantiated, *Principal component analysis*, *Diffusion maps*,
    *tSNE on Diffusion maps*, and *MAGIC imputation* data objects will be created
    using the `palantir` default parameters.

    If running `palantir` using default parameters is not satisfactory,
    `d.palantir` methods can be used to override and substitute the individual
    outputs already embedded into `adata`.

    **Plotting**

    *tSNE visualization*

    >>> fig, ax = d.palantir.plot.plot_tsne(d.tsne)
    >>> fig, ax = d.palantir.plot.plot_tsne_by_cell_sizes(d.data_df, d.tsne)

    *Gene expression can be visualized on tSNE maps*

    >>> d.palantir.plot.plot_gene_expression(d.imp_df, d.tsne, ['CD34', 'MPO', 'GATA1', 'IRF8'])

    *Diffusion maps*

    >>> d.palantir.plot.plot_diffusion_components(d.tsne, d.dm_res)

    **Visualizing Palantir results**

    Palantir can be run by specifying an approximate early cell. While Palantir
    automatically determines the terminal states, they can also be specified using the
    `termine_states` parameter.

    >>> start_cell = 'Run5_164698952452459'
    >>> pr_res = d.palantir.core.run_palantir(d.ms_data, start_cell, num_waypoints=500)
    >>> palantir.plot.plot_palantir_results(pr_res, d.tsne)

    .. note::
       A `start_cell` must be defined for every data set. The start cell for
       this dataset was chosen based on high expression of CD34.

    For further demonstration of palantir visualizations please follow this notebook
    `Palantir_sample_notebook.ipynb <https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst other things.
    """

    logg.info('Palantir diffusion maps')

    try:
        import palantir
    except ImportError:
        raise ImportError(
            '\nplease install palantir: \n\n'
            '\tgit clone git://github.com/dpeerlab/Palantir.git\n'
            '\tcd Palantir\n'
            '\tsudo -H pip3 install .'
        )

    # Palantir normalizations

    if not inplace:
        adata = adata.copy()
    data_df = adata.to_df()
    if normalize:
        data_df = palantir.preprocess.normalize_counts(data_df)
        logg.info('data normalized ...')
    if log_transform:
        data_df = palantir.preprocess.log_transform(data_df)
        logg.info('data log transformed ...')
    if filter_low:
        data_df = palantir.preprocess.filter_counts_data(data_df)
        logg.info(
            'data filtered for low counts:\n'
            '\tcell_min_molecules=1000\n'
            '\tgenes_min_cells=10'
        )
    if normalize or log_transform or filter_low:
        adata.uns['palantir_norm_data'] = data_df

    # Processing

    logg.info('PCA in progress ...')
    pca_projections, var_r = palantir.utils.run_pca(data_df)
    adata.uns['palantir_pca_results'] = dict(
        pca_projections=pca_projections,
        variance_ratio=var_r,
    )

    logg.info('Diffusion maps in progress ...')
    dm_res = adata.uns['palantir_diff_maps'] = \
        palantir.utils.run_diffusion_maps(pca_projections)
    ms_data = adata.uns['palantir_ms_data'] = \
        palantir.utils.determine_multiscale_space(dm_res)

    logg.info('tSNE in progress ...')
    adata.uns['palantir_tsne'] = palantir.utils.run_tsne(ms_data)

    logg.info('imputation in progress ...')
    adata.uns['palantir_imp_df'] = \
        palantir.utils.run_magic_imputation(data_df, dm_res)

    logg.info('End of processing, start plotting.')
    return None if inplace else adata

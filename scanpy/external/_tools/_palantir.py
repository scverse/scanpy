"""Run Diffusion maps using the adaptive anisotropic kernel
"""

from anndata import AnnData
from scanpy import logging as logg


def palantir(adata: AnnData):
    """
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
        The `ms_data` - Multi scale data matrix

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

    *Load sample data*

    >>> adata = sc.read_csv(filename="Palantir/data/marrow_sample_scseq_counts.csv.gz")

    **Pre-processing**

    The provided adata will be used as input to the embedded `palantir` methods:

    >>> d = sce.tl.palantir( adata=adata )

    At this point, a new class object, `d`, will be instantiated. If the data
    needs pre-processing - filtering low genes/cells counts, or normalization,
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

    - note that a `start_cell` must be defined for every data set. The start cell for
    this dataset was chosen based on high expression of CD34.

    For further demonstration of palantir visualizations please follow this notebook
    `Palantir_sample_notebook.ipynb <https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst other things.

    """

    logg.info('Palantir diffusion maps', r=True)

    class _wrapper_cls(object):
        """
        A wrapper class to instantiate a new object that wraps `palantir` as an
        attribute reference attached to the class, together with other attribute
        references. The class uses instance variables, to preprocess and generate
        data using the embedded palantir package.
        Pre-processing of data is important step before start using the palantir
        methods.
        palantir accepts as input a Counts matrix: Cells x Genes.

        Methods used are:
            - instantiation initiation
            - instance function to embed palantir
            - pre-processing of input data
        """

        def __init__(self ,
                     adata,
                     func=None ,
                     normalize = False,
                     log_transform = False,
                     filter_low = False
                    ):
            """
            Parameters
            ----------
            adata : AnnData, or Dataframe of cells X genes
            func : function wrapper to import palantir (not to be used)
            normalize : `bool` (default: `False`)
                property setter passed to palantir to normalize using palantir method
                `palantir.preprocess.normalize_counts`.
            log_transform : `bool` (default: `False`)
                property setter passed to palantir. Some datasets show better signal in the log
                scale. Applied using `palantir.preprocess.log_transform`
            filter_low : `bool` (default: `False`)
                property setter passed to palantir to remove low molecule count cells and low detection genes
            """

            # instantiate variables
            self.func = func
            self.adata = adata
            self._normalize = normalize
            self._log_transform = log_transform
            self._filter_low = filter_low

            try:
                # for AnnData
                self.data_df = self.adata.to_df()
            except AttributeError:
                # assume the data is a cell X genes Dataframe
                logg.info('Assuming the data is a cell X genes Dataframe',
                	      r=True)

            # load palantir
            self.__call__()
            logg.info('palantir loaded ...', r=True)

        def __call__(self):
            """
            Call for function to import palantir and instantiate it as a class
            attribute
            """
            self.palantir = self.func()

        def process(self):

            """
            A method to run `palantir` on input Data Frame
            """

            # Principal component analysis
            logg.info('PCA in progress ...', r=True)

            self.pca_projections, self.var_r = self.palantir.utils.run_pca(self.data_df)

            adata.uns['palantir_pca_results'] = {}
            adata.uns['palantir_pca_results']['pca_projections'] = self.pca_projections
            adata.uns['palantir_pca_results']['variance_ratio'] = self.var_r

            # Diffusion maps
            logg.info('Diffusion maps in progress ...', r=True)

            self.dm_res = self.palantir.utils.run_diffusion_maps(self.pca_projections)
            self.ms_data = self.palantir.utils.determine_multiscale_space(self.dm_res)

            adata.uns['palantir_diff_maps'] = self.dm_res
            adata.uns['palantir_ms_data'] = self.ms_data

            # tSNE visualization
            logg.info('tSNE in progress ...', r=True)

            self.tsne = self.palantir.utils.run_tsne(self.ms_data)

            adata.uns['palantir_tsne'] = self.tsne

            # MAGIC imputation
            logg.info('imputation in progress ...', r=True)

            self.imp_df = self.palantir.utils.run_magic_imputation(self.data_df, self.dm_res)

            adata.uns['palantir_imp_df'] = self.imp_df

            logg.info('End of processing, start plotting.', r=True)

        @property
        def normalize(self):
            return self._normalize
        @normalize.setter
        def normalize(self , value):
            if value is True:
                self.data_df = self.palantir.preprocess.normalize_counts(self.data_df)
                adata.uns['palantir_norm_data'] = self.data_df
                logg.info('data normalized ...', r=True)

        @property
        def log_transform(self):
            return self._log_transform
        @log_transform.setter
        def log_transform(self , value):
            if value is True:
                self.data_df = self.palantir.preprocess.log_transform(self.data_df)
                adata.uns['palantir_norm_data'] = self.data_df
                logg.info('data log transformed ...', r=True)

        @property
        def filter_low(self):
            return self._filter_low
        @filter_low.setter
        def filter_low(self , value):
            if value is True:
                self.data_df = self.palantir.preprocess.filter_counts_data(self.data_df)
                adata.uns['palantir_norm_data'] = self.data_df
                logg.info('data filtered for low counts:\n\t' +\
                          'cell_min_molecules=1000\n\tgenes_min_cells=10',
                          r=True)


    def wrapper_cls(adata, func=None):
        """
        Class wrapper to pass a function to the class alongside positional argument
        """
        if func:
            return _wrapper_cls(func)
        else:
            def wrapper(func):
                return _wrapper_cls(adata, func)
            return wrapper

    # import palantir and wrap it in a function passed to the wrapper class
    # this method allows passing positional argument of adata to `_wrapper_cls`
    @wrapper_cls(adata)
    def _run():
        import importlib
        try:
            palantir = importlib.import_module('palantir')
        except ImportError:
            raise ImportError(
                '\nplease install palantir: \n\n\t'
                'git clone git://github.com/dpeerlab/Palantir.git\n\t'
                'cd Palantir\n\t'
                'sudo -H pip3 install .\n')
        return palantir
    return _run

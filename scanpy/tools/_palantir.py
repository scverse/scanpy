"""Run Diffusion maps using the adaptive anisotropic kernel
"""
from .. import logging as logg

def palantir( adata, **kargs ):
    
    """
    Run Diffusion maps using the adaptive anisotropic kernel [Setty27]_.

    Palantir is an algorithm to align cells along differentiation trajectories. 
    Palantir models differentiation as a stochastic process where stem cells 
    differentiate to terminally differentiated cells by a series of steps through
    a low dimensional phenotypic manifold. Palantir effectively captures the continuity 
    in cell states and the stochasticity in cell fate determination. 
    Palantir has been designed to work with multidimensional single cell data from diverse
    technologies such as Mass cytometry and single cell RNA-seq.

    Full documentation can be found here https://github.com/dpeerlab/Palantir

    :param adata: :class:`~anndata.AnnData`, or Dataframe of cells X genes\n

    :param normalize: `bool` (default: `False`) provide argument
                      for raw counts to normalize using palantir method, 
                      `palantir.preprocess.normalize_counts`

    :param log_transform: `bool` (default: `False`) some datasets show 
                      better signal in the log scale, applied using 
                      `palantir.preprocess.log_transform`

    :return:
              `data_df` - DataFrame of normalized adata\n
              `pca_projections` - PCA projections of the data\n
              `var_r` - PCA explained variance ratio\n
              `dm_res` - Diffusion components, corresponding eigen values 
                         and the diffusion operator\n
              `ms_data` - Multi scale data matrix\n
              `tsne` - tSNE embedding of the data\n
              `imp_df` - Imputed data matrix\n

    Return objects will be pushed to adata.

    Example
    -------

    >>> import scanpy.external as sce
    >>> import scanpy as sc

    A sample data is available at https://github.com/dpeerlab/Palantir/tree/master/data
    To view the plots, it is recommended to run Jupyter notebook
    
    *Load sample data*
    
    >>> adata = sc.read_csv(filename="Palantir/data/marrow_sample_scseq_counts.csv.gz")
    
    **Data preprocessing**
    
    At this point, a new class object, `d`, will be instantiated. The provided adata will be 
    processed and normalized, and then passed to palantir methods to generate the return objects 
    listed above. The generated objects will be pushed to the `adata` and stored for further use.
    Once instantiated, *Principal component analysis*, *Diffusion maps*, *tSNE on Diffusion maps*,
    and *MAGIC imputation* data objects will be created using the `palantir` default parameters.
    
    >>> d = sce.tl.palantir( adata=adata, normalize = True, log_transform = True )
    
    The created object `d.palantir` can be used to override the default parameters
    and re-run palantir different methods with the preferred parameters. 
    
    **Plotting**
    
    *tSNE visualization*
    
    >>> fig, ax = d.palantir.plot.plot_tsne(d.tsne)
    >>> fig, ax = d.palantir.plot.plot_tsne_by_cell_sizes(d.adata.X, d.tsne)
    
    *Gene expression can be visualized on tSNE maps*
    
    >>> d.palantir.plot.plot_gene_expression(d.imp_df, d.tsne, ['CD34', 'MPO', 'GATA1', 'IRF8'])
    
    *Diffusion maps*
    
    >>> d.palantir.plot.plot_diffusion_components(d.tsne, d.dm_res)
    
    For further demonstration of palantir visualizations please follow the link below:
    https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
    that also provides a guide to draw *gene expression trends* amongst other things.
    
    """
    import numpy as np

    logg.info('Palantir diffusion maps', r=True)

    class _wrapper_cls( object ):
        """
        A wrapper class to instantiate a new object that wraps `palantir` as an attribute 
        reference attached to the class, together with other attribute references.
        The class uses instance variables, to preprocess and generate data using the 
        embedded palantir package.
        Pre-processing of data is important step before start using the palantir methods. 
        palantir accepts as input a Counts matrix: Cells x Genes.
        
        Methods used are:
            - instantiation initiation
            - instance function to embed palantir
            - pre-processing of input data
        """

        def __init__( self ,
                                adata,
                                func=None ,
                                normalize = False,
                                log_transform = False,
                    ):
            
            """
            :input adata:         AnnData, or Dataframe of cells X genes
            :input func:          function wrapped, in this case a function to import palantir (not to be used)
            :input normalize:     `bool` (default: `False`), passed to palantir to normalize raw counts
            :input log_transform: `bool` (default: `False`), passed to palantir to log transform raw counts
            """
            # instantiate variables
            self.func = func
            self.adata = adata

            # load palantir
            self.__call__()
            logg.info('palantir loaded ...', r=True)

            # load data and normalize if necessary
            self.preprocessing( self.adata,
                                normalize = normalize,
                                log_transform = log_transform )

            adata.obsm['palantir_norm_data'] = np.array(self.data_df)

            # Principal component analysis
            logg.info('PCA in progress ...', r=True)
            self.pca_projections, self.var_r = self.palantir.utils.run_pca(self.data_df)
            adata.uns['palantir_pca'] = {}
            adata.uns['palantir_pca']['pca_projections'] = self.pca_projections
            adata.uns['palantir_pca']['variance_ratio'] = self.var_r

            # Diffusion maps
            logg.info('Diffusion maps in progress ...', r=True)
            self.dm_res = self.palantir.utils.run_diffusion_maps(self.pca_projections)
            self.ms_data = self.palantir.utils.determine_multiscale_space(self.dm_res)
            adata.uns['palantir_diff_maps'] = {}
            adata.uns['palantir_diff_maps']['dm_res'] = self.dm_res
            adata.obsm['X_palantir_diffmap'] = np.array(self.ms_data)

            # tSNE visualization
            logg.info('tSNE in progress ...', r=True)
            self.tsne = self.palantir.utils.run_tsne(self.ms_data)
            adata.obsm['X_palantir_tsne'] = np.array(self.tsne)

            # MAGIC imputation
            logg.info('imputation in progress ...', r=True)
            self.imp_df = self.palantir.utils.run_magic_imputation(self.data_df, self.dm_res)
            adata.obsm['X_palantir_imputation'] = np.array(self.imp_df)

            logg.info('End of processing, start plotting.', r=True)
            

        def __call__( self ):
            """
            Call for function to import palantir and instantiate it as a class attribute
            """
            self.palantir = self.func()
            

        def preprocessing( self, data_df = None,
                                 normalize = False,
                                 log_transform = False):
            try:
                # for AnnData
                self.data_df = self.adata.to_df()
            except AttributeError:
                # assume the data is a cell X genes Dataframe
                logg.info('Assuming the data is a cell X genes Dataframe', r=True)
                
            # Remove low molecule count cells and low detection genes
            self.data_df = self.palantir.preprocess.filter_counts_data(self.data_df)
            logg.info('Data loaded ...', r=True)
           
            if normalize:
                try:
                    self.data_df = self.palantir.preprocess.normalize_counts(self.data_df)
                    logg.info('data normalized ...', r=True)
                except AttributeError as e:
                    raise AttributeError( "Missing data_df: " + str(e) )
                    logg.error("Missing data_df: " + str(e), r=True)
                    
            if log_transform:
                try:
                    self.data_df = self.palantir.preprocess.log_transform(self.data_df)
                    print("data log transformed ...")
                except AttributeError as e:
                    raise AttributeError( "Missing data_df: " + str(e) )
                    logg.error("Missing data_df: " + str(e), r=True)


    def wrapper_cls( adata, func=None , **kargs):
        """
        Class wrapper to pass a function and instantiate arguments passed to the class
        """
        if func:
            return _wrapper_cls( func )
        else:
            def wrapper( func ):
                return _wrapper_cls( adata, func , **kargs)
            return wrapper

    # import palantir and wrap it in a function passed to the wrapper class
    # this method allows passing positional argument of adata, and additional
    # arguments used by the different methods of `_wrapper_cls`
    @wrapper_cls( adata , **kargs )
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

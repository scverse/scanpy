"""Harmony for data visualization with augmented affinity matrix at discrete timepoints
"""
from .. import logging as logg
import pandas as pd
from collections import OrderedDict
from itertools import chain
import numpy as np
import os

def harmony( adata, **kargs ):

    """Harmony for data visualization with augmented affinity matrix at discrete timepoints [Setty28]_.

    Harmony is a framework for data visualization, trajectory detection and interpretation
    for scRNA-seq data measured at discrete timepoints. Harmony constructs an augmented
    affinity matrix by augmenting the kNN graph affinity matrix with mutually nearest
    neighbors between successive time points. This augmented affinity matrix forms the
    basis for generated a force directed layout for visualization and also serves as input
    for computing the diffusion operator which can be used for trajectory detection using
    `Palantir`.

    Full documentation can be found here https://github.com/dpeerlab/Harmony

    :param adata: `list` of :class:`~anndata.AnnData`, or `.csv` filenames of count matrices
                  from two time points. Please ensure that replicates of the same time
                  point are consecutive in order
    :param timepoints: `list` (default: `None`), timepoints at which each cell was measured.
                       Important for Harmony augmented affinity matrix
    :param sample_names: `list` (default: `None`), sample names in the form of timepoint_rep#
                         of length `adata`

    :return:
        `harmony_adata` - :class:`~anndata.AnnData` of normalized data with the following objects:

            `harmony_adata.obsm['layout']` --> force directed layout

            `harmony_adata.uns['aff']` --> affinity matrix

            `harmony_adata.uns['aug_aff']` --> augmented affinity matrix

            `harmony_adata.uns['timepoint_connections']` --> timepoint connections

    Example
    -------

    >>> import scanpy.external as sce
    >>> import scanpy as sc
    >>> import os

    A sample data is available at https://github.com/dpeerlab/Harmony/tree/master/data

    To view the plots, it is recommended to run jupyter notebook

    *Load sample data*

    >>> # harmony_dir = os.path.expanduser('~/Harmony/data/')
    >>> csv_files = [   harmony_dir + 'Lib1-3_E3.5.csv',
                        harmony_dir + 'Lib1-4_E3.5.csv',
                        harmony_dir + 'Lib1-1_E4.5.csv',
                        harmony_dir + 'Lib1-2_E4.5.csv'    ]
    >>> sample_names = ['sa1_Rep1', 'sa1_Rep2',
                        'sa2_Rep1', 'sa2_Rep2']
    >>> # timepoints = ['sa1', 'sa2', 'sa3', 'sa4']
    >>> d = sce.tl.harmony(adata=csv_files , sample_names=sample_names)

    At this point, a new class object, `d`, will be instantiated

    **Pre-processing**

    Data normalization is performed internally once `sce.tl.harmony` is instantiated.

    The data can be log transformed by setting `log_transform` to `True`

    >>> d.log_transform = True

    The created object `d.harmony` can be used to override the default parameters as well
    as running the different methods embedded in it for plotting

    Follow the next step to pass the data to harmony methods, to generate the
    return objects listed above.

    **Run Harmony**

    >>> d.process()

    By calling this method `harmony` will run and generate the **Harmony augmented affinity matrix**.
    The generated objects will be pushed to `harmony_adata` and stored for further use.

    **Visualization using force directed layouts**

    Use Harmony methods for plotting. Example:

    >>> d.harmony.plot.plot_timepoints(d.layout, d.tp)

    For further demonstration of Harmony visualizations please follow the link
    below:

    https://github.com/dpeerlab/Harmony/blob/master/notebooks/Harmony_sample_notebook.ipynb

    that provides a comprehensive guide to draw *gene expression trends*, amongst other things.

    """

    logg.info('Harmony augmented affinity matrix', r=True)

    class _wrapper_cls( object ):

        """
        A wrapper class to instantiate a new object that wraps `Harmony` as an
        attribute reference attached to the class, together with other attribute
        references. The class uses instance variables, to preprocess and generate
        data using the embedded harmony package.

        harmony accepts as input a list of Counts matrices: Cells x Genes.

        Methods used are:
            - instantiation initiation
            - instance function to embed harmony
            - processing of input data
        """

        def __init__( self ,
                                adata,
                                func = None ,
                                timepoints = None,
                                sample_names = None,
                                log_transform = False
                    ):
            """
            :param adata: `list` of :class:`~anndata.AnnData`, or `.csv` filenames of count
                          matrices from two time points. Please ensure that replicates of
                          the same time point are consecutive in order
            :input func: function wrapper to import harmony (not to be used)
            :input timepoints: `list` (default: `None`), timepoints at which each cell was
                               measured. Important for Harmony augmented affinity matrix
            :input sample_names: `list` (default: `None`), sample names in the form of
                                 timepoint_rep# of length `adata`
            :input log_transform: `bool` (default: `False`), property setter passed to harmony

            """

            # instantiate variables
            self.func = func
            self.adata = adata
            self.timepoints = timepoints
            self.sample_names = sample_names
            self._log_transform = log_transform

            # load harmony
            self.__call__()
            logg.info('harmony loaded ...', r=True)

            # determine if adata is a .csv file or a scanpy.AnnData object
            # and load the counts matrices into the instantiated attribute
            if self.sample_names is None:
                print( '"sample_names" not provided, constructed internally')
                rep = range(len(self.adata))
                self.sample_names = ['_'.join(['sample{}'.format(n+1),str(n+1)]) for n in rep]
                print( '"sample_names": {}'.format(self.sample_names) )

            try:
                assert(self.timepoints)
            except AssertionError:
                try:
                    self.timepoints = [i.split("_")[0] for i in self.sample_names]
                except:
                    raise AssertionError('"timepoints" must be a list of Time Points!!\n\n\t' +\
                                         '"timepoints" = {}'.format( self.timepoints ) +\
                                         '\n\tset the values of Time Points and try again!')
            # keep the order of unique time points
            self.timepoints = list(OrderedDict.fromkeys(self.timepoints))

            try:
                assert(os.path.exists(self.adata[0]))
                counts = self.harmony.utils.load_from_csvs(self.adata, self.sample_names)
                logg.info('Input data is a list of .csv files', r=True)
            except:
                try:
                    if not isinstance(adata[0], pd.DataFrame):
                        assert( isinstance( adata[0].to_df(), pd.DataFrame) )
                    else:
                        assert( isinstance( adata[0], pd.DataFrame) )
                    counts = self.load_from_AnnData(self.adata, self.sample_names)
                    logg.info('Input data is a list of scanpy.AnnData objects', r=True)
                except:
                    raise RuntimeError('Input data must be a list of .csv files,' +\
                                       ' scanpy.AnnData, or pandas DataFrames')

            # PREPROCESS
            # normalize data
            norm_df = self.harmony.utils.normalize_counts(counts)

            # select highly variable genes
            hvg_genes = self.harmony.utils.hvg_genes(norm_df)

            self.data_df = norm_df.loc[:,hvg_genes]
            logg.info('data normalized and returned as annData ...', r=True)

        def __call__( self ):
            """
            Call for function to import harmony and instantiate it as a class
            attribute
            """
            self.harmony = self.func()


        def process( self ):

            """
            A method to run `harmony` on input Data Frame
            """

            # Harmony augmented affinity matrix
            logg.info('Harmony augmented affinity matrix ...', r=True)

            self.tp = pd.Series(index=self.data_df.index)
            for t in self.timepoints:
                cells = self.data_df.index[self.data_df.index.str.contains(t)]
                self.tp[cells] = t
            self.timepoint_connections = pd.DataFrame(columns=[0, 1])
            index = 0
            for i in range(len(self.timepoints)-1):
                self.timepoint_connections.loc[index, :] = self.timepoints[i:i+2]
                index += 1

            # compute the augmented and non-augmented affinity matrices
            self.aug_aff, self.aff = self.harmony.core.augmented_affinity_matrix(
                                                    self.data_df,
                                                    self.tp,
                                                    self.timepoint_connections
                                                    )

            # Visualization using force directed layouts
            self.layout = self.harmony.plot.force_directed_layout(self.aug_aff,
                                                                  self.data_df.index)

            # push outputs to a new scanpy.AnnDana
            from scanpy import AnnData
            self.harmony_adata = AnnData( self.data_df )
            self.harmony_adata.obsm['layout'] = np.array(self.layout)
            self.harmony_adata.uns['tp'] = self.tp
            self.harmony_adata.uns['aff'] = self.aff
            self.harmony_adata.uns['aug_aff'] = self.aug_aff
            self.harmony_adata.uns['sample_names'] = self.sample_names
            self.harmony_adata.uns['timepoints'] = self.timepoints
            self.harmony_adata.uns['timepoint_connections'] = self.timepoint_connections

            logg.info('End of processing, start plotting.', r=True)

            return self.harmony_adata


        def load_from_AnnData(self,
                                    adata_list,
                                    sample_names = None,
                                    min_cell_count = 10
                             ):
            """
            Read in scRNA-seq data from :class:`~anndata.AnnData` list

            :param adata_list: `List` of :class:`~anndata.AnnData`, or Dataframes of cells X genes
            :param sample_names: `List` of prefixes to attach to the cell names
            :param min_cell_count: Minimum number of cells a gene should be detected in
            :return: Pandas data frame representing the count matrix

            """

            # Load counts
            print('Loading count matrices...')
            counts_dict = OrderedDict()

            for adata, sample in zip(adata_list, sample_names):
                print(sample)
                try:
                    # for AnnData
                    counts = adata.to_df()
                except AttributeError:
                    # assume the data is a cell X genes Dataframe
                    logg.info('Assuming the data is a cell X genes Dataframe',
                	      r=True)
                    counts = adata

                # Update sample names
                counts.index = sample + '_' + counts.index.astype(str)

                # Remove zero count genes
                counts = counts.loc[:, counts.sum() > 0]
                counts = counts.astype(np.int16)

                # Update dictionary
                counts_dict[sample] = counts

            # Concatenate cells and genes
            print('Concatenating data..')
            all_cells = list(chain(*[list(counts_dict[sample].index)
                                     for sample in sample_names]))
            all_genes = list(
                chain(*[list(counts_dict[sample].columns) for sample in sample_names]))
            all_genes = list(set(all_genes))

            # Counts matrix
            counts = pd.DataFrame(0, index=all_cells,
                                  columns=all_genes, dtype=np.int16)
            for sample in sample_names:
                sample_counts = counts_dict[sample]
                counts.loc[sample_counts.index, sample_counts.columns] = sample_counts

            # Filter out low detection genes
            gs = counts.sum()
            counts = counts.loc[:, counts.columns[gs > min_cell_count]]

            return counts

        @property
        def log_transform( self ):
            return self._log_transform
        @log_transform.setter
        def log_transform( self , value ):
            if value is True:
                self.data_df = self.harmony.utils.log_transform(self.data_df)
                logg.info('data log transformed ...', r=True)



    def wrapper_cls( adata, func=None, **kargs ):
        """
        Class wrapper to pass a function to the class alongside positional argument
        """
        if func:
            return _wrapper_cls( func )
        else:
            def wrapper( func ):
                return _wrapper_cls( adata, func, **kargs )
            return wrapper

    # import palantir and wrap it in a function passed to the wrapper class
    # this method allows passing positional argument of adata to `_wrapper_cls`
    @wrapper_cls( adata , **kargs )
    def _run():
        import importlib
        try:
            import palantir
        except:
            raise ImportError(
                ' please install palantir: \n\n\t'
                'git clone git://github.com/dpeerlab/Palantir.git\n\t'
                'cd Palantir\n\t'
                'sudo -H pip3 install .\n')
        try:
            harmony = importlib.import_module('harmony')
        except ImportError:
            raise ImportError(
                '\nplease install harmony: \n\n\t'
                'git clone git://github.com/dpeerlab/Harmony.git\n\t'
                'cd Harmony\n\t'
                'sudo -H pip3 install .\n')
        return harmony
    return _run

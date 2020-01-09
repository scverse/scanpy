"""Harmony time series for data visualization with augmented affinity matrix at discrete timepoints
"""

from scanpy import logging as logg
import pandas as pd
from collections import OrderedDict
from itertools import chain
import numpy as np
import os


def harmony_timeseries(data: list, **kargs):
    """
    Harmony time series for data visualization with augmented affinity matrix at discrete timepoints [Nowotschin18i]_.

    Harmony time series is a framework for data visualization, trajectory detection and interpretation
    for scRNA-seq data measured at discrete timepoints. Harmony constructs an augmented
    affinity matrix by augmenting the kNN graph affinity matrix with mutually nearest
    neighbors between successive time points. This augmented affinity matrix forms the
    basis for generated a force directed layout for visualization and also serves as input
    for computing the diffusion operator which can be used for trajectory detection using
    `Palantir <https://github.com/dpeerlab/Palantir>`__.

    .. note::
        More information and bug reports `here <https://github.com/dpeerlab/Harmony>`__.

    Parameters
    ----------
    data
        list of :class:`~anndata.AnnData`, or cell x gene `.csv` filenames of count matrices from two time
        points. Please ensure that replicates of the same time point are consecutive in order

    timepoints
        defaults to None, a list of timepoints at which each cell was measured. Important for Harmony augmented affinity matrix

    sample_names
        defaults to None, a list sample names in the form of timepoint_rep# of length `data`

    log_transform
        defaults to False, property setter passed to harmony


    Returns
    -------
    `harmony_adata` of :class:`~anndata.AnnData` of normalized data:

    - `.obsm['layout']` force directed layout

    - `.uns['aff']` affinity matrix

    - `.uns['aug_aff']` augmented affinity matrix

    - `.uns['timepoint_connections']` timepoint connections

    Example
    -------

    >>> import scanpy.external as sce
    >>> import scanpy as sc
    >>> import os

    There are two ways to load data used with harmony_timeseries. If you have `csv`
    files, follow these steps. If you have `AnnData` sets, follow the next snippet.

    *Load sample data from csv files*

    A sample data is available `here <https://github.com/dpeerlab/Harmony/tree/master/data>`_.

    To view the plots, it is recommended to run jupyter notebook

    >>> # harmony_dir = os.path.expanduser('~/Harmony/data/')
    >>> csv_files = [
    ...     harmony_dir + i for i in [
    ...         'Lib1-3_E3.5.csv',
    ...         'Lib1-4_E3.5.csv',
    ...         'Lib1-1_E4.5.csv',
    ...         'Lib1-2_E4.5.csv'
    ...     ]
    ... ]
    >>> sample_names = [
    ...     'sa1_Rep1',
    ...     'sa1_Rep2',
    ...     'sa2_Rep1',
    ...     'sa2_Rep2'
    ... ]
    >>> # timepoints = ['sa1', 'sa2', 'sa3', 'sa4']
    >>> d = sce.tl.harmony_timeseries(data=csv_files, sample_names=sample_names)

    At this point, a new class object, `d`, will be instantiated

    *Loading sample data from `AnnData` sets*

    Assuming you have 4 `AnnData` sets from 2 time points with 2 replicates,
    instead of `csv` files, the following snippet shows how to feed a list of
    `AnnData` sets to harmony_timeseries.

    Loading sample data from Scanpy builtin datasets:

    >>> adata_ref = sc.datasets.pbmc3k()

    Splitting `adata_ref` into 4 `adata` sets (just to show how it works). Normally,
    you would have different `AnnData` sets from different time points ready to use.

    >>> n = int(adata_ref.shape[0] / 4)
    >>> datalist = [adata_ref[n * i : n * (i + 1)] for i in range(4)]
    >>> sample_names = [
    ...     'sa1_Rep1',
    ...     'sa1_Rep2',
    ...     'sa2_Rep1',
    ...     'sa2_Rep2'
    ... ]
    >>> # Instantiate harmony_timeseries:
    >>> d = sce.tl.harmony_timeseries(data=datalist, sample_names=sample_names)

    The next steps apply on both cases of loaded data.

    **Pre-processing**

    Data normalization is performed internally once `sce.tl.harmony_timeseries` is instantiated.

    The data can be log transformed by setting `log_transform` to `True`

    >>> d.log_transform = True

    The created object `d.harmony_timeseries` can be used to override the default parameters as well
    as running the different methods embedded in it for plotting

    Follow the next step to pass the data to harmony methods, to generate the
    return objects listed above.

    **Run harmony_timeseries**

    >>> d.process(n_components=902)

    By calling this method `harmony_timeseries` will run and generate the **Harmony augmented affinity matrix**.
    The generated objects will be pushed to `harmony_adata` and stored for further use.

    **Visualization using force directed layouts**

    Use Harmony methods for plotting. Example:

    >>> d.harmony_timeseries.plot.plot_timepoints(d.layout, d.tp)

    For further demonstration of Harmony visualizations please follow this notebook
    `Harmony_sample_notebook.ipynb <https://github.com/dpeerlab/Harmony/blob/master/notebooks/Harmony_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst other things.
    """

    logg.info("Harmony augmented affinity matrix")

    class _wrapper_cls(object):

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

        def __init__(
            self,
            data: list,
            func=None,
            timepoints: list = None,
            sample_names: list = None,
            log_transform: bool = False,
        ):
            """
            Parameters
            ----------

            data
                list of :class:`~anndata.AnnData`, or `.csv` filenames of count matrices from
                two time points. Please ensure that replicates of the same time point are
                consecutive in order
            func
                function wrapper to import harmony (not to be used)
            timepoints
                default: `None`, timepoints at which each cell was measured.
                Important for Harmony augmented affinity matrix
            sample_names
                default: `None`, sample names in the form of timepoint_rep# of length `data`
            log_transform
                default: `False`, property setter passed to harmony
            """

            # instantiate variables
            self.func = func
            self.data = data
            self.timepoints = timepoints
            self.sample_names = sample_names
            self._log_transform = log_transform

            # load harmony_timeseries
            self.__call__()
            logg.info("harmony_timeseries loaded ...")

            # determine if data is a .csv file or a scanpy.AnnData object
            # and load the counts matrices into the instantiated attribute
            if self.sample_names is None:
                print('"sample_names" not provided, constructed internally')
                rep = range(len(self.data))
                self.sample_names = [
                    "_".join(["sample{}".format(n + 1), str(n + 1)])
                    for n in rep
                ]
                print('"sample_names": {}'.format(self.sample_names))

            try:
                assert self.timepoints
            except AssertionError:
                try:
                    self.timepoints = [
                        i.split("_")[0] for i in self.sample_names
                    ]
                except:
                    msg = (
                        '"timepoints" must be a list of Time Points!!\n\n\t'
                        + '"timepoints" = {}'.format(self.timepoints)
                        + "\n\tset the values of Time Points and try again!"
                    )
                    raise AssertionError(msg)
            # keep the order of unique time points
            self.timepoints = list(OrderedDict.fromkeys(self.timepoints))

            try:
                assert os.path.exists(self.data[0])
                counts = self.harmony_timeseries.utils.load_from_csvs(
                    self.data, self.sample_names
                )
                logg.info("Input data is a list of .csv files")
            except:
                try:
                    if not isinstance(data[0], pd.DataFrame):
                        assert isinstance(data[0].to_df(), pd.DataFrame)
                    else:
                        assert isinstance(data[0], pd.DataFrame)
                    counts = self.load_from_AnnData(
                        self.data, self.sample_names
                    )
                    logg.info("Input data is a list of scanpy.AnnData objects")
                except:
                    raise RuntimeError(
                        "Input data must be a list of .csv files,"
                        + " scanpy.AnnData, or pandas DataFrames"
                    )

            # PREPROCESS
            # normalize data
            norm_df = self.harmony_timeseries.utils.normalize_counts(counts)

            # select highly variable genes
            hvg_genes = self.harmony_timeseries.utils.hvg_genes(norm_df)

            self.data_df = norm_df.loc[:, hvg_genes]
            logg.info("data normalized and returned as annData ...")

        def __call__(self):
            """
            Call for function to import harmony_timeseries and instantiate it as a class
            attribute
            """
            self.harmony_timeseries = self.func()

        def process(self, n_components: int = 1000):

            """
            A method to run `harmony_timeseries` on input Data Frame
            :param n_components: int
            """

            # harmony_timeseries augmented affinity matrix
            logg.info("harmony_timeseries augmented affinity matrix ...")

            self.tp = pd.Series(index=self.data_df.index)
            for t in self.timepoints:
                cells = self.data_df.index[self.data_df.index.str.contains(t)]
                self.tp[cells] = t
            self.timepoint_connections = pd.DataFrame(columns=[0, 1])
            index = 0
            for i in range(len(self.timepoints) - 1):
                self.timepoint_connections.loc[index, :] = self.timepoints[
                    i : i + 2
                ]
                index += 1

            # compute the augmented and non-augmented affinity matrices
            (
                self.aug_aff,
                self.aff,
            ) = self.harmony_timeseries.core.augmented_affinity_matrix(
                self.data_df,
                self.tp,
                self.timepoint_connections,
                pc_components=n_components,
            )

            # Visualization using force directed layouts
            self.layout = self.harmony_timeseries.plot.force_directed_layout(
                self.aug_aff, self.data_df.index
            )

            # push outputs to a new scanpy.AnnDana
            from scanpy import AnnData

            self.harmony_adata = AnnData(self.data_df)
            self.harmony_adata.obsm["layout"] = np.array(self.layout)
            self.harmony_adata.uns["tp"] = self.tp
            self.harmony_adata.uns["aff"] = self.aff
            self.harmony_adata.uns["aug_aff"] = self.aug_aff
            self.harmony_adata.uns["sample_names"] = self.sample_names
            self.harmony_adata.uns["timepoints"] = self.timepoints
            self.harmony_adata.uns[
                "timepoint_connections"
            ] = self.timepoint_connections

            logg.info("End of processing, start plotting.")

            return self.harmony_adata

        def load_from_AnnData(
            self,
            data_list: list,
            sample_names: list = None,
            min_cell_count: int = 10,
        ):
            """
            Read in scRNA-seq data from :class:`~anndata.AnnData` list

            Parameters
            ----------
            data_list
                of :class:`~anndata.AnnData`, or Dataframes of cells X genes
            sample_names
                default: `None`, list of prefixes to attach to the cell names
            min_cell_count
                default: 10, Minimum number of cells a gene should be detected in

            Returns
            ----------
            Pandas data frame representing the count matrix
            """

            # Load counts
            print("Loading count matrices...")
            counts_dict = OrderedDict()

            for data, sample in zip(data_list, sample_names):
                print(sample)
                try:
                    # for AnnData
                    counts = data.to_df()
                except AttributeError:
                    # assume the data is a cell X genes Dataframe
                    logg.info("Assuming the data is a cell X genes Dataframe")
                    counts = data

                # Update sample names
                counts.index = sample + "_" + counts.index.astype(str)

                # Remove zero count genes
                counts = counts.loc[:, counts.sum() > 0]
                counts = counts.astype(np.int16)

                # Update dictionary
                counts_dict[sample] = counts

            # Concatenate cells and genes
            print("Concatenating data..")
            all_cells = list(
                chain(
                    *[
                        list(counts_dict[sample].index)
                        for sample in sample_names
                    ]
                )
            )
            all_genes = list(
                chain(
                    *[
                        list(counts_dict[sample].columns)
                        for sample in sample_names
                    ]
                )
            )
            all_genes = list(set(all_genes))

            # Counts matrix
            counts = pd.DataFrame(
                0, index=all_cells, columns=all_genes, dtype=np.int16
            )
            for sample in sample_names:
                sample_counts = counts_dict[sample]
                counts.loc[
                    sample_counts.index, sample_counts.columns
                ] = sample_counts

            # Filter out low detection genes
            gs = counts.sum()
            counts = counts.loc[:, counts.columns[gs > min_cell_count]]

            return counts

        @property
        def log_transform(self):
            return self._log_transform

        @log_transform.setter
        def log_transform(self, value):
            if value is True:
                self.data_df = self.harmony_timeseries.utils.log_transform(
                    self.data_df
                )
                logg.info("data log transformed ...")

    def wrapper_cls(data, func=None, **kargs):
        """
        Class wrapper to pass a function to the class alongside positional argument
        """
        if func:
            return _wrapper_cls(func)
        else:

            def wrapper(func):
                return _wrapper_cls(data, func, **kargs)

            return wrapper

    # import palantir and wrap it in a function passed to the wrapper class
    # this method allows passing positional argument of data to `_wrapper_cls`
    @wrapper_cls(data, **kargs)
    def _run():
        import importlib

        try:
            import palantir
        except ImportError:
            raise ImportError(
                " please install palantir: \n\n\t"
                "git clone git://github.com/dpeerlab/Palantir.git\n\t"
                "cd Palantir\n\t"
                "sudo -H pip3 install .\n"
            )
        try:
            harmony = importlib.import_module("harmony")
        except ImportError:
            raise ImportError(
                "\nplease install harmony: \n\n\t"
                "git clone git://github.com/dpeerlab/Harmony.git\n\t"
                "cd Harmony\n\t"
                "sudo -H pip3 install .\n"
            )
        return harmony

    return _run

"""Harmony time series for data visualization with augmented affinity matrix at
    discrete time points
"""

import pandas as pd
from anndata import AnnData

from ... import logging as logg


def harmony_timeseries(data: AnnData, tp: str):
    """\
    Harmony time series for data visualization with augmented affinity matrix
    at discrete time points [Nowotschin18i]_.

    Harmony time series is a framework for data visualization, trajectory
    detection and interpretation for scRNA-seq data measured at discrete
    time points. Harmony constructs an augmented affinity matrix by augmenting
    the kNN graph affinity matrix with mutually nearest neighbors between
    successive time points. This augmented affinity matrix forms the basis for
    generated a force directed layout for visualization and also serves as input
    for computing the diffusion operator which can be used for trajectory
    detection using Palantir_.

    .. _Palantir: https://github.com/dpeerlab/Palantir

    More about **Palantir** can be found here:
    https://github.com/dpeerlab/Palantir.

    .. note::
        More information and bug reports `here
        <https://github.com/dpeerlab/Harmony>`__.

    Parameters
    ----------
    data
        Annotated data matrix of shape n_obs `×` n_vars. Rows correspond to
        cells and columns to genes. Rows represent two or more time points,
        where replicates of the same time point are consecutive in order.

    tp
        key name of observation annotation `.obs` representing time points

    Returns
    -------
    Updates `.uns` with `timepoint_connections`

    Example
    -------

    >>> import scanpy as sc
    >>> import scanpy.external as sce

    **Load** `AnnData`

    A sample with real data is available `here
    <https://github.com/dpeerlab/Harmony/tree/master/data>`_.

    Random data sets of three time points with two replicates each:

    >>> adata_ref = sc.datasets.pbmc3k()
    >>> start = [596, 615, 1682, 1663, 1409, 1432]
    >>> adatas = [adata_ref[i : i + 1000] for i in start]
    >>> sample_names = [
            "sa1_Rep1",
            "sa1_Rep2",
            "sa2_Rep1",
            "sa2_Rep2",
            "sa3_Rep1",
            "sa3_Rep2",
        ]
    >>> timepoints = [i.split("_")[0] for i in sample_names]
    >>> for ad, sn, tp in zip(adatas, sample_names, timepoints):
            ad.obs["time_points"] = tp
            ad.obs["sample_name"] = sn
    >>> adata = adatas[0].concatenate(*adatas[1:], join="outer")

    Normalize and filter for highly expressed genes

    >>> sc.pp.normalize_total(adata, target_sum=10000)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    Run harmony_timeseries

    >>> d = sce.tl.harmony_timeseries(data=adata, tp="time_points")

    **Harmony augmented affinity matrix**

    >>> aug_aff, aff = d.harmony_timeseries.core.augmented_affinity_matrix(
            data_df=adata.to_df(),
            timepoints=adata.obs["time_points"],
            timepoint_connections=adata.uns["timepoint_connections"],
        )

    **Visualization using force directed layouts**

    >>> layout = d.harmony_timeseries.plot.force_directed_layout(
            affinity_matrix=aug_aff, cell_names=adata.obs.index
        )

    Use any of Scanpy or Harmony methods for plotting. Example:

    >>> d.harmony_timeseries.plot.plot_timepoints(
            layout=layout, timepoints=adata.obs.time_points
        )

    For further demonstration of Harmony visualizations please follow this
    notebook
    `Harmony_sample_notebook.ipynb
    <https://github.com/dpeerlab/Harmony/blob/master/notebooks/
    Harmony_sample_notebook.ipynb>`_.
    It provides a comprehensive guide to draw *gene expression trends*, amongst
    other things.
    """

    logg.info("Harmony augmented affinity matrix")

    class _wrapper_cls(object):

        """/
        A wrapper class to instantiate a new object that wraps `Harmony` as an
        attribute reference attached to the class, together with other attribute
        references. The class uses instance variables, to preprocess and
        generate data using the embedded harmony package.

        harmony accepts as input AnnData: Cells x Genes.

        Methods used are:
            - instantiation initiation
            - instance function to embed harmony
            - processing of input data
        """

        def __init__(
            self, adata: AnnData, tp: str, func=None,
        ):
            """/
            Parameters
            ----------

            adata
                Annotated data matrix of shape n_obs `×` n_vars. Rows correspond
                to cells and columns to genes. Rows represent two or more time
                points, where replicates of the same time point are consecutive
                in order.
            tp
                key name of observation annotation `.obs` representing time
                points
            func
                function wrapper to import harmony (not to be used)
            """

            # instantiate variables
            self.func = func

            timepoint_connections = pd.DataFrame(columns=[0, 1])
            index = 0
            timepoints = adata.obs[tp].unique().tolist()
            for i in range(len(timepoints) - 1):
                timepoint_connections.loc[index, :] = timepoints[i : i + 2]
                index += 1
            adata.uns["timepoint_connections"] = timepoint_connections

            # load harmony_timeseries
            self.__call__()
            logg.info("harmony_timeseries loaded ...")

        def __call__(self):
            """/
            Call for function to import harmony_timeseries and instantiate it
            as a class attribute
            """
            self.harmony_timeseries = self.func()

    def wrapper_cls(df, tpoints, func=None):
        """/
        Class wrapper to pass a function to the class alongside positional
        argument
        """
        if func:
            return _wrapper_cls(func)
        else:

            def wrapper(f):
                return _wrapper_cls(df, tpoints, f)

            return wrapper

    # import Harmony and wrap it in a function passed to the wrapper class this
    # method allows passing positional argument of data to `_wrapper_cls`
    @wrapper_cls(data, tp)
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

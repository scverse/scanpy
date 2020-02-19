from typing import Union, List, Optional, Any, Tuple, Collection

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib.axes import Axes
from numpy.core.umath_tests import inner1d

from .._utils import _doc_params
from ..plotting import embedding
from ..plotting._docs import (
    doc_adata_color_etc,
    doc_edges_arrows,
    doc_scatter_embedding,
    doc_show_save_ax,
)
from ..plotting._tools.scatterplots import _wraps_plot_scatter
from ..plotting import _utils


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def phate(adata, **kwargs) -> Union[List[Axes], None]:
    """\
    Scatter plot in PHATE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False`, a list of :class:`~matplotlib.axes.Axes` objects.
    Every second element corresponds to the 'right margin'
    drawing area for color bars and legends.

    Examples
    --------
    >>> from anndata import AnnData
    >>> import scanpy.external as sce
    >>> import phate
    >>> data, branches = phate.tree.gen_dla(
    ...     n_dim=100,
    ...     n_branch=20,
    ...     branch_length=100,
    ... )
    >>> data.shape
    (2000, 100)
    >>> adata = AnnData(data)
    >>> adata.obs['branches'] = branches
    >>> sce.tl.phate(adata, k=5, a=20, t=150)
    >>> adata.obsm['X_phate'].shape
    (2000, 2)
    >>> sce.pl.phate(
    ...     adata,
    ...     color='branches',
    ...     color_map='tab20',
    ... )
    """
    return embedding(adata, 'phate', **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def trimap(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in TriMap basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return embedding(adata, 'trimap', **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def harmony_timeseries(
    adata, *, show: bool = True, return_fig: bool = False, **kwargs
) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in Harmony force-directed layout basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `return_fig` is True, a :class:`~matplotlib.figure.Figure`.
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """

    tp_name = adata.uns["harmony_timepoint_var"]
    tps = adata.obs[tp_name].unique()

    fig, axes = plt.subplots(1, len(tps))
    for i, tp in enumerate(tps):
        p = embedding(
            adata,
            'harmony',
            color=tp_name,
            groups=tp,
            title=tp,
            show=False,
            ax=axes[i],
            legend_loc='none',
        )
        p.set_axis_off()
    if return_fig:
        return fig
    elif not show:
        return axes


def sam(
    adata: AnnData,
    projection: Union[str, np.ndarray] = 'X_umap',
    c: Optional[Union[str, np.ndarray]] = None,
    cmap: str = 'Spectral_r',
    linewidth: float = 0.0,
    edgecolor: str = 'k',
    axes: Optional[Axes] = None,
    colorbar: bool = True,
    s: float = 10.0,
    **kwargs: Any,
) -> Axes:
    """\
    Scatter plot using the SAM projection or another input projection.

    Parameters
    ----------
    projection
        A case-sensitive string indicating the projection to display (a key
        in adata.obsm) or a 2D numpy array with cell coordinates. If None,
        projection defaults to UMAP.
    c
        Cell color values overlaid on the projection. Can be a string from adata.obs
        to overlay cluster assignments / annotations or a 1D numpy array.
    axes
        Plot output to the specified, existing axes. If None, create new
        figure window.
    kwargs
        all keyword arguments in matplotlib.pyplot.scatter are eligible.
    """

    if isinstance(projection, str):
        try:
            dt = adata.obsm[projection]
        except KeyError:
            raise ValueError(
                'Please create a projection first using run_umap or run_tsne'
            )
    else:
        dt = projection

    if axes is None:
        plt.figure()
        axes = plt.gca()

    if c is None:
        axes.scatter(
            dt[:, 0], dt[:, 1], s=s, linewidth=linewidth, edgecolor=edgecolor, **kwargs
        )
        return axes

    if isinstance(c, str):
        try:
            c = np.array(list(adata.obs[c]))
        except KeyError:
            pass

    if isinstance(c[0], (str, np.str_)) and isinstance(c, (np.ndarray, list)):
        from SAM import ut

        i = ut.convert_annotations(c)
        ui, ai = np.unique(i, return_index=True)
        cax = axes.scatter(
            dt[:, 0],
            dt[:, 1],
            c=i,
            cmap=cmap,
            s=s,
            linewidth=linewidth,
            edgecolor=edgecolor,
            **kwargs,
        )

        if colorbar:
            cbar = plt.colorbar(cax, ax=axes, ticks=ui)
            cbar.ax.set_yticklabels(c[ai])
    else:
        if not isinstance(c, (np.ndarray, list)):
            colorbar = False
        i = c

        cax = axes.scatter(
            dt[:, 0],
            dt[:, 1],
            c=i,
            cmap=cmap,
            s=s,
            linewidth=linewidth,
            edgecolor=edgecolor,
            **kwargs,
        )

        if colorbar:
            plt.colorbar(cax, ax=axes)
    return axes


def wishbone_marker_trajectory(
    adata: AnnData, no_bins: int = 150, smoothing_factor: int = 1,
):
    """\
    Plot marker trends along trajectory

    Parameters
    ----------
    adata
        Annotated data matrix.
    no_bins
        Number of bins for calculating marker density
    smoothing_factor
        Parameter controlling the degree of smoothing

    Returns
    -------
    Updates `adata` with the following fields:

    `weights_wishbone` : :class:`pandas.DataFrame` (`adata.uns`)
        Computed gaussian weights for points at each location
    `branch_point_bin` : 'int' (`adata.uns`)
        Identifies the bin with the branch point. In case of no branching,
        defaults to `no_bins`
    `bins_wishbone` : :class:`numpy.ndarray` (`adata.uns`)
        Computed bin locations and bin memberships
    `branches_wishbone`: :class:`numpy.ndarray` (`adata.uns`)
        In case of branching, returns a list of branches [2, 3].
    """

    # Compute bin locations and bin memberships
    # Sort trajectory
    trajectory = adata.obs['trajectory_wishbone'].sort_values()
    bins = np.linspace(np.min(trajectory), np.max(trajectory), no_bins)

    # Compute gaussian weights for points at each location
    # Standard deviation estimated from Silverman's approximation
    stdev = np.std(trajectory) * 1.34 * len(trajectory) ** (-1 / 5) * smoothing_factor
    weights = np.exp(
        -((np.tile(trajectory, [no_bins, 1]).T - bins) ** 2 / (2 * stdev ** 2))
    ) * (1 / (2 * np.pi * stdev ** 2) ** 0.5)

    # Adjust weights if data has branches
    if 'branch_wishbone' in adata.obs.keys():
        # Branch of the trunk
        trunk = adata.obs['branch_wishbone'][trajectory.index[0]]
        adata.uns['branches_wishbone'] = np.array([2, 3])
        # Counts of branch cells in each bin
        branch_counts = pd.DataFrame(np.zeros([len(bins) - 1, 3]), columns=[1, 2, 3])
        for j in branch_counts.columns:
            branch_counts[j] = pd.Series(
                [
                    pd.Series(
                        adata.obs['branch_wishbone'][
                            trajectory.index[
                                (trajectory > bins[i - 1]) & (trajectory < bins[i])
                            ]
                        ]
                        == j
                    ).sum()
                    for i in range(1, len(bins))
                ]
            )
        # Frequencies
        branch_counts = branch_counts.divide(branch_counts.sum(axis=1), axis=0)

        # Identify the bin with the branch point by looking at the weights
        weights = pd.DataFrame(weights, index=trajectory.index, columns=range(no_bins))
        bp_bin = weights.columns[np.where(branch_counts[trunk] < 0.9)[0][0]] + 0
        if bp_bin < 0:
            bp_bin = 3

    else:
        bp_bin = no_bins

    # return weight to adata
    adata.uns['weights_wishbone'] = weights
    adata.uns['branch_point_bin'] = bp_bin
    adata.uns['bins_wishbone'] = bins


@_doc_params(show_save_ax=doc_show_save_ax)
def wishbone(
    adata: AnnData,
    markers: List[str],
    show_variance: bool = False,
    no_bins: int = 150,
    smoothing_factor: int = 1,
    min_delta: float = 0.1,
    linestyle: List[Union[str, Tuple[int, Collection[int]]]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    legend_loc: Union[int, str] = 'upper left',
    legend_fontsize: Optional[Union[int, float, str]] = None,
    legend_prop: Optional[dict] = None,
    return_fig: bool = False,
    show: bool = True,
    save: Optional[Union[str, bool]] = None,
    **kargs,
) -> Union[Axes, None]:
    """\
    Plot marker trends along trajectory.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Iterable of markers/genes to be plotted.
    show_variance
        Logical indicating if the trends should be accompanied with variance.
    no_bins
        Number of bins for calculating marker density.
    smoothing_factor
        Parameter controlling the degree of smoothing.
    min_delta
        Minimum difference in marker expression after normalization to show
        separate trends for the two branches.
    linestyle
        A pair of line styles for the branching trajectory lines, passed to
        :func:`~matplotlib.pyplot.plot`.
    figsize
        width, height
    legend_loc
        The location of the legend.
    legend_fontsize
        Controls the font size of the legend. If the value is numeric the
        size will be the absolute font size in points. String values, are
        relative to the current default font size. This argument is only
        used if `legend_prop` is not specified.
    legend_prop
        The font properties of the legend. If None, the current
        :data:`matplotlib.rcParams` will be used.
    return_fig:
        Return the matplotlib figure.
    {show_save_ax}
    **kargs
        additional arguments passed to :func:`~matplotlib.pyplot.plot`.

    Returns
    -------
    Updates `adata` with the following fields:

    `Trunk_wishbone` : :class:`pandas.DataFrame` (`adata.uns`)
        Computed values before branching
    `Branch1_wishbone` : :class:`pandas.DataFrame` (`adata.uns`)
        Computed values for the first branch
    `Branch2_wishbone` : :class:`pandas.DataFrame` (`adata.uns`)
        Computed values for the second branch.
    """
    # Compute bin locations and bin memberships
    trajectory = adata.obs['trajectory_wishbone'].sort_values()
    wishbone_marker_trajectory(
        adata=adata, no_bins=no_bins, smoothing_factor=smoothing_factor
    )
    weights = adata.uns['weights_wishbone'].copy()

    # Plot marker tsne_res
    xaxis = adata.uns['bins_wishbone']

    # Set up return object
    ret_values = dict()
    ret_values["Trunk"] = pd.DataFrame(
        xaxis[0 : adata.uns['branch_point_bin']], columns=["x"]
    )
    ret_values["Branch1"] = pd.DataFrame(
        xaxis[(adata.uns['branch_point_bin'] - 2) :], columns=["x"]
    )
    ret_values["Branch2"] = pd.DataFrame(
        xaxis[(adata.uns['branch_point_bin'] - 2) :], columns=["x"]
    )

    # Marker colors
    colors = sns.color_palette("Set1", len(markers))
    scaling_factor = 2

    if figsize is None:
        width = 2 * len(markers)
        height = 0.75 * len(markers)
    else:
        width, height = figsize

    fig = plt.figure(figsize=(width, height))
    ax = plt.gca()
    sns.set(context="paper", style="ticks", font_scale=1.5, font="Bitstream Vera Sans")

    for marker, color in zip(markers, colors):

        # Marker expression repeated no bins times
        y = adata[trajectory.index, marker].X
        y = y.toarray()
        y = y.reshape(y.size)
        rep_mark = np.tile(y, [no_bins, 1]).T

        # Normalize y
        y_min = np.percentile(y, 1)
        y = (y - y_min) / (np.percentile(y, 99) - y_min)
        y[y < 0] = 0
        y[y > 1] = 1
        norm_rep_mark = pd.DataFrame(np.tile(y, [no_bins, 1])).T

        if 'branch_wishbone' not in adata.obs:
            # Weight and plot
            vals = (rep_mark * weights) / sum(weights)

            # Normalize
            vals = vals.sum(axis=0)
            vals = vals - np.min(vals)
            vals = vals / np.max(vals)

            # Plot
            plt.plot(xaxis, vals, label=marker, color=color, **kargs)

            # Show errors if specified
            if show_variance:
                # Scale the marks based on y and values to be plotted
                temp = ((norm_rep_mark - vals - np.min(y)) / np.max(y)) ** 2
                # Calculate standard deviations
                wstds = (
                    inner1d(np.asarray(temp).T, np.asarray(weights).T) / weights.sum()
                )

                plt.fill_between(
                    xaxis,
                    vals - scaling_factor * wstds,
                    vals + scaling_factor * wstds,
                    alpha=0.2,
                    color=color,
                )

            # Return values
            ret_values["Trunk"][marker] = vals[0 : adata.uns['branch_point_bin']]
            ret_values["Branch1"][marker] = vals[(adata.uns['branch_point_bin'] - 2) :]
            ret_values["Branch2"][marker] = vals[(adata.uns['branch_point_bin'] - 2) :]

        else:  # Branching trajectory
            rep_mark = pd.DataFrame(
                rep_mark, index=trajectory.index, columns=range(no_bins)
            )

            # Plot trunk first
            weights = adata.uns['weights_wishbone'].copy()

            plot_vals = ((rep_mark * weights) / np.sum(weights)).sum()
            trunk_vals = plot_vals[0 : adata.uns['branch_point_bin']]

            branch_vals = []
            for br in adata.uns['branches_wishbone']:
                # Mute weights of the branch cells and plot
                weights = adata.uns['weights_wishbone'].copy()
                weights.loc[adata.obs.index[adata.obs['branch_wishbone'] == br]] = 0

                plot_vals = ((rep_mark * weights) / np.sum(weights)).sum()
                branch_vals.append(plot_vals[(adata.uns['branch_point_bin'] - 1) :])

            # Min and max
            temp = trunk_vals.append(branch_vals[0]).append(branch_vals[1])
            min_val = np.min(temp)
            max_val = np.max(temp)

            # Plot the trunk
            plot_vals = ((rep_mark * weights) / np.sum(weights)).sum()
            plot_vals = (plot_vals - min_val) / (max_val - min_val)
            plt.plot(
                xaxis[0 : adata.uns['branch_point_bin']],
                plot_vals[0 : adata.uns['branch_point_bin']],
                label=marker,
                color=color,
                **kargs,
            )

            if show_variance:
                # Calculate weighted stds for plotting
                # Scale the marks based on y and values to be plotted
                temp = ((norm_rep_mark - plot_vals - np.min(y)) / np.max(y)) ** 2
                # Calculate standard deviations
                wstds = (
                    inner1d(np.asarray(temp).T, np.asarray(weights).T) / weights.sum()
                )

                # Plot
                plt.fill_between(
                    xaxis[0 : adata.uns['branch_point_bin']],
                    plot_vals[0 : adata.uns['branch_point_bin']]
                    - scaling_factor * wstds[0 : adata.uns['branch_point_bin']],
                    plot_vals[0 : adata.uns['branch_point_bin']]
                    + scaling_factor * wstds[0 : adata.uns['branch_point_bin']],
                    alpha=0.1,
                    color=color,
                )

            # Add values to return values
            ret_values["Trunk"][marker] = plot_vals[0 : adata.uns['branch_point_bin']]

            # Identify markers which need a split
            if (
                sum(
                    abs(pd.Series(branch_vals[0]) - pd.Series(branch_vals[1]))
                    > min_delta
                )
                < 5
            ):
                # Split not necessary, plot the trunk values
                plt.plot(
                    xaxis[(adata.uns['branch_point_bin'] - 1) :],
                    plot_vals[(adata.uns['branch_point_bin'] - 1) :],
                    color=color,
                    **kargs,
                )

                # Add values to return values
                ret_values["Branch1"][marker] = list(
                    plot_vals[(adata.uns['branch_point_bin'] - 2) :]
                )
                ret_values["Branch2"][marker] = list(
                    plot_vals[(adata.uns['branch_point_bin'] - 2) :]
                )

                if show_variance:
                    # Calculate weighted stds for plotting
                    # Scale the marks based on y and values to be plotted
                    temp = ((norm_rep_mark - plot_vals - np.min(y)) / np.max(y)) ** 2
                    wstds = (
                        inner1d(np.asarray(temp).T, np.asarray(weights).T)
                        / weights.sum()
                    )
                    # Plot
                    plt.fill_between(
                        xaxis[(adata.uns['branch_point_bin'] - 1) :],
                        plot_vals[(adata.uns['branch_point_bin'] - 1) :]
                        - scaling_factor * wstds[(adata.uns['branch_point_bin'] - 1) :],
                        plot_vals[(adata.uns['branch_point_bin'] - 1) :]
                        + scaling_factor * wstds[(adata.uns['branch_point_bin'] - 1) :],
                        alpha=0.1,
                        color=color,
                    )
            else:
                if linestyle is None:
                    linestyle = [":", "--"]
                linestyle = pd.Series(linestyle, index=adata.uns['branches_wishbone'])
                # Plot the two branches separately
                for br_ind, br in enumerate(adata.uns['branches_wishbone']):
                    # Mute weights of the branch cells and plot
                    weights = adata.uns['weights_wishbone'].copy()

                    # Smooth weigths
                    smooth_bins = 10
                    if adata.uns['branch_point_bin'] < smooth_bins:
                        smooth_bins = adata.uns['branch_point_bin'] - 1
                    for i in range(smooth_bins):
                        weights.loc[
                            adata.obs.index[adata.obs['branch_wishbone'] == br],
                            adata.uns['branch_point_bin'] + i - smooth_bins,
                        ] *= ((smooth_bins - i) / smooth_bins) * 0.25
                    weights.loc[
                        adata.obs.index[adata.obs['branch_wishbone'] == br],
                        (adata.uns['branch_point_bin']) : weights.shape[1],
                    ] = 0

                    # Calculate values to be plotted
                    plot_vals = ((rep_mark * weights) / np.sum(weights)).sum()
                    plot_vals = (plot_vals - min_val) / (max_val - min_val)
                    plt.plot(
                        xaxis[(adata.uns['branch_point_bin'] - 2) :],
                        plot_vals[(adata.uns['branch_point_bin'] - 2) :],
                        linestyle=linestyle[br],
                        color=color,
                        **kargs,
                    )

                    if show_variance:
                        # Calculate weighted stds for plotting
                        # Scale the marks based on y and values to be plotted
                        temp = (
                            (norm_rep_mark - plot_vals - np.min(y)) / np.max(y)
                        ) ** 2
                        # Calculate standard deviations
                        wstds = (
                            inner1d(np.asarray(temp).T, np.asarray(weights).T)
                            / weights.sum()
                        )

                        # Plot
                        plt.fill_between(
                            xaxis[(adata.uns['branch_point_bin'] - 1) :],
                            plot_vals[(adata.uns['branch_point_bin'] - 1) :]
                            - scaling_factor
                            * wstds[(adata.uns['branch_point_bin'] - 1) :],
                            plot_vals[(adata.uns['branch_point_bin'] - 1) :]
                            + scaling_factor
                            * wstds[(adata.uns['branch_point_bin'] - 1) :],
                            alpha=0.1,
                            color=color,
                        )

                    # Add values to return values
                    ret_values["Branch%d" % (br_ind + 1)][marker] = list(
                        plot_vals[(adata.uns['branch_point_bin'] - 2) :]
                    )

    # Clean up the plotting
    # Clean xlim
    plt.legend(
        loc=legend_loc,
        bbox_to_anchor=(1, 1),
        prop=legend_prop,
        fontsize=legend_fontsize,
    )
    # Annotations
    # Add trajectory as underlay
    cm = matplotlib.cm.Spectral_r
    yval = 0
    plt.scatter(
        trajectory,
        np.repeat(yval - 0.1, len(trajectory)),
        c=trajectory,
        cmap=cm,
        edgecolors="none",
        s=8,
    )
    sns.despine()
    plt.xticks(np.arange(0, 1.1, 0.1))

    # Clean xlim
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.2, 1.1])
    plt.xlabel("Wishbone trajectory")
    plt.ylabel("Normalized expression")

    adata.uns['Trunk_wishbone'] = ret_values['Trunk']
    adata.uns['Branch1_wishbone'] = ret_values['Branch1']
    adata.uns['Branch2_wishbone'] = ret_values['Branch2']

    _utils.savefig_or_show('wishbone_trajectory', show=show, save=save)

    if return_fig:
        return fig
    elif not show:
        return ax

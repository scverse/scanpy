from matplotlib import pyplot as pl
from pandas.api.types import is_categorical_dtype
import numpy as np
from matplotlib import rcParams
from .. import utils
from ...utils import sanitize_anndata, doc_params
from ... import settings
from ..docs import doc_adata_color_etc, doc_edges_arrows, doc_scatter_bulk, doc_show_save_ax
from ... import logging as logg


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def umap(adata, **kwargs):
    """\
    Scatter plot in UMAP basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    return plot_scatter(adata, basis='umap', **kwargs)


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def tsne(adata, **kwargs):
    """\
    Scatter plot in tSNE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    return plot_scatter(adata, basis='tsne', **kwargs)


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def phate(adata, **kwargs):
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
    If `show==False`, a list of `matplotlib.Axis` objects. Every second element
    corresponds to the 'right margin' drawing area for color bars and legends.

    Examples
    --------
    >>> import scanpy.api as sc
    >>> import phate
    >>> data, branches = phate.tree.gen_dla(n_dim=100,
                                            n_branch=20,
                                            branch_length=100)
    >>> data.shape
    (2000, 100)
    >>> adata = sc.AnnData(data)
    >>> adata.obs['branches'] = branches
    >>> sc.tl.phate(adata, k=5, a=20, t=150)
    >>> adata.obsm['X_phate'].shape
    (2000, 2)
    >>> sc.pl.phate(adata,
                    color='branches',
                    color_map='tab20')
    """
    return plot_scatter(adata, basis='phate', **kwargs)


@doc_params(adata_color_etc=doc_adata_color_etc, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def diffmap(adata, **kwargs):
    """\
    Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """

    return plot_scatter(adata, basis='diffmap', **kwargs)


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def draw_graph(adata, layout=None, **kwargs):
    """\
    Scatter plot in graph-drawing basis.

    Parameters
    ----------
    {adata_color_etc}
    layout : {{'fa', 'fr', 'drl', ...}}, optional (default: last computed)
        One of the `draw_graph` layouts, see
        :func:`~scanpy.api.tl.draw_graph`. By default, the last computed layout
        is used.
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    if layout is None:
        layout = str(adata.uns['draw_graph']['params']['layout'])
    basis = 'draw_graph_' + layout
    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError('Did not find {} in adata.obs. Did you compute layout {}?'
                         .format('draw_graph_' + layout, layout))

    return plot_scatter(adata, basis=basis, **kwargs)


@doc_params(adata_color_etc=doc_adata_color_etc, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def pca(adata, **kwargs):
    """\
    Scatter plot in PCA coordinates.

    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    return plot_scatter(adata, basis='pca', **kwargs)


def plot_scatter(adata,
                 color=None,
                 use_raw=True,
                 sort_order=True,
                 edges=False,
                 edges_width=0.1,
                 edges_color='grey',
                 arrows=False,
                 arrows_kwds=None,
                 basis=None,
                 groups=None,
                 components=None,
                 projection='2d',
                 color_map=None,
                 palette=None,
                 right_margin=None,
                 left_margin=None,
                 size=None,
                 frameon=None,
                 legend_fontsize=None,
                 legend_fontweight=None,
                 legend_loc='right margin',
                 panels_per_row=4,
                 title=None,
                 show=None,
                 save=None,
                 ax=None, **kwargs):

    sanitize_anndata(adata)
    if color_map is not None:
        print("instead of color_map use cmap")
        kwargs['cmap'] = color_map
    if size is not None:
        print("instead of size use 's'")
        kwargs['s'] = size

    # TODO
    # are this options really needed?
    if right_margin is not None:
        print("right_margin is currently not used")
    if left_margin is not None:
        print("left_margin is currently not used")

    if projection == '3d':
        from mpl_toolkits.mplot3d import Axes3D
        args_3d = {'projection': '3d'}
    else:
        args_3d = {}

    ####
    # get the points position and the components list (only if components is not 'None)
    data_points, components_list = _get_data_points(adata, basis, projection, components)

    ###
    # setup layout. Most of the code is for the case when multiple plots are required
    # 'color' is a list of names that want to be plotted. Eg. ['Gene1', 'louvain', 'Gene2']
    if isinstance(color, list) and len(color) > 1:
        if ax is not None:
            raise ValueError("When plotting multiple panels (each for a given value of 'color' "
                             "a given ax can not be used")

        multi_panel = True
        # each plot needs to be its own panel
        from matplotlib import gridspec
        # set up the figure
        n_panels_x = panels_per_row
        n_panels_y = np.ceil(len(color) / n_panels_x).astype(int)
        # each panel will have the size of rcParams['figure.figsize']
        fig = pl.figure(figsize=(n_panels_x * rcParams['figure.figsize'][0],
                                 n_panels_y * rcParams['figure.figsize'][1]))
        left = 0.2 / n_panels_x
        bottom = 0.13 / n_panels_y
        gs = gridspec.GridSpec(nrows=n_panels_y,
                               ncols=n_panels_x,
                               left=left,
                               right=1-(n_panels_x-1)*left-0.01/n_panels_x,
                               bottom=bottom,
                               top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
                               hspace=0.2,
                               wspace=0.2)
    else:
        # this case handles color='variable' and color=['variable'], which are the same
        if isinstance(color, str) or color is None:
            color = [color]
        multi_panel = False
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111, **args_3d)

    ###
    # make the plots
    axs = []
    for count, value_to_plot in enumerate(color):

        ###
        # when plotting, the color of the dots is determined for each plot
        # the data is either categorical or continuous and the data could be in
        # 'obs' or in 'var'
        categorical = False
        if value_to_plot is None:
            color_vector = 'lightgray'
        # check if value to plot is in obs
        elif value_to_plot in adata.obs.columns:
            if is_categorical_dtype(adata.obs[value_to_plot]):
                categorical = True

                try:
                    adata.uns[value_to_plot + "_colors"]
                except KeyError:
                    utils.add_colors_for_categorical_sample_annotation(adata, value_to_plot)

                if palette:
                    utils.add_colors_for_categorical_sample_annotation(adata, value_to_plot,
                                                                       palette=palette, force_update_colors=True)

                # for categorical data, colors should be
                # stored in adata.uns
                # Obtain color vector by converting every category
                # into its respective color

                color_vector = [adata.uns[value_to_plot + '_colors'][x] for x in adata.obs[value_to_plot].cat.codes]
                if groups is not None:
                    if isinstance(groups, str):
                        groups = [groups]
                    color_vector = np.array(color_vector, dtype='<U15')
                    # set color to 'light gray' for all values
                    # that are not in the groups
                    color_vector[~adata.obs[value_to_plot].isin(groups)] = "lightgray"
            else:
                color_vector = adata.obs[value_to_plot]
        # check if value to plot is in var
        elif use_raw is False and value_to_plot in adata.var_names:
            # TODO
            # can this be done in a different, probably in a faster way?
            color_vector = adata[:,value_to_plot].X

        elif use_raw is True and value_to_plot in adata.raw.var_names:
            color_vector = adata.raw[:,value_to_plot].X
        else:
            raise ValueError("Given 'color': {} is not a valid observation "
                             "or var. Valid observations are: {}".format(value_to_plot, adata.obs.columns))

        # check if higher value points should be plot on top
        if value_to_plot is not None and categorical is False and sort_order is True:
            order = np.argsort(color_vector)
            color_vector = color_vector[order]
            data_points = data_points[order, :]

        # if plotting multiple panels, get the ax from the grid spec
        # else use the ax value (either user given or created previously)
        if multi_panel is True:
            ax = pl.subplot(gs[count], **args_3d)
            axs.append(ax)
        if frameon is False:
            ax.axis('off')
        if title is None and value_to_plot is not None:
            ax.set_title(value_to_plot)
        else:
            ax.set_title(title)

        if 's' not in kwargs:
            kwargs['s'] = 120000 / data_points.shape[0]

        # make the scatter plot
        cax= ax.scatter(data_points[:, 0], data_points[:, 1],
                        marker=".", c=color_vector, rasterized=settings._vector_friendly,
                        **kwargs)

        # remove y and x ticks
        ax.set_yticks([])
        ax.set_xticks([])
        if projection == '3d':
            ax.set_zticks([])

        # set default axis_labels
        name = _basis2name(basis)
        if components is not None:
            axis_labels = [name + str(x + 1) for x in components_list]
        elif projection == '3d':
            axis_labels = [name + str(x + 1) for x in range(3)]

        else:
            axis_labels = [name + str(x + 1) for x in range(2)]

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        if projection == '3d':
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[2], labelpad=-7)
        ax.autoscale_view()

        if edges:
            plot_edges(ax, adata, basis, edges_width, edges_color)
        if arrows:
            plot_arrows(ax, adata, basis, arrows_kwds)

        if value_to_plot is None:
            # if only dots were plotted without an associated value
            # there is not need to plot a legend or a colorbar
            continue

        _add_legend_or_colorbar(adata, ax, cax, categorical, value_to_plot, legend_loc,
                                data_points, legend_fontweight, legend_fontsize, groups,
                                multi_panel)

    axs = axs if multi_panel else ax
    utils.savefig_or_show(basis, show=show, save=save)
    if show == False:
        return axs


def _get_data_points(adata, basis, projection, components):
    """
    Returns the data points corresponding to the selected basis, projection and/or components
    """
    n_dims = 2
    if projection == '3d':
        # check if the data has a third dimension
        if adata.obsm['X_' + basis].shape[1] == 2:
            if settings._low_resolution_warning:
                logg.warn('Selected projections is "3d" but only two dimensions '
                          'are available. Only these two dimensions will be plotted')
        else:
            n_dims = 3

    if components is not None:
        # components should be a string of coma separated integers
        if isinstance(components, str):
            components_list = [int(x.strip()) - 1 for x in components.split(',')]
        elif isinstance(components, list):
            components_list = [int(x) - 1 for x in components]
        else:
            raise ValueError("Given components: '{}' are not valid. Please check. "
                             "A valid example is `components='2,3'`")
        # check if the components are present in the data
        try:
            data_points = adata.obsm['X_' + basis][:, components_list]
        except:
            raise ValueError("Given components: '{}' are not valid. Please check. "
                             "A valid example is `components='2,3'`")

    else:
        data_points = adata.obsm['X_' + basis][:, :n_dims]
        components_list = None
    return data_points, components_list


def _add_legend_or_colorbar(adata, ax, cax, categorical, value_to_plot, legend_loc,
                            scatter_array, legend_fontweight, legend_fontsize, groups,
                            multi_panel):
    """
    Adds a color bar or a legend to the given ax. A legend is added when the
    data is categorical and a color bar is added when a continuous value was used.

    """
    # add legends or colorbars
    if categorical is True:
        # add legend to figure
        categories = list(adata.obs[value_to_plot].cat.categories)
        colors = adata.uns[value_to_plot + '_colors']

        if groups is not None:
            # only label groups with the respective color
            colors = [colors[categories.index(x)] for x in groups]
            categories = groups

        if legend_loc == 'right margin':
            for idx, label in enumerate(categories):
                color = colors[idx]
                # use empty scatter to set labels
                ax.scatter([], [], c=color, label=label)
            if multi_panel is True:
                # Shrink current axis by 20% to fit legend
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(
                frameon=False, loc='center left',
                bbox_to_anchor=(1, 0.5),
                ncol=(1 if len(categories) <= 14
                      else 2 if len(categories) <= 30 else 3),
                fontsize=legend_fontsize)

        if legend_loc == 'on data':
            # identify centroids to put labels
            for label in categories:
                _scatter = scatter_array[adata.obs[value_to_plot] == label, :]
                x_pos, y_pos = np.median(_scatter, axis=0)

                ax.text(x_pos, y_pos, label,
                        weight=legend_fontweight,
                        verticalalignment='center',
                        horizontalalignment='center',
                        fontsize=legend_fontsize)
    else:
        # add colorbar to figure
        pl.colorbar(cax, ax=ax)


def _basis2name(basis):
    """
    converts the 'basis' into the proper name.
    """

    component_name = (
        'DC' if basis == 'diffmap'
        else 'tSNE' if basis == 'tsne'
        else 'UMAP' if basis == 'umap'
        else 'PC' if basis == 'pca'
        else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis
        else basis)
    return component_name


def plot_edges(ax, adata, basis, edges_width, edges_color):
    import networkx as nx

    if 'neighbors' not in adata.uns:
        raise ValueError('`edges=True` requires `pp.neighbors` to be run before.')
    g = nx.Graph(adata.uns['neighbors']['connectivities'])
    edge_collection = nx.draw_networkx_edges(
        g, adata.obsm['X_' + basis],
        ax=ax, width=edges_width, edge_color=edges_color)
    edge_collection.set_zorder(-2)
    edge_collection.set_rasterized(settings._vector_friendly)


def plot_arrows(ax, adata, basis, arrows_kwds=None):
    if 'Delta_' + basis not in adata.obsm.keys():
        raise ValueError('`arrows=True` requires \'Delta_\' + basis from velocyto.')
    X = adata.obsm['X_' + basis]
    V = adata.obsm['Delta_' + basis]
    quiver_kwds = arrows_kwds if arrows_kwds is not None else {}
    ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], **quiver_kwds,
              rasterized=settings._vector_friendly)

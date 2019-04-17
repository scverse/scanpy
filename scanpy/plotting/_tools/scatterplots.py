from collections import abc
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple

import numpy as np
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pandas.api.types import is_categorical_dtype
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import is_color_like, Colormap

from .. import _utils as utils
from .._docs import doc_adata_color_etc, doc_edges_arrows, doc_scatter_bulk, doc_show_save_ax
from ..._settings import settings
from ...utils import sanitize_anndata, doc_params
from ... import logging as logg


def plot_scatter(
    adata: AnnData,
    basis: str,
    *,
    color: Union[str, Sequence[str], None] = None,
    gene_symbols: Optional[str] = None,
    use_raw: Optional[bool] = None,
    sort_order: bool = True,
    edges: bool = False,
    edges_width: float = 0.1,
    edges_color: Union[str, Sequence[float], Sequence[str]] = 'grey',
    arrows: bool = False,
    arrows_kwds: Optional[Mapping[str, Any]] = None,
    groups: Optional[str] = None,
    components: Union[str, Sequence[str]] = None,
    layer: Optional[str] = None,
    projection: str = '2d',
    color_map: Union[Colormap, str, None] = None,
    palette: Union[str, Sequence[str], Cycler, None] = None,
    size: Optional[float] = None,
    frameon: Optional[bool] = None,
    legend_fontsize: Optional[int] = None,
    legend_fontweight: str = 'bold',
    legend_loc: str = 'right margin',
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: Optional[float] = None,
    title: Union[str, Sequence[str], None] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Axes] = None,
    return_fig: Optional[bool] = None,
    **kwargs
) -> Union[Figure, Axes, None]:
    sanitize_anndata(adata)
    if color_map is not None:
        kwargs['cmap'] = color_map
    if size is not None:
        kwargs['s'] = size
    if 'edgecolor' not in kwargs:
        # by default turn off edge color. Otherwise, for
        # very small sizes the edge will not reduce its size
        # (https://github.com/theislab/scanpy/issues/293)
        kwargs['edgecolor'] = 'none'

    if projection == '3d':
        from mpl_toolkits.mplot3d import Axes3D
        args_3d = {'projection': '3d'}
    else:
        args_3d = {}

    if use_raw is None:
        # check if adata.raw is set
        use_raw = adata.raw is not None

    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    if adata.raw is None and use_raw:
        raise ValueError(
            "`use_raw` is set to True but AnnData object does not have raw. "
            "Please check."
        )
    # turn color into a python list
    color = [color] if isinstance(color, str) or color is None else list(color)
    if title is not None:
        # turn title into a python list if not None
        title = [title] if isinstance(title, str) else list(title)

    ####
    # get the points position and the components list (only if components is not 'None)
    data_points, components_list = _get_data_points(adata, basis, projection, components)

    ###
    # setup layout. Most of the code is for the case when multiple plots are required
    # 'color' is a list of names that want to be plotted. Eg. ['Gene1', 'louvain', 'Gene2'].
    # component_list is a list of components [[0,1], [1,2]]
    if (isinstance(color, abc.Sequence) and len(color) > 1) or len(components_list) > 1:
        if ax is not None:
            raise ValueError(
                "When plotting multiple panels (each for a given value of 'color') "
                "a given ax can not be used"
            )
        if len(components_list) == 0:
            components_list = [None]

        multi_panel = True
        # each plot needs to be its own panel
        from matplotlib import gridspec
        # set up the figure
        num_panels = len(color) * len(components_list)
        n_panels_x = min(ncols, num_panels)
        n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
        # each panel will have the size of rcParams['figure.figsize']
        fig = pl.figure(figsize=(n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
                                 n_panels_y * rcParams['figure.figsize'][1]))
        left = 0.2 / n_panels_x
        bottom = 0.13 / n_panels_y
        gs = gridspec.GridSpec(
            nrows=n_panels_y, ncols=n_panels_x,
            left=left, right=1-(n_panels_x-1)*left-0.01/n_panels_x,
            bottom=bottom, top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
            hspace=hspace, wspace=wspace,
        )
    else:
        if len(components_list) == 0:
            components_list = [None]
        multi_panel = False
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111, **args_3d)

    ###
    # make the plots
    axs = []
    import itertools
    idx_components = range(len(components_list))

    # use itertools.product to make a plot for each color and for each component
    # For example if color=[gene1, gene2] and components=['1,2, '2,3'].
    # The plots are: [color=gene1, components=[1,2], color=gene1, components=[2,3],
    #                 color=gene2, components = [1, 2], color=gene2, components=[2,3]]
    for count, (value_to_plot, component_idx) in enumerate(itertools.product(color, idx_components)):
        color_vector, categorical = _get_color_values(
            adata, value_to_plot, layer=layer,
            groups=groups, palette=palette,
            use_raw=use_raw, gene_symbols=gene_symbols,
        )

        # check if higher value points should be plot on top
        if sort_order is True and value_to_plot is not None and categorical is False:
            order = np.argsort(color_vector)
            color_vector = color_vector[order]
            _data_points = data_points[component_idx][order, :]

            # check if 'size' is given (stored in kwargs['s']
            # and reorder it.
            import pandas.core.series
            if 's' in kwargs and kwargs['s'] is not None \
                and isinstance(kwargs['s'],(list, pandas.core.series.Series, np.ndarray)) \
                and len(kwargs['s']) == len(color_vector):
                kwargs['s'] = np.array(kwargs['s'])[order]
        else:
            _data_points = data_points[component_idx]

        # if plotting multiple panels, get the ax from the grid spec
        # else use the ax value (either user given or created previously)
        if multi_panel is True:
            ax = pl.subplot(gs[count], **args_3d)
            axs.append(ax)
        if not (settings._frameon if frameon is None else frameon):
            ax.axis('off')
        if title is None:
            if value_to_plot is not None:
                ax.set_title(value_to_plot)
            else:
                ax.set_title('')
        else:
            try:
                ax.set_title(title[count])
            except IndexError:
                logg.warn("The title list is shorter than the number of panels. Using 'color' value instead for"
                          "some plots.")
                ax.set_title(value_to_plot)

        if 's' not in kwargs:
            kwargs['s'] = 120000 / _data_points.shape[0]

        # make the scatter plot
        if projection == '3d':
            cax = ax.scatter(
                _data_points[:, 0], _data_points[:, 1], _data_points[:, 2],
                marker=".", c=color_vector, rasterized=settings._vector_friendly,
                **kwargs,
            )
        else:
            cax = ax.scatter(
                _data_points[:, 0], _data_points[:, 1],
                marker=".", c=color_vector, rasterized=settings._vector_friendly,
                **kwargs,
            )

        # remove y and x ticks
        ax.set_yticks([])
        ax.set_xticks([])
        if projection == '3d':
            ax.set_zticks([])

        # set default axis_labels
        name = _basis2name(basis)
        if components is not None:
            axis_labels = [name + str(x + 1) for x in components_list[component_idx]]
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
            utils.plot_edges(ax, adata, basis, edges_width, edges_color)
        if arrows:
            utils.plot_arrows(ax, adata, basis, arrows_kwds)

        if value_to_plot is None:
            # if only dots were plotted without an associated value
            # there is not need to plot a legend or a colorbar
            continue

        _add_legend_or_colorbar(
            adata, ax, cax, categorical, value_to_plot, legend_loc,
            _data_points, legend_fontweight, legend_fontsize, groups, multi_panel,
        )

    if return_fig is True:
        return fig
    axs = axs if multi_panel else ax
    utils.savefig_or_show(basis, show=show, save=save)
    if show is False:
        return axs


def _wraps_plot_scatter(wrapper):
    annots_orig = {
        k: v for k, v in wrapper.__annotations__.items()
        if k not in {'adata', 'kwargs'}
    }
    annots_scatter = {
        k: v for k, v in plot_scatter.__annotations__.items()
        if k != 'basis'
    }
    wrapper.__annotations__ = {**annots_scatter, **annots_orig}
    wrapper.__wrapped__ = plot_scatter
    return wrapper


# API


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def umap(adata, **kwargs) -> Union[Axes, List[Axes], None]:
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
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return plot_scatter(adata, 'umap', **kwargs)


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def tsne(adata, **kwargs) -> Union[Axes, List[Axes], None]:
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
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return plot_scatter(adata, 'tsne', **kwargs)


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
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
    If `show==False`, a list of :class:`~matplotlib.axes.Axes` objects. Every second element
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
    return plot_scatter(adata, 'phate', **kwargs)


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def diffmap(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return plot_scatter(adata, 'diffmap', **kwargs)


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def draw_graph(adata, layout=None, **kwargs) -> Union[Axes, List[Axes], None]:
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
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    if layout is None:
        layout = str(adata.uns['draw_graph']['params']['layout'])
    basis = 'draw_graph_' + layout
    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError('Did not find {} in adata.obs. Did you compute layout {}?'
                         .format('draw_graph_' + layout, layout))

    return plot_scatter(adata, basis, **kwargs)


@_wraps_plot_scatter
@doc_params(adata_color_etc=doc_adata_color_etc, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def pca(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in PCA coordinates.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return plot_scatter(adata, 'pca', **kwargs)


# Helpers


def _get_data_points(adata, basis, projection, components) -> Tuple[List[np.ndarray], List[Tuple[int, int]]]:
    """
    Returns the data points corresponding to the selected basis, projection and/or components.

    Because multiple components are given (eg components=['1,2', '2,3'] the
    returned data are lists, containing each of the components. When only one component is plotted
    the list length is 1.

    Returns
    -------
    data_points : list
        Each entry is a numpy array containing the data points
    components : list
        The cleaned list of components. Eg. [(0,1)] or [(0,1), (1,2)]
        for components = [1,2] and components=['1,2', '2,3'] respectively
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

    if components == 'all':
        components = ['{},{}'.format(*((i, i+1) if i % 2 == 1 else (i+1, i)))
            for i in range(1, adata.obsm['X_{}'.format(basis)].shape[1])]

    components_list = []
    offset = 0
    if basis == 'diffmap': offset = 1
    if components is not None:
        # components have different formats, either a list with integers, a string
        # or a list of strings.

        if isinstance(components, str):
            # eg: components='1,2'
            components_list.append(tuple(int(x.strip()) - 1 + offset for x in components.split(',')))

        elif isinstance(components, abc.Sequence):
            if isinstance(components[0], int):
                # components=[1,2]
                components_list.append(tuple(int(x) - 1 + offset for x in components))
            else:
                # in this case, the components are str
                # eg: components=['1,2'] or components=['1,2', '2,3]
                # More than one component can be given and is stored
                # as a new item of components_list
                for comp in components:
                    components_list.append(tuple(int(x.strip()) - 1 + offset for x in comp.split(',')))

        else:
            raise ValueError("Given components: '{}' are not valid. Please check. "
                             "A valid example is `components='2,3'`")
        # check if the components are present in the data
        try:
            data_points = []
            for comp in components_list:
                data_points.append(adata.obsm['X_' + basis][:, comp])
        except:
            raise ValueError("Given components: '{}' are not valid. Please check. "
                             "A valid example is `components='2,3'`")

        if basis == 'diffmap':
            # remove the offset added in the case of diffmap, such that
            # plot_scatter can print the labels correctly.
            components_list = [tuple(number-1 for number in comp) for comp in components_list]
    else:
        data_points = [adata.obsm['X_' + basis][:, offset:offset+n_dims]]
        components_list = []
    return data_points, components_list


def _add_legend_or_colorbar(adata, ax, cax, categorical, value_to_plot, legend_loc,
                            scatter_array, legend_fontweight, legend_fontsize,
                            groups, multi_panel):
    """
    Adds a color bar or a legend to the given ax. A legend is added when the
    data is categorical and a color bar is added when a continuous value was used.

    """
    # add legends or colorbars
    if categorical is True:
        # add legend to figure
        categories = list(adata.obs[value_to_plot].cat.categories)
        colors = adata.uns[value_to_plot + '_colors']

        if multi_panel is True:
            # Shrink current axis by 10% to fit legend and match
            # size of plots that are not categorical
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.91, box.height])

        if groups is not None:
            # only label groups with the respective color
            colors = [colors[categories.index(x)] for x in groups]
            categories = groups

        if legend_loc == 'right margin':
            for idx, label in enumerate(categories):
                color = colors[idx]
                # use empty scatter to set labels
                ax.scatter([], [], c=color, label=label)
            ax.legend(
                frameon=False, loc='center left',
                bbox_to_anchor=(1, 0.5),
                ncol=(1 if len(categories) <= 14
                      else 2 if len(categories) <= 30 else 3),
                fontsize=legend_fontsize)

        if legend_loc == 'on data':
            # identify centroids to put labels
            all_pos = np.zeros((len(categories), 2))
            for ilabel, label in enumerate(categories):
                _scatter = scatter_array[adata.obs[value_to_plot] == label, :]
                x_pos, y_pos = np.median(_scatter, axis=0)

                ax.text(x_pos, y_pos, label,
                        weight=legend_fontweight,
                        verticalalignment='center',
                        horizontalalignment='center',
                        fontsize=legend_fontsize)

                all_pos[ilabel] = [x_pos, y_pos]
            # this is temporary storage for access by other tools
            utils._tmp_cluster_pos = all_pos
    else:
        # add colorbar to figure
        pl.colorbar(cax, ax=ax, pad=0.01, fraction=0.08, aspect=30)


def _set_colors_for_categorical_obs(adata, value_to_plot, palette):
    """
    Sets the adata.uns[value_to_plot + '_colors'] according to the given palette

    Parameters
    ----------
    adata
        annData object
    value_to_plot
        name of a valid categorical observation
    palette
        Palette should be either a valid :func:`~matplotlib.pyplot.colormaps` string,
        a list of colors (in a format that can be understood by matplotlib,
        eg. RGB, RGBS, hex, or a cycler object with key='color'

    Returns
    -------
    None
    """
    from matplotlib.colors import to_hex
    from cycler import Cycler, cycler

    categories = adata.obs[value_to_plot].cat.categories
    # check is palette is a valid matplotlib colormap
    if isinstance(palette, str) and palette in pl.colormaps():
        # this creates a palette from a colormap. E.g. 'Accent, Dark2, tab20'
        cmap = pl.get_cmap(palette)
        colors_list = [to_hex(x) for x in cmap(np.linspace(0, 1, len(categories)))]

    else:
        # check if palette is a list and convert it to a cycler, thus
        # it doesnt matter if the list is shorter than the categories length:
        if isinstance(palette, abc.Sequence):
            if len(palette) < len(categories):
                logg.warn("Length of palette colors is smaller than the number of "
                          "categories (palette length: {}, categories length: {}. "
                          "Some categories will have the same color."
                          .format(len(palette), len(categories)))
            # check that colors are valid
            _color_list = []
            for color in palette:
                if not is_color_like(color):
                    # check if the color is a valid R color and translate it
                    # to a valid hex color value
                    if color in utils.additional_colors:
                        color = utils.additional_colors[color]
                    else:
                        raise ValueError("The following color value of the given palette is not valid: {}".format(color))
                _color_list.append(color)

            palette = cycler(color=_color_list)
        if not isinstance(palette, Cycler):
            raise ValueError("Please check that the value of 'palette' is a "
                             "valid matplotlib colormap string (eg. Set2), a "
                             "list of color names or a cycler with a 'color' key.")
        if 'color' not in palette.keys:
            raise ValueError("Please set the palette key 'color'.")

        cc = palette()
        colors_list = [to_hex(next(cc)['color']) for x in range(len(categories))]

    adata.uns[value_to_plot + '_colors'] = colors_list


def _set_default_colors_for_categorical_obs(adata, value_to_plot):
    """
    Sets the adata.uns[value_to_plot + '_colors'] using default color palettes

    Parameters
    ----------
    adata : annData object
    value_to_plot : name of a valid categorical observation

    Returns
    -------
    None
    """
    from .. import palettes

    categories = adata.obs[value_to_plot].cat.categories
    length = len(categories)

    # check if default matplotlib palette has enough colors
    if len(rcParams['axes.prop_cycle'].by_key()['color']) >= length:
        cc = rcParams['axes.prop_cycle']()
        palette = [next(cc)['color'] for _ in range(length)]

    else:
        if length <= 28:
            palette = palettes.default_26
        elif length <= len(palettes.default_64):  # 103 colors
            palette = palettes.default_64
        else:
            palette = ['grey' for i in range(length)]
            logg.info('the obs value: "{}" has more than 103 categories. Uniform '
                      '\'grey\' color will be used for all categories.')

    adata.uns[value_to_plot + '_colors'] = palette[:length]


def _get_color_values(adata, value_to_plot, groups=None, palette=None, use_raw=False,
                      gene_symbols=None, layer=None):
    """
    Returns the value or color associated to each data point.
    For categorical data, the return value is list of colors taken
    from the category palette or from the given `palette` value.

    For non-categorical data, the values are returned
    """

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

            if palette:
                # use category colors base on given palette
                _set_colors_for_categorical_obs(adata, value_to_plot, palette)
            else:
                if value_to_plot + '_colors' not in adata.uns or \
                    len(adata.uns[value_to_plot + '_colors']) < len(adata.obs[value_to_plot].cat.categories):
                    #  set a default palette in case that no colors or few colors are found
                    _set_default_colors_for_categorical_obs(adata, value_to_plot)
                else:
                    # check that the colors in 'uns' are valid
                    _palette = []
                    for color in adata.uns[value_to_plot + '_colors']:
                        if not is_color_like(color):
                            # check if the color is a valid R color and translate it
                            # to a valid hex color value
                            if color in utils.additional_colors:
                                color = utils.additional_colors[color]
                            else:
                                logg.warn("The following color value found in adata.uns['{}'] "
                                          " is not valid: '{}'. Default colors are used.".format(value_to_plot + '_colors', color))
                                _set_default_colors_for_categorical_obs(adata, value_to_plot)
                                _palette = None
                                break
                        _palette.append(color)
                    if _palette is not None:
                        adata.uns[value_to_plot + '_colors'] = _palette
            # for categorical data, colors should be
            # stored in adata.uns[value_to_plot + '_colors']
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
            color_vector = adata.obs[value_to_plot].values
    # when value_to_plot is not in adata.obs
    else:
        if gene_symbols is not None and gene_symbols in adata.var.columns:
            if value_to_plot not in adata.var[gene_symbols].values:
                logg.error("Gene symbol {!r} not found in given gene_symbols "
                           "column: {!r}".format(value_to_plot, gene_symbols))
                return
            value_to_plot = adata.var[adata.var[gene_symbols] == value_to_plot].index[0]
        if layer is not None and value_to_plot in adata.var_names:
            if layer not in adata.layers.keys():
                raise KeyError('Selected layer: {} is not in the layers list. The list of '
                               'valid layers is: {}'.format(layer, adata.layers.keys()))
            color_vector = adata[:, value_to_plot].layers[layer]
        elif use_raw and value_to_plot in adata.raw.var_names:
            color_vector = adata.raw[:, value_to_plot].X
        elif value_to_plot in adata.var_names:
            color_vector = adata[:, value_to_plot].X
        else:
            raise ValueError("The passed `color` {} is not a valid observation annotation "
                             "or variable name. Valid observation annotation keys are: {}"
                             .format(value_to_plot, adata.obs.columns))

    return color_vector, categorical


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

"""\
Shared docstrings for plotting function parameters.
"""


doc_adata_color_etc = """\
adata
    Annotated data matrix.
color
    Keys for annotations of observations/cells or variables/genes, e.g.,
    `'ann1'` or `['ann1', 'ann2']`.
gene_symbols
    Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
    refer to the index column of the `.var` DataFrame. Setting this option allows
    alternative names to be used.
use_raw
    Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,
    defaults to `True` if `layer` isn't provided and `adata.raw` is present.
layer
    Name of the AnnData object layer that wants to be plotted. By default
    adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.
    If `layer` is set to a valid layer name, then the layer is plotted. `layer`
    takes precedence over `use_raw`.\
"""

doc_edges_arrows = """\
edges
    Show edges.
edges_width
    Width of edges.
edges_color
    Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.
arrows
    Show arrows (requires to run :func:`scvelo.tl.velocity_embedding` before).
    Deprecated in favor of :func:`scvelo.pl.velocity_embedding` and friends.
arrows_kwds
    Passed to :meth:`~matplotlib.axes.Axes.quiver`\
"""

# Docs for pl.scatter
doc_scatter_basic = """\
sort_order
    For continuous annotations used as color parameter, plot data points
    with higher values on top of others.
groups
    Restrict to a few categories in categorical observation annotation.
    The default is not to restrict to any groups.
components
    For instance, `['1,2', '2,3']`. To plot all available components use
    `components='all'`.
projection
    Projection of plot (default: `'2d'`).
legend_loc
    Location of legend, either `'on data'`, `'right margin'` or a valid keyword
    for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
legend_fontsize
    Numeric size in pt or string describing the size.
    See :meth:`~matplotlib.text.Text.set_fontsize`.
legend_fontweight
    Legend font weight. A numeric value in range 0-1000 or a string.
    Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
    See :meth:`~matplotlib.text.Text.set_fontweight`.
legend_fontoutline
    Line width of the legend font outline in pt. Draws a white outline using
    the path effect :class:`~matplotlib.patheffects.withStroke`.
size
    Point size. If `None`, is automatically computed as 120000 / n_cells.
    Can be a sequence containing the size for each cell. The order should be
    the same as in adata.obs. If `img_key` not `None`, size is the scaling factor
    for the spot size.
color_map
    Color map to use for continous variables. Can be a name or a
    :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
    or `mpl.cm.cividis`), see :func:`~matplotlib.cm.get_cmap`.
    If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
    The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
palette
    Colors to use for plotting categorical annotation groups.
    The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
    (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, or a sequence of
    matplotlib colors like `['red', '#ccdd11', (0.1, 0.2, 1)]`
    (see :func:`~matplotlib.colors.is_color_like`).
    If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
    variable already has colors stored in `adata.uns["{var}_colors"]`.
    If provided, values of `adata.uns["{var}_colors"]` will be set.
frameon
    Draw a frame around the scatter plot. Defaults to value set in
    :func:`~scanpy.set_figure_params`, defaults to `True`.
title
    Provide title for panels either as string or list of strings,
    e.g. `['title1', 'title2', ...]`.
img_key
    Key for image data, stored in `adata.uns`.
crop_coord
    Coordinates to use for cropping the image (left, right, top, bottom).
alpha_img
    Alpha value for image.
bw
    Plot Image in gray scale.\
"""

doc_vminmax = """\
vmin
    Minimum value to plot. Values smaller than vmin are plotted with the same color as vmin.
    vmin can be a number, a string, a function or `None`. If vmin is a string and has the format `pN`,
    this is interpreted as a vmin=percentile(N). For example vmin='p1.5' is interpreted as
    the 1.5 percentile. If vmin is function, then vmin is interpreted as the return value
    of the function over the list of values to plot. For example to set vmin tp the mean of
    the values to plot, `def my_vmin(values): return np.mean(values)` and then
    set `vmin=my_vmin`. If vmin is None (default) an automatic minimum value is used
    as defined by matplotlib `scatter` function. When making multiple plots, vmin can
    be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
vmax
    Maximum value to plot. The format is the same as for `vmin`\
"""

doc_outline = """\
add_outline
    If set to True, this will add a thin border around groups of dots. In some situations
    this can enhance the aesthetics of the resulting image
outline_color
    Tuple with two valid color names used to adjust the add_outline. The first color is the
    border color (default: black), while the second color is a gap color between the
    border color and the scatter dot (default: white).
outline_width
    Tuple with two width numbers used to adjust the outline. The first value is the width
    of the border color as a fraction of the scatter dot size (default: 0.3). The second value is
    width of the gap color (default: 0.05).\
"""

doc_panels = """\
ncols
    Number of panels per row.
wspace
    Adjust the width of the space between multiple panels.
hspace
    Adjust the height of the space between multiple panels.
return_fig
    Return the matplotlib figure.\
"""

# Docs for pl.pca, pl.tsne, … (everything in _tools.scatterplots)
doc_scatter_embedding = f"""\
{doc_scatter_basic}
{doc_vminmax}
{doc_outline}
{doc_panels}
kwargs
    Arguments to pass to :func:`matplotlib.pyplot.scatter`,
    for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).\
"""

doc_show_save_ax = """\
show
     Show the plot, do not return axis.
save
    If `True` or a `str`, save the figure.
    A string is appended to the default filename.
    Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
ax
    A matplotlib axes object. Only works if plotting a single component.\
"""


doc_common_plot_args = """\
adata
    Annotated data matrix.
var_names
    `var_names` should be a valid subset of `adata.var_names`.
    If `var_names` is a mapping, then the key is used as label
    to group the values (see `var_group_labels`). The mapping values
    should be sequences of valid `adata.var_names`. In this
    case either coloring or 'brackets' are used for the grouping
    of var names depending on the plot. When `var_names` is a mapping,
    then the `var_group_labels` and `var_group_positions` are set.
groupby
    The key of the observation grouping to consider.
use_raw
    Use `raw` attribute of `adata` if present.
log
    Plot on logarithmic axis.
num_categories
    Only used if groupby observation is not categorical. This value
    determines the number of groups into which the groupby observation
    should be subdivided.
figsize
    Figure size when `multi_panel=True`.
    Otherwise the `rcParam['figure.figsize]` value is used.
    Format is (width, height)
dendrogram
    If True or a valid dendrogram key, a dendrogram based on the hierarchical
    clustering between the `groupby` categories is added.
    The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
    If `tl.dendrogram` has not been called previously the function is called
    with default parameters.
gene_symbols
    Column name in `.var` DataFrame that stores gene symbols.
    By default `var_names` refer to the index column of the `.var` DataFrame.
    Setting this option allows alternative names to be used.
var_group_positions
    Use this parameter to highlight groups of `var_names`.
    This will draw a 'bracket' or a color block between the given start and end
    positions. If the parameter `var_group_labels` is set, the corresponding
    labels are added on top/left. E.g. `var_group_positions=[(4,10)]`
    will add a bracket between the fourth `var_name` and the tenth `var_name`.
    By giving more positions, more brackets/color blocks are drawn.
var_group_labels
    Labels for each of the `var_group_positions` that want to be highlighted.
var_group_rotation
    Label rotation degrees.
    By default, labels larger than 4 characters are rotated 90 degrees.
layer
    Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
    If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
    then the layer is plotted. `layer` takes precedence over `use_raw`.\
"""

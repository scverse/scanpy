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
neighbors_key
    Where to look for neighbors connectivities.
    If not specified, this looks .obsp['connectivities'] for connectivities
    (default storage place for pp.neighbors).
    If specified, this looks
    .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
arrows
    Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).
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
dimensions
    0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
    Unlike `components`, this argument is used in the same way as `colors`, e.g. is
    used to specify a single plot at a time. Will eventually replace the components
    argument.
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
colorbar_loc
    Where to place the colorbar for continous variables. If `None`, no colorbar
    is added.
size
    Point size. If `None`, is automatically computed as 120000 / n_cells.
    Can be a sequence containing the size for each cell. The order should be
    the same as in adata.obs.
color_map
    Color map to use for continous variables. Can be a name or a
    :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
    or `mpl.cm.cividis`), see :func:`~matplotlib.cm.get_cmap`.
    If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
    The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
palette
    Colors to use for plotting categorical annotation groups.
    The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
    (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
    categories to colors, or a sequence of colors. Colors must be valid to
    matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
    If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
    variable already has colors stored in `adata.uns["{var}_colors"]`.
    If provided, values of `adata.uns["{var}_colors"]` will be set.
na_color
    Color to use for null or masked values. Can be anything matplotlib accepts as a
    color. Used for all points if `color=None`.
na_in_legend
    If there are missing values, whether they get an entry in the legend. Currently
    only implemented for categorical legends.
frameon
    Draw a frame around the scatter plot. Defaults to value set in
    :func:`~scanpy.set_figure_params`, defaults to `True`.
title
    Provide title for panels either as string or list of strings,
    e.g. `['title1', 'title2', ...]`.
"""

doc_vbound_percentile = """\
vmin
    The value representing the lower limit of the color scale. Values smaller than vmin are plotted
    with the same color as vmin. vmin can be a number, a string, a function or `None`. If
    vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
    For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
    vmin is interpreted as the return value of the function over the list of values to plot.
    For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
    np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
    minimum value is used as defined by matplotlib `scatter` function. When making multiple
    plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
vmax
    The value representing the upper limit of the color scale. The format is the same as for `vmin`.
vcenter
    The value representing the center of the color scale. Useful for diverging colormaps.
    The format is the same as for `vmin`.
    Example: sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')\
"""

doc_vboundnorm = """\
vmin
    The value representing the lower limit of the color scale. Values smaller than vmin are plotted
    with the same color as vmin.
vmax
    The value representing the upper limit of the color scale. Values larger than vmax are plotted
    with the same color as vmax.
vcenter
    The value representing the center of the color scale. Useful for diverging colormaps.
norm
    Custom color normalization object from matplotlib. See
    `https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.\
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
{doc_vbound_percentile}
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
categories_order
    Order in which to show the categories. Note: add_dendrogram or add_totals
    can change the categories order.
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

doc_rank_genes_groups_plot_args = """\
adata
    Annotated data matrix.
groups
    The groups for which to show the gene ranking.
n_genes
    Number of genes to show. This can be a negative number to show for
    example the down regulated genes. eg: num_genes=-10. Is ignored if
    `gene_names` is passed.
gene_symbols
    Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
    refer to the index column of the `.var` DataFrame. Setting this option allows
    alternative names to be used.
groupby
    The key of the observation grouping to consider. By default,
    the groupby is chosen from the rank genes groups parameter but
    other groupby options can be used.  It is expected that
    groupby is a categorical. If groupby is not a categorical observation,
    it would be subdivided into `num_categories` (see :func:`~scanpy.pl.dotplot`).
min_logfoldchange
    Value to filter genes in groups if their logfoldchange is less than the
    min_logfoldchange
key
    Key used to store the ranking results in `adata.uns`.\
"""

doc_rank_genes_groups_values_to_plot = """\
values_to_plot
    Instead of the mean gene value, plot the values computed by `sc.rank_genes_groups`.
    The options are: ['scores', 'logfoldchanges', 'pvals', 'pvals_adj',
    'log10_pvals', 'log10_pvals_adj']. When plotting logfoldchanges a divergent
    colormap is recommended. See examples below.
var_names
    Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes)
    to check their fold changes or p-values, instead of the top/bottom genes. The
    var_names could be a dictionary or a list as in :func:`~scanpy.pl.dotplot` or
    :func:`~scanpy.pl.matrixplot`. See examples below.\
"""

doc_scatter_spatial = """\
library_id
    library_id for Visium data, e.g. key in `adata.uns["spatial"]`.
img_key
    Key for image data, used to get `img` and `scale_factor` from `"images"`
    and `"scalefactors"` entires for this library. To use spatial coordinates,
    but not plot an image, pass `img_key=None`.
img
    image data to plot, overrides `img_key`.
scale_factor
    Scaling factor used to map from coordinate space to pixel space.
    Found by default if `library_id` and `img_key` can be resolved.
    Otherwise defaults to `1.`.
spot_size
    Diameter of spot (in coordinate space) for each point. Diameter
    in pixels of the spots will be `size * spot_size * scale_factor`.
    This argument is required if it cannot be resolved from library info.
crop_coord
    Coordinates to use for cropping the image (left, right, top, bottom).
    These coordinates are expected to be in pixel space (same as `basis`)
    and will be transformed by `scale_factor`.
    If not provided, image is automatically cropped to bounds of `basis`,
    plus a border.
alpha_img
    Alpha value for image.
bw
    Plot image data in gray scale.\
"""

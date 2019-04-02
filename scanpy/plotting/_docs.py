"""Shared docstrings for plotting function parameters.
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
    Use `.raw` attribute of `adata` for coloring with gene expression. If
    `None`, uses `.raw` if present.
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
edges_color : matplotlib color(s), optional (default: 'grey')
    Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.
arrows
    Show arrows (requires to run :func:`~scanpy.api.tl.rna_velocity` before).
arrows_kwds
    Passed to :func:`~matplotlib.axes.Axes.quiver`\
"""


_doc_scatter_common = """\
sort_order
    For continuous annotations used as color parameter, plot data points
    with higher values on top of others.
groups
    Restrict to a few categories in categorical observation annotation.
    The default is not to restrict to any groups.
components
    For instance, `['1,2', '2,3']`. To plot all available components use
    `components='all'`.
projection : {'2d', '3d'}, optional (default: '2d')
    Projection of plot.
legend_loc
    Location of legend, either 'on data', 'right margin' or valid keywords for
    `matplotlib.legend`.
legend_fontsize
    Legend font size.
legend_fontweight : {'normal', 'bold', ...}, optional (default: `None`)
    Legend font weight. Defaults to 'bold' if `legend_loc == 'on data'`,
    otherwise to 'normal'. Available are `['light', 'normal', 'medium',
    'semibold', 'bold', 'heavy', 'black']`.
size
    Point size. If `None`, is automatically computed.
color_map
    Color map to use for continous variables. Anything that works for `cmap`
    argument of `pyplot.scatter` should work here (e.g. `"magma"`, `"viridis"`,
    `mpl.cm.cividis`). If `None` value of `mpl.rcParams["image.cmap"]` is used.
palette
    Colors to use for plotting categorical annotation groups. The palette can be
    a valid `matplotlib.pyplot.colormap` name like `'Set2'` or `'tab20'`, a list
    of colors like `['red', '#ccdd11', (0.1, 0.2, 1)]` or a Cycler object. If
    `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
    variable already has colors stored in `adata.uns["{var}_colors"]`. If
    provided, values of `adata.uns["{var}_colors"]` will be set by this palette.
frameon
    Draw a frame around the scatter plot. Defaults to value set in
    :func:`~scanpy.api.tl.set_figure_params`, defaults to `True`.
"""
_doc_scatter_panels = """\
ncols
    Number of panels per row.
wspace
    Adjust the width of the space between multiple panels.
hspace
    Adjust the height of the space between multiple panels.
"""
_doc_scatter_meta = """\
title
    Provide title for panels either as string or list of strings,
    e.g. `['title1', 'title2', ...]`.
kwargs : further keyword arguments, optional
    Arguments to pass to :func:`matplotlib.pyplot.scatter`,
    for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).
return_fig
    Return the matplotlib figure.\
"""

# temporarily add a special variable doc_scatter_temp for pl.scatter
# because currently pl.scatter does not accept ncols, wspace, and hspace
doc_scatter_temp = _doc_scatter_common + _doc_scatter_meta
doc_scatter_bulk = _doc_scatter_common + _doc_scatter_panels + _doc_scatter_meta


doc_show_save_ax = """\
show
     Show the plot, do not return axis.
save
    If `True` or a `str`, save the figure. A string is appended to the default
    filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
ax
    A matplotlib axes object. Only works if plotting a single component.\
"""


doc_common_plot_args = """\
adata : :class:`~anndata.AnnData`
    Annotated data matrix.
var_names : `str` or list of `str`
    `var_names` should be a valid subset of  `adata.var_names`.
groupby : `str` or `None`, optional (default: `None`)
    The key of the observation grouping to consider.
log : `bool`, optional (default: `False`)
    Plot on logarithmic axis.
use_raw : `bool`, optional (default: `None`)
    Use `raw` attribute of `adata` if present.
num_categories : `int`, optional (default: `7`)
    Only used if groupby observation is not categorical. This value
    determines the number of groups into which the groupby observation
    should be subdivided.
figsize : (`float`, `float`), optional (default: `None`)
    Figure size when multi_panel = True. Otherwise the rcParam['figure.figsize] value is used.
    Format is (width, height)
dendrogram: `bool` or `str`, optional (default, `False`)
    If True or a valid dendrogram key, a dendrogram based on the hierarchical clustering 
    between the `groupby` categories is added. The dendrogram information is computed
    using :ref:`scanpy.tl.dendrogram`. If `tl.dendrogram` has not been called previously
    the function is called with default parameters.
gene_symbols : string, optional (default: `None`)
    Column name in `.var` DataFrame that stores gene symbols. By default `var_names` 
    refer to the index column of the `.var` DataFrame. Setting this option allows
    alternative names to be used.
var_group_positions :  list of `tuples`.
    Use this parameter to highlight groups of `var_names`.
    This will draw a 'bracket' or a color block between the given start and end positions. If the
    parameter `var_group_labels` is set, the corresponding labels are added on
    top/left. E.g. var_group_positions = [(4,10)] will add a bracket
    between the fourth var_name and the tenth var_name. By giving more
    positions, more brackets/color blocks are drawn.
var_group_labels : list of `str`
    Labels for each of the var_group_positions that want to be highlighted.
var_group_rotation : `float` (default: `None`)
    Label rotation degrees. By default, labels larger than 4 characters are rotated 90 degrees
layer: `str`, (default `None`)
    Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
    If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
    then the layer is plotted. `layer` takes precedence over `use_raw`.\
"""

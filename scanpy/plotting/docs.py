"""Shared docstrings for plotting function parameters.
"""


doc_adata_color_etc = """\
adata : :class:`~anndata.AnnData`
    Annotated data matrix.
color : string or list of strings, optional (default: `None`)
    Keys for annotations of observations/cells or variables/genes, e.g.,
    `'ann1'` or `['ann1', 'ann2']`.
use_raw : `bool`, optional (default: `None`)
    Use `.raw` attribute of `adata` for coloring with gene expression. If
    `None`, uses `.raw` if present.\
"""


doc_edges_arrows = """\
edges : `bool`, optional (default: `False`)
    Show edges.
edges_width : `float`, optional (default: 0.1)
    Width of edges.
edges_color : matplotlib color, optional (default: 'grey')
    Color of edges.
arrows : `bool`, optional (default: `False`)
    Show arrows (requires to run :func:`~scanpy.api.tl.rna_velocity` before).\
"""


doc_scatter_bulk = """\
sort_order : `bool`, optional (default: `True`)
    For continuous annotations used as color parameter, plot data points
    with higher values on top of others.
groups : `str`, optional (default: `all groups`)
    Restrict to a few categories in categorical observation annotation.
components : `str` or list of `str`, optional (default: '1,2')
    For instance, `['1,2', '2,3']`. To plot all available components use
    `components='all'`.
projection : {'2d', '3d'}, optional (default: '2d')
    Projection of plot.
legend_loc : `str`, optional (default: 'right margin')
    Location of legend, either 'on data', 'right margin' or valid keywords for
    `matplotlib.legend`.
legend_fontsize : `int`, optional (default: `None`)
    Legend font size.
legend_fontweight : {'normal', 'bold', ...}, optional (default: `None`)
    Legend font weight. Defaults to 'bold' if `legend_loc == 'on data'`,
    otherwise to 'normal'. Available are `['light', 'normal', 'medium',
    'semibold', 'bold', 'heavy', 'black']`.
size : `float` (default: `None`)
    Point size. If `None`, is automatically computed.
palette : `str`, list of `str`, or `Cycler` optional (default: `None`)
    Colors to use for plotting categorical annotation groups. The palette can be
    a valid `matplotlib.pyplot.colormap` name like `'Set2'` or `'tab20'`, a list
    of colors like `['red', '#ccdd11', (0.1, 0.2, 1)]` or a Cycler object.
frameon : `bool`, optional (default: `None`)
    Draw a frame around the scatter plot. Defaults to value set in
    :func:`~scanpy.api.tl.set_figure_params`, defaults to `True`.
ncols : `int` (default: 4)
    Number of panels per row.
wspace : `float` (default: 0.1)
    Adjust the width of the space between multiple panels.
hspace : `float` (default: 0.25)
    Adjust the height of the space between multiple panels.
title : `str` or list of `str`, optional (default: `None`)
    Provide title for panels either as, e.g. `['title1', 'title2', ...]`.
kwargs : further keyword arguments, optional
    Arguments to pass to `matplotlib.pyplot.scatter`, for instance, the color
    map (e.g. `cmap='viridis'`) or the maximum and minimum values
    (e.g. `vmin=-2, vmax=5`).\
"""


doc_show_save_ax = """\
show : `bool`, optional (default: `None`)
     Show the plot, do not return axis.
save : `bool` or `str`, optional (default: `None`)
    If `True` or a `str`, save the figure. A string is appended to the default
    filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
ax : `matplotlib.Axes`, optional (default: `None`)
    A matplotlib axes object. Only works if plotting a single component.\
"""

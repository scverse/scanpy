"""Settings
"""

verbosity = 1
"""Set global verbosity level.

Level 0: only show 'error' messages.
Level 1: also show 'warning' messages.
Level 2: also show 'info' messages.
Level 3: also show 'hint' messages.
Level 4: also show very detailed progress.
Level 5: also show even more detailed progress.
etc.
"""

plot_suffix = ''
"""Global suffix that is appended to figure filenames.
"""

file_format_data = 'h5ad'
"""File format for saving AnnData objects.

Allowed are 'txt', 'csv' (comma separated value file) for exporting and 'h5ad'
(hdf5) for lossless saving.
"""

file_format_figs = 'pdf'
"""File format for saving figures.

For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
`matplotlib.pyplot.savefig`).
"""

autosave = False
"""Save plots/figures as files in directory 'figs'.

Do not show plots/figures interactively.
"""

autoshow = True
"""Show all plots/figures automatically if autosave == False.

There is no need to call the matplotlib pl.show() in this case.
"""

writedir = './write/'
"""Directory where the function scanpy.write writes to by default.
"""

cachedir = './cache/'
"""Default cache directory.
"""

figdir = './figures/'
"""Directory where plots are saved.
"""

max_memory = 15
"""Maximal memory usage in Gigabyte.

Is currently not well respected....
"""

n_jobs = 1
"""Default number of jobs/ CPUs to use for parallel computing.
"""

logfile = ''
"""Name of logfile. By default is set to '' and writes to standard output."""

categories_to_ignore = ['N/A', 'dontknow', 'no_gate', '?']
"""Categories that are omitted in plotting etc.
"""

_frameon = True
"""See set_figure_params.
"""


# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------


def set_figure_params(scanpy=True, dpi=80, dpi_save=150, frameon=True, vector_friendly=True,
                      fontsize=14, color_map=None, format='pdf',
                      transparent=False, ipython_format='png2x'):
    """Set resolution/size, styling and format of figures.

    Parameters
    ----------
    scanpy : `bool`, optional (default: `True`)
        Init default values for ``matplotlib.rcParams`` suited for Scanpy.
    dpi : `int`, optional (default: `80`)
        Resolution of rendered figures - this influences the size of figures in notebooks.
    dpi_save : `int`, optional (default: `150`)
        Resolution of saved figures. This should typically be higher to achieve
        publication quality.
    frameon : `bool`, optional (default: `True`)
        Add frames and axes labels to scatter plots.
    vector_friendly : `bool`, optional (default: `True`)
        Plot scatter plots using `png` backend even when exporting as `pdf` or `svg`.
    fontsize : `int`, optional (default: 14)
        Set the fontsize for several `rcParams` entries. Ignored if `scanpy=False`.
    color_map : `str`, optional (default: `None`)
        Convenience method for setting the default color map. Ignored if `scanpy=False`.
    format : {'png', 'pdf', 'svg', etc.}, optional (default: 'pdf')
        This sets the default format for saving figures: `file_format_figs`.
    transparent : `bool`, optional (default: `True`)
        Save figures with transparent back ground. Sets
        `rcParams['savefig.transparent']`.
    ipython_format : list of `str`, optional (default: 'png2x')
        Only concerns the notebook/IPython environment; see
        `IPython.core.display.set_matplotlib_formats` for more details.
    """
    try:
        import IPython
        IPython.core.display.set_matplotlib_formats(ipython_format)
    except:
        pass
    from matplotlib import rcParams
    global _vector_friendly
    _vector_friendly = vector_friendly
    global file_format_figs
    file_format_figs = format
    if dpi is not None:
        rcParams['figure.dpi'] = dpi
    if dpi_save is not None:
        rcParams['savefig.dpi'] = dpi_save
    if transparent is not None:
        rcParams['savefig.transparent'] = transparent
    if scanpy:
        from .plotting._rcmod import set_rcParams_scanpy
        set_rcParams_scanpy(fontsize=fontsize, color_map=color_map)
    global _frameon
    _frameon = frameon



# ------------------------------------------------------------------------------
# Private global variables & functions
# ------------------------------------------------------------------------------

_vector_friendly = False
"""Set to true if you want to include pngs in svgs and pdfs.
"""

_low_resolution_warning = True
"""Print warning when saving a figure with low resolution."""

def _set_start_time():
    from time import time
    return time()

_start = _set_start_time()
"""Time when the settings module is first imported."""

_previous_time = _start
"""Variable for timing program parts."""

_previous_memory_usage = -1
"""Stores the previous memory usage."""


def _is_run_from_ipython():
    """Determines whether run from Ipython.

    Only affects progress bars.
    """
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

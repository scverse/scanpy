"""Settings, Logging and Timing

Sets global variables like verbosity, manages logging and timing.
"""
# Note
# ----
# The very first version (tracking cpu time) of this was based on
# http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution

import time
from functools import reduce

# directory for configuration files etc.
# if not os.path.exists('.scanpy/'):
#     os.makedirs('.scanpy/')

# --------------------------------------------------------------------------------
# Global Settings Attributes
# --------------------------------------------------------------------------------

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

run_name = ''
"""Run name. Often associated with a certain way of preprocessing of parameter combination.

All files generated during the run have this name as prefix.
"""

_run_basename = ''
"""Basename of the run.

This usually associated with a certain way of preprocessing the data and running
several tools after that.

Determines the naming of all output files.
"""

_run_suffix = ''
"""Global suffix that is appended to project identifier.
"""

plot_suffix = ''
"""Global suffix that is appended to figure filenames.

Is needed when the computation parameters remain unchanged, but only plotting
parameters are changed.
"""

file_format_data = 'h5ad'
"""File format for saving AnnData objects.

Allowed are 'txt', 'csv' (comma separated value file) for exporting and 'h5'
(hdf5) and 'npz' for importing and exporting.
"""

file_format_figs = 'png'
"""File format for saving figures.

For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
matplotlib.pyplot.savefig).
"""

autosave = False  # This should be renamed "autosave"
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
"""

n_jobs = 2
"""Maximal number of jobs/ CPUs to use for parallel computing.
"""

logfile = ''
"""Name of logfile. By default is set to '' and writes to standard output."""

import __main__ as main
is_run_from_file = not hasattr(main, '__file__')
"""Determines whether run from file.

Currently only affects whether total computation
time since importing this module is output after leaving the session.
"""

def _is_run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
"""Determines whether run from Ipython.

Only affects progress bars
"""
is_run_from_ipython = _is_run_from_ipython()


_dpi = 400
"""Resolution of png figures.

We need this global variable as, for example, Seaborn resets rcParams['savefig.dpi'].
"""

categories_to_ignore = ['N/A', 'dontknow', 'no_gate', '?']
"""Categories that are omitted in plotting etc.
"""

_low_resolution_warning = True
"""Print warning when saving a figure with low resolution."""


# --------------------------------------------------------------------------------
# Global Setting Functions
# --------------------------------------------------------------------------------


def set_figure_params(dpi=None, scanpy=True, figure_formats=['png2x']):
    """Set resolution/size, styling and format of figures.

    Parameters
    ----------
    dpi : `int`, optional
        Resolution of png output in dots per inch.
    scanpy : `bool`, optional (default: `True`)
        Run :func:`scanpy.api.pl.set_rcParams_scanpy` to init Scanpy rcParams.
    figure_formats : list of strings
        Only concerns the IPython environment; see
        `IPython.core.display.set_matplotlib_formats` for more details. For
        setting the default format for saving figures, directly set
        `file_format_figs`.
    """
    try:
        import IPython
        IPython.core.display.set_matplotlib_formats(*figure_formats)
    except:
        pass
    from matplotlib import rcParams
    global _dpi
    if dpi is not None:
        _dpi = dpi
        rcParams['savefig.dpi'] = _dpi
        rcParams['figure.dpi'] = _dpi
    if scanpy:
        from .plotting.rcmod import set_rcParams_scanpy
        set_rcParams_scanpy()


# ------------------------------------------------------------------------------
# Private global variables
# ------------------------------------------------------------------------------

_start = time.time()
"""Time when the settings module is first imported."""

_previous_time = _start
"""Variable for timing program parts."""

_previous_memory_usage = -1
"""Stores the previous memory usage."""


# --------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------


def mi(*msg, end='\n'):
    """Write message to log output, ignoring the verbosity level.

    Parameters
    ----------
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if logfile == '':
        # in python 3, the following works
        print(*msg, end=end)
        # we do not bother for compat anymore?
        # due to compatibility with the print statement in python 2 we choose
        # print(' '.join([str(m) for m in msg]))
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with open(logfile, 'a') as f:
            f.write(out + end)


def m(v=0, *msg):
    """Write message to log output, depending on verbosity level.

    Now is deprecatd. See logging.msg().

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if v < verbosity - 2:
        mi(*msg)


def mt(v=0, *msg, start=False):
    """Write message to log output and show computation time.

    Depends on chosen verbosity level.

    Now is deprecatd. See logging.m().

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg : str
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if v < verbosity - 2:
        global _previous_time
        now = time.time()
        if start:
            _previous_time = now
        elapsed = now - _previous_time
        _previous_time = now
        mi(_sec_to_str(elapsed), '-', *msg)


def _sec_to_str(t):
    """Format time in seconds.

    Parameters
    ----------
    t : int
        Time in seconds.
    """
    return "%d:%02d:%02d.%02d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*100,), 100, 60, 60])

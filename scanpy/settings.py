# Author: F. Alex Wolf (http://falexwolf.de)
"""Settings, Logging and Timing

Sets global variables like verbosity, manages logging and timing.
"""
# Note
# ----
# The very first version (tracking cpu time) of this was based on
# http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution

import atexit
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

file_format_data = 'h5'
"""File format for saving AnnData objects.

Allowed are 'txt', 'csv' (comma separated value file) for exporting and 'h5'
(hdf5) and 'npz' for importing and exporting.
"""

file_format_figs = 'png'
"""File format for saving figures.

For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
matplotlib.pyplot.savefig).
"""

recompute = 'read'
"""If set to 'none', use the results of previous computations.

Recompute and overwrite previous result and preprocessing files.
"""

savefigs = False
"""Save plots/figures as files in directory 'figs'.

Do not show plots/figures interactively.
"""

autoshow = True
"""Show all plots/figures automatically if savefigs == False.

There is no need to call the matplotlib pl.show() in this case.
"""

writedir = './write/'
"""Directory where the function scanpy.write writes to by default.
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


def set_figure_params(dpi=None, figure_formats=['png2x']):
    """Set resolution and format of figures.

    Parameters
    ----------
    dpi : int, optional
        Resolution of png output in dots per inch.
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
    if dpi is not None: _dpi = dpi
    # need to set the following two lines as older Jupyter notebooks seem to use
    # 'savefig.dpi' and more rescent ones 'figure.dpi'
    rcParams['savefig.dpi'] = _dpi
    rcParams['figure.dpi'] = _dpi


set_dpi = set_figure_params
"""Deprecated: merely for backwards compatibility. See `set_figure_params` instead."""

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
# Command-line arguments for global attributes
# --------------------------------------------------------------------------------


def add_args(p):
    """
    Add arguments that affect the global variables in settings.

    Parameters
    ----------
    p : argparse.ArgumentParser
        Input parser.

    Returns
    -------
    p : argparse.ArgumentParser
        Updated parser.
    """
    aa = p.add_argument_group('Save figures').add_argument
    aa('-s', '--savefigs',
       type=str, default='', const=file_format_figs, nargs='?', metavar='ext',
       help='Save figures either as "png", "svg" or "pdf". Just providing '
            '"--savefigs" will save to "png" (default: do not save figures).')
    aa('--figdir',
       type=str, default=figdir, metavar='dir',
       help='Change figure directory and save figures (default: %(default)s).')
    aa = p.add_argument_group('Run a tool repeatedly to try out different parameters').add_argument
    aa('-r', '--recomp',
       type=str, default='none', const='tool', nargs='?', metavar='x',
       help='Recompute and overwrite result files of previous calculations. '
            'Just providing "-r" only recomputes the tool, "-r pp" also '
            'recomputes preprocessing, "-r read" also rereads data files '
            '(default: do not recompute).')
    aa('--suffix',
       type=str, default='', metavar='suffix',
       help='Is appended to ppkey in result filename (default: "").')
    aa('--psuffix',
       type=str, default='', metavar='psuffix',
       help='Is appended to suffix. Useful when only changing plotting '
            'parameters (default: "").')
    aa = p.add_argument_group('Do subsampling to speedup computations').add_argument
    aa('-ss', '--subsample',
       type=int, default=1, metavar='i',
       help='Pass integer i > 1 if you want to use a fraction of 1/i '
            'of the data (default: %(default)d).')
    aa = p.add_argument_group('General settings').add_argument
    aa('-v', '--verbosity',
       type=int, default=verbosity, metavar='v',
       help='Pass v = 2 (no hints "-->", only info, warnings, errors) for less output, '
            'v = 1 (no info), v = 0 (no warnings, only errors) '
            'and v > 3 for more output (default: %(default)d).')
    aa('--max_memory',
       type=int, default=max_memory, metavar='m',
       help='Pass maximal memory usage in GB (default: %(default)s).')
    aa('--n_jobs',
       type=int, default=n_jobs, metavar='n',
       help='Maximal number of CPUs to use (default: %(default)s).')
    aa('-ff', '--fileformat',
       type=str, default=file_format_data, metavar='ext',
       help='Pass file format for exporting results, either "csv", '
            '"txt", "h5" or "npz" (default: %(default)s).')
    aa('--writedir',
       type=str, default=writedir, metavar='dir',
       help='Change write directory (default: %(default)s).')
    aa('-l', '--logfile',
       action='store_const', default=False, const=True,
       help='Write to logfile instead of standard output.')
    aa = p.add_argument_group('Other').add_argument
    aa('-h', '--help',
       action='help',
       help='Show this help message and exit.')
    return p


def process_args(args):
    """
    Init global variables from dictionary of arguments.

    Parameters
    ----------
    args : dict
        Dictionary of command-line arguments defined in add_args.
    """

    # set the arguments as global variables
    global _run_suffix
    _run_suffix = args['suffix']
    args.pop('suffix')

    global plot_suffix
    plot_suffix = args['psuffix']
    args.pop('psuffix')

    global recompute
    recompute = args['recomp']
    args.pop('recomp')

    global verbosity
    verbosity = args['verbosity']
    args.pop('verbosity')

    global savefigs
    global file_format_figs
    if args['savefigs'] == '':
        savefigs = False
    else:
        savefigs = True
        file_format_figs = args['savefigs']
        if args['savefigs'] == 'png':
            import matplotlib
            matplotlib.use('Agg')
            print('... you passed `--savefigs png`, now '
                  'using matplotlib "Agg" backend')
    args.pop('savefigs')

    global figdir
    if figdir != args['figdir']:
        savefigs = True
    figdir = args['figdir']
    if figdir[-1] != '/':
        figdir += '/'
    from os import path, makedirs
    if not path.exists(figdir):
        print('creating directory', figdir, 'for saving figures')
        makedirs(figdir)
    args.pop('figdir')

    global file_format_data
    file_format_data = args['fileformat']
    args.pop('fileformat')

    global writedir
    writedir = args['writedir']
    if writedir[-1] != '/':
        writedir += '/'
    args.pop('writedir')

    global max_memory
    max_memory = args['max_memory']
    args.pop('max_memory')

    global n_jobs
    n_jobs = args['n_jobs']
    args.pop('n_jobs')

    global _run_basename, run_name
    _run_basename = args['run_name']
    run_name = _run_basename + _run_suffix

    global logfile
    if args['logfile']:
        logfile = writedir + run_name + '_log.txt'

    return args


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

    Now is deprecatd. See logging.m().

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
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*1000,), 1000, 60, 60])


def _terminate():
    """Function called when program terminates.

    Similar to mt, but writes total runtime.
    """
    if verbosity > 0:
        now = time.time()
        elapsed_since_start = now - _start
        mi(29*"_")
        mi(_sec_to_str(elapsed_since_start), '- total wall time')


# report total runtime upon shutdown
if not is_run_from_file:
    atexit.register(_terminate)

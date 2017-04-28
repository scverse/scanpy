# Author: F. Alex Wolf (http://falexwolf.de)
"""
Settings, Logging and Timing

Sets global variables like verbosity, manages logging and timing.
"""
# Note
# ----
# The very first version (tracking cpu time) of this was based on
# http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution

import os
if not os.path.exists('.scanpy/'):  # directory for configuration files etc.
    os.makedirs('.scanpy/')
import matplotlib
if 'DISPLAY' not in os.environ:  # login via ssh but no xserver
    matplotlib.use('Agg')
    print('... WARNING: did not find DISPLAY variable needed for interactive plotting\n'
          '--> try ssh with `-X` or `-Y`')
import atexit
import time
from functools import reduce
from matplotlib import rcParams

#--------------------------------------------------------------------------------
# Global Settings Attributes
#--------------------------------------------------------------------------------

verbosity = 1
"""Set global verbosity level, choose from {0,...,6}.
"""

exkey = ''
"""Global example key.
"""

suffix = ''
"""Global suffix that is appended to basekey of output and figure files.
"""

plotsuffix = ''
"""Global suffix that is appended to figure filenames.

Is needed when the computation parameters remain unchanged, but only plotting
parameters are changed.
"""

file_format_data = 'h5'
"""File format for saving AnnData objects.

Allowed are 'h5' (hdf5), 'xlsx' (Excel) or 'csv' (comma separated value
file).
"""

file_format_figures = 'png'
"""File format for saving figures.

For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
matplotlib.pyplot.savefig).
"""

recompute = 'none'
"""Don't use the results of previous calculations.

Recompute and overwrite previous result and preprocessing files.
"""

savefigs = False
"""Save plots/figures as files in directory 'figs'.

Do not show plots/figures interactively.
"""

autoshow = True
"""Show all plots/figures automatically if savefigs == False.

There is no need to call sc.show() in this case.
"""

writedir = './write/'
"""Directory where the function scanpy.write writes to by default.
"""

figdir = './figs/'
"""Directory where plots are saved.
"""

basekey = ''
"""Basename for file reading and writing.
"""

max_memory = 15
"""Maximal memory usage in Gigabyte.
"""


#--------------------------------------------------------------------------------
# Global Setting Functions
#--------------------------------------------------------------------------------


def set_logfile(filename=''):
    """
    Define filename of logfile.

    If not defined, log output will be to the standard output.

    Parameters
    ----------
    filename : str
        Filename of
    """
    global _logfilename, verbosity
    _logfilename = filename
    # if providing a logfile name, automatically set verbosity to a very high level
    verbosity = 5


def set_dpi(dpi=200):
    """
    Set resolution of png figures.

    Parameters
    ----------
    dpi : int, optional
        Resolution of png output in dots per inch.
    """
    # default setting as in scanpy.plot
    rcParams['savefig.dpi'] = dpi


def set_jupyter():
    """
    Update figure resolution for use in jupyter notebook.

    Avoids that figures get displayed too large. To set a specific value for the
    resolution, use the dpi function.
    """
    dpi(60)
    global autoshow
    autoshow = True


# ------------------------------------------------------------------------------
# Private global variables
# ------------------------------------------------------------------------------

import __main__ as main
_start = time.time()
"""Time when the settings module is first imported."""
_intermediate = _start
"""Variable for timing program parts."""
_logfilename = ''
"""Name of logfile."""
_is_interactive = not hasattr(main, '__file__')
"""Determines whether run as file or imported as package."""


#--------------------------------------------------------------------------------
# Command-line arguments for global attributes
#--------------------------------------------------------------------------------


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
       type=str, default='', const=file_format_figures, nargs='?', metavar='ext',
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
       help='Is appended to exkey in result filename (default: "").')
    aa('--psuffix',
       type=str, default='', metavar='psuffix',
       help='Is appended to suffix. Useful when only changing plotting '
            'parameters (default: "").')
    aa = p.add_argument_group('Do subsampling to speedup computations').add_argument
    aa('-ss', '--subsample',
       type=int, default=1, metavar='i',
       help='Specify integer i > 1 if you want to use a fraction of 1/i'
            ' of the data (default: %(default)d).')
    aa = p.add_argument_group('General settings (all saved in .scanpy/config)').add_argument
    aa('-v', '--verbosity',
       type=int, default=1, metavar='v',
       help='Specify v = 0 for no output and v > 1 for more output'
            ' (default: %(default)d).')
    aa('--max_memory',
       type=int, default=16, metavar='m',
       help='Specify maximal memory usage in GB (default: %(default)s).')
    aa('--fileformat',
       type=str, default=file_format_data, metavar='ext',
       help='Specify file format for saving results, either "h5", "csv", '
            '"txt" or "npz" (default: %(default)s).')
    aa('--writedir',
       type=str, default=writedir, metavar='dir',
       help='Change write directory (default: %(default)s).')
    aa('--logfile',
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
    global suffix
    suffix = args['suffix']
    args.pop('suffix')

    global plotsuffix
    plotsuffix = args['psuffix']
    args.pop('psuffix')

    global recompute
    recompute = args['recomp']
    args.pop('recomp')

    global verbosity
    verbosity = args['verbosity']
    args.pop('verbosity')

    global savefigs
    global file_format_figures
    if args['savefigs'] == '':
        savefigs = False
    else:
        savefigs = True
        file_format_figures = args['savefigs']
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

    return args


#--------------------------------------------------------------------------------
# Logging
#--------------------------------------------------------------------------------


def m(v=0,*msg):
    """
    Write message to log output, depending on verbosity level.

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if verbosity > v:
        mi(*msg)


def mi(*msg):
    """
    Write message to log output, ignoring the verbosity level.

    Parameters
    ----------
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if _logfilename == '':
        # in python 3, the following works
        # print(*msg)
        # due to compatibility with the print statement in python 2 we choose
        print(' '.join([str(m) for m in msg]))
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with open(_logfilename) as f:
            f.write(out + '\n')


def mt(v=0,*msg, start=False):
    """
    Write message to log output and show computation time.

    Depends on chosen verbosity level.

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg : str
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if verbosity > v:
        global _intermediate
        now = time.time()
        if start:
            _intermediate = now
        elapsed_since_start = now - _start
        elapsed = now - _intermediate
        _intermediate = now
        mi(_sec_to_str(elapsed),'-',*msg)


def _sec_to_str(t):
    """
    Format time in seconds.

    Parameters
    ----------
    t : int
        Time in seconds.
    """
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

def _terminate():
    """
    Function called when program terminates.

    Similar to mt, but writes total runtime.
    """
    if verbosity > 0:
        now = time.time()
        elapsed_since_start = now - _start
        mi(27*"_")
        mi(_sec_to_str(elapsed_since_start),'- total wall time')


# report total runtime upon shutdown
if not _is_interactive:
    atexit.register(_terminate)


def _jupyter_deprecated(do=True):
    """
    Update figure params for particular environments like jupyter.
    """

    fscale = 1
    fontsize = 14
    rcParams['savefig.dpi'] = 100

    if do:
        fscale = 0.375
        fontsize = 6

    # figure unit length and figure scale
    ful = fscale*4
    fontsize = fscale*14

    rcParams['lines.linewidth'] = fscale*1.5
    rcParams['lines.markersize'] = fscale**2*6
    rcParams['lines.markeredgewidth'] = fscale**2*1

    rcParams['figure.figsize'] = (1.25*ful,ful)
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = 0.92*fontsize
    rcParams['axes.titlesize'] = fontsize

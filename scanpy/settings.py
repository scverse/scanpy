# coding: utf-8
"""
Settings and Logfile 
====================

From package Scanpy (https://github.com/theislab/scanpy).
Written in Python 3 (compatible with 2).
Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).

Sets global variables like verbosity, manages log output and timing.

Note
----
The very first version (tracking cpu time) of this was based on
http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution
"""   

import atexit
from time import clock
from functools import reduce
from matplotlib import rcParams

#--------------------------------------------------------------------------------
# Global Settings
#--------------------------------------------------------------------------------

verbosity = 1
""" Set global verbosity level, choose from {0,...,6}. """

exkey = ''
""" Global example key.
"""

suffix = ''
""" Global suffix, which is appended to basekey of output and figure files.
"""

plotsuffix = ''
""" Global suffix which is appended to figure filenames.

Is needed when the computation parameters remain unchanged, but only plotting
parameters are changed.
"""

extd = 'h5'
""" Global file extension format for data storage. 

Allowed are 'h5' (hdf5), 'xlsx' (Excel) or 'csv' (comma separated value
file).
"""

extf = 'png'
""" Global file extension for saving figures.

Recommended are 'png' and 'pdf'. Many other formats work as well (see
matplotlib.pyplot.savefig).
"""

recompute = 'none'
""" Don't use the results of previous calculations.

Recompute and overwrite previous result and preprocessing files.  
"""

savefigs = False
""" Save plots/figures as files in directory 'figs'.

Do not show plots/figures interactively.
"""

autoshow = True
""" Show all plots/figures automatically if savefigs == False. 

There is no need to call sc.show() in this case.
"""

writedir = 'write/'
""" Directory where the function scanpy.write writes to by default.
"""

figdir = 'figs/'
""" Directory where plots are saved.
"""

fsig = ''
""" File signature.
"""

basekey = ''
""" Basename for file reading and writing.
"""

#--------------------------------------------------------------------------------
# Command-line arguments for global variables in settings
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
    aa = p.add_argument_group('Subsampling speeds up computations').add_argument
    aa('-s', '--subsample',
       type=int, default=1, metavar='s',
       help='Specify integer s > 1 if you want to use a fraction of 1/s'
            ' of the data (default: %(default)d).')
    aa = p.add_argument_group('Save figures').add_argument
    aa('--savefigs',
       type=str, default='', const=extf, nargs='?', metavar='extf',
       help='Save figures to files. With the exception of interactive sessions, '
            'and do not show interactive plots anymore. '
            'Specify the file format via the extension, e.g. "pdf", '
            '"svg", "png". Just providing "--savefigs" will save to "png" '
            '(default: do not save figures).')
    aa('--figdir',
       type=str, default=figdir, metavar='dir',
       help='Change figure directory (default: %(default)s).')
    aa = p.add_argument_group('Run a tool repeatedly, to try out different parameters').add_argument
    aa('-r', '--recompute',
       type=str, default='none', const='tool', nargs='?', metavar='r',
       help='Recompute and overwrite result files of previous calculations. '
            'Just providing "-r" recomputes the tool, "-r all" also recomputes preprocessing.'
            '(default: do not recompute).'
       )
    aa('--suffix',
       type=str, default='', metavar='suffix',
       help='Specify suffix to append to example key'
            ' (default: "").')
    aa('--plotsuffix',
       type=str, default='', metavar='psuffix',
       help='Specify suffix to append to suffix, when you simply change plotting parameters'
            ' (default: "").')
    aa = p.add_argument_group('General settings').add_argument
    aa('-h', '--help',
       action='help',
       help='Show this help message and exit.')
    aa('-v', '--verbosity',
       type=int, default=1, metavar='v',
       help='Specify v = 0 for no output and v > 1 for more output.'
            ' (default: %(default)d).')
    aa('--logfile',
       action='store_const', default=False, const=True,
       help='Write to logfile instead of to standard output.')
    aa('--fileform',
       type=str, default=extd, metavar='extd',
       help='Specify file extension and by that'
            ' file format for saving results, either "h5" or "xlsx".'
            ' (default: %(default)s).')
    aa('--writedir',
       type=str, default=writedir, metavar='dir',
       help='Change outfile directory (default: %(default)s).')

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
    plotsuffix = args['plotsuffix']
    args.pop('plotsuffix')

    global recompute
    recompute = args['recompute']
    args.pop('recompute')

    global verbosity
    verbosity = args['verbosity']
    args.pop('verbosity')

    global savefigs
    global extf
    if args['savefigs'] == '':
        savefigs = False
    else:   
        savefigs = True
        extf = args['savefigs']
    args.pop('savefigs')

    global figdir
    figdir = args['figdir'] + '/'
    if figdir[-1] != '/':
        figdir += '/'
    from os import path, makedirs
    if not path.exists(figdir):
        makedirs(figdir)    
    args.pop('figdir')

    global extd
    extd = args['fileform']
    args.pop('fileform')

    global writedir
    writedir = args['writedir'] 
    if writedir[-1] != '/':
        writedir += '/'
    args.pop('writedir')

    # from these arguments, init further global variables
    global exkey
    global basekey
    if 'exkey' in args:
        exkey = args['exkey']
        basekey = args['exkey'] + suffix
    else:
        basekey = 'test' + suffix

    # file signature to be appended to each filename
    global fsig
    if args['subsample'] != 1:
        fsig += '_s{:02}'.format(args['subsample'])

    return args

#--------------------------------------------------------------------------------
# Output
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
    if logfilename == '':
        # in python 3, the following works
        # print(*msg)
        # due to compatibility with the print statement in python 2 we choose
        print(' '.join([str(m) for m in msg]))
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with open(logfilename) as f:
            f.write(out + '\n')

def mt(v=0,*msg):
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
        global intermediate
        now = clock()
        elapsed_since_start = now - start
        elapsed = now - intermediate
        intermediate = now
        mi(_sec_to_str(elapsed),'-',*msg)

def logfile(filename=''):
    """ 
    Define filename of logfile. 

    If not defined, log output will be to the standard output.

    Parameters
    ----------
    filename : str
        Filename of 
    """
    global logfilename, verbosity
    logfilename = filename
    # if providing a logfile name, automatically set verbosity to a very high level
    verbosity = 5


def dpi(dpi=200):
    """
    Set resolution of png figures.

    Parameters
    ----------
    dpi : int, optional
        Resolution of png output in dots per inch.
    """
    # default setting as in scanpy.plot
    rcParams['savefig.dpi'] = dpi

def jupyter():
    """
    Update figure resolution for use in jupyter notebook. 

    Avoids that figures get displayed too large. To set a specific value for the
    resolution, use the dpi function.
    """
    dpi(60)
    global autoshow
    autoshow = True

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
        now = clock()
        elapsed_since_start = now - start
        mi(27*"_")
        mi(_sec_to_str(elapsed_since_start),'- total runtime')

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

# further global variables
start = clock() # clock() is deprecated since version python version 3.3
intermediate = start
logfilename = ''
separator = 80*"-"

atexit.register(_terminate)


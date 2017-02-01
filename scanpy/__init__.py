# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Scanpy - Single-Cell Analysis in Python
=======================================

Reference
---------
Wolf, Angerer & Theis, bioRxiv doi:... (2017)
"""

# standard modules
import os
# scientific modules
from .compat.matplotlib import pyplot as pl
# scanpy modules
from . import settings as sett
from . import tools
from . import preprocess
from . import utils
from .tools import get_tool
# reexports
from .utils import read, write, read_params, transpose_ddata
from .exs import exdata, examples, example, annotate
from .preprocess import preprocess, pp
from .preprocess.advanced import subsample
from .tools.diffmap import diffmap
from .tools.tsne import tsne
from .tools.dpt import dpt
from .tools.difftest import difftest
from .tools.sim import sim

__all__ = [
    # example use cases
    'example', # call example
    'exdata', # show available example data
    'exs', # show available example use cases
    # help
    'help', # show help for a given tool
    # elementary operations
    'read',
    'write',
    'annotate',
    # preprocessing
    'preprocess', 'pp',
    'subsample',
    # visualization
    'diffmap',
    'tsne',
    # subgroup identification
    'dpt',
    'ctpaths',
    # differential expression testing
    'difftest',
    # simulation
    'sim'
    # plotting
    'plot',
    'show', # show plots
    # management
    'run'
    # utils
    'transpose_ddata'
]

def plot(*dtools,**kwargs):
    """
    Plot the result of a computation with a tool.

    Parameters
    ----------
    *dtools : dicts
       An arbitrary number of tool dictionaries.
    """
    if sett.savefigs:
        if 'writekey' not in dtools[0]:
            raise ValueError('Need key "writekey" in dict d' 
                             + dtools[0]['type'],' - call sc.write first')

    toolkey = dtools[0]['type']

    # TODO: this is not a good solution
    if toolkey == 'sim':
        dtools = [dtools[0]]

    get_tool(toolkey).plot(*dtools, **kwargs)

def help(toolkey,string=False):
    """
    Display help for tool.
    """
    doc = get_tool(toolkey, func=True).__doc__.replace('\n    ','\n')
    if string:
        return doc
    print(doc)

def show():
    """
    Show plots.
    """
    pl.show()

def run(toolkey, exkey, **kwargs):
    """
    Run example with specified tool, do preprocessing and read/write outfiles.

    Output files store the dictionary returned by the tool. File type is
    determined by variable sett.extd allowed are 'h5' (hdf5), 'xlsx' (Excel) or
    'csv' (comma separated value file).

    If called twice with the same settings the dictionary that is read from the 
    existing output file is returned.

    Parameters
    ----------
    toolkey : str
        Name of the tool.
    exkey : str
        Identifies the example.
    kwargs : keyword arguments
    """
    kwargs['exkey'] = exkey
    run_args(toolkey, kwargs)

def run_args(toolkey, args):
    """
    Run specified tool, do preprocessing and read/write outfiles.

    Output files store the dictionary returned by the tool. File type is
    determined by variable sett.extd allowed are 'h5' (hdf5), 'xlsx' (Excel) or
    'csv' (comma separated value file).

    If called twice with the same settings the existing output file is returned.

    Parameters
    ----------
    toolkey : str
        Name of the tool.
    args : dict containing
        exkey : str
            String that identifies the example use key.

    Returns
    -------
    dfunc : dict of type toolkey
    dadd : dict 
         Additional dict used for plotting in a later step.
    """
    if args['plotparams']:
        if args['plotparams'][0] == 'help':
            from sys import exit
            exit(get_tool(toolkey).plot.__doc__)

    writekey = sett.basekey + '_' + toolkey + sett.fsig
    resultfile = sett.writedir + writekey + '.' + sett.extd
    paramsfile = sett.writedir + writekey + '_params.txt'
    if args['logfile']:
        logfile = sett.writedir + writekey + '_log.txt'
        sett.logfile(logfile)

    if toolkey == 'sim':
        if args['paramsfile'] != '':
            params = read_params(args['paramsfile'])
        else:
            paramsfile_sim = 'sim/' + args['exkey'] + '_params.txt'
            params = read_params(paramsfile_sim)
            sett.m(0,'--> you can specify your custom params file using the option\n'
                     '    "--paramsfile" or provide parameters directly via "--params"')
        if 'writedir' not in params:
            params['writedir'] = sett.writedir + sett.basekey + '_' + toolkey
    else:
        ddata, exmodule = example(args['exkey'], return_module=True)
        params = {}
        # try to load tool parameters from dexamples
        try:
            did_not_find_params_in_exmodule = False
            dexample = exmodule.dexamples[args['exkey']]
            params = {}
            for key in dexample.keys():
                if toolkey in key:
                    params = dexample[key]
        except:
            did_not_find_params_in_exmodule = True
            pass
        # if parameters have been specified in a parameter file
        # update the current param dict with these
        if args['paramsfile'] != '':
            add_params = read_params(args['paramsfile'])
            params = utils.update_params(params, add_params)
        # same if parameters have been specified on the command line
        if args['params']:
            add_params = utils.get_params_from_list(args['params'])
            params = utils.update_params(params, add_params)
        elif did_not_find_params_in_exmodule and args['paramsfile'] != '':
            sett.m(0, 'using default parameters, change them using "--params"')

    # subsampling
    if args['subsample'] != 1:
        ddata = subsample(ddata, args['subsample'])

    # previous tool
    if args['prev'] != '':
        prevkey = sett.basekey + '_' + args['prev'] + sett.fsig
        dprev = read(prevkey)
    # all tools that require a previous tool
    elif toolkey in ['difftest', 'scdg']:
        print('Error: need to provide a tool to option --prev')
        print('--> presumably one for identification of subgroups')
        from sys import exit
        exit(0)
    
    # simply load resultfile
    if os.path.exists(resultfile) and not sett.recompute:
        dtool = read(writekey)
    # call the tool resultfile
    else:
        # TODO: solve this in a nicer way, also get an ordered dict for params
        from inspect import getcallargs
        tool = get_tool(toolkey, func=True)
        if toolkey == 'sim':
            dtool = tool(**params)
            params = getcallargs(tool, **params)
        elif args['prev'] != '':
            dtool = tool(dprev, ddata, **params)
            params = getcallargs(tool, dprev, ddata, **params)
            # TODO: Would be good to name the first argument dprev_or_ddata
            #       in difftest, but this doesn't work
            del params['dprev']
        else:
            dtool = tool(ddata, **params)
            params = getcallargs(tool, ddata, **params)
        if 'ddata' in params:
            del params['ddata']
        elif 'ddata_or_X' in params:
            del params['ddata_or_X']
        dtool['writekey'] = writekey
        write(writekey, dtool)
        sett.m(0, 'wrote result to', utils.get_filename_from_key(writekey))
        # save a copy of the parameters to a file
        utils.write_params(paramsfile, params)

    # plotting and postprocessing
    plotparams = {}
    if args['plotparams']:
        plotparams = utils.get_params_from_list(args['plotparams'])
    if toolkey == 'sim':
        plot(dtool, plotparams)
    else:
        # post-processing specific to example and tool
        postprocess = args['exkey'] + '_' + toolkey
        if postprocess in dir(exmodule) and args['subsample'] == 1:
            dtool = getattr(exmodule, postprocess)(dtool)
            write(writekey, dtool)
        if args['plotkey'] != '':
            plotwritekey = sett.exkey + '_' +  args['plotkey'] + sett.fsig
            dplot = read(plotwritekey)
            sett.m(0, '--> using result', plotwritekey, 'for plotting')
            plotargs = [dtool, ddata, dplot]
        else: 
            plotargs = [dtool, ddata]
        if args['prev'] != '':
            plotargs.append(dprev)
        plot(*tuple(plotargs), **plotparams)

def read_args_run_tool(toolkey):
    """
    Read arguments and run tool specified by toolkey.
    """
    args = utils.read_args_tool(toolkey, exs.dexamples())
    return run_args(toolkey, args)

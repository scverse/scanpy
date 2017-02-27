# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Scanpy - Single-Cell Analysis in Python

Reference
---------
Wolf, Angerer & Theis, bioRxiv doi:... (2017)
"""

import numpy as np
import os

from . import settings as sett
from . import tools
from . import preprocess
from . import utils
from .tools import get_tool
from .classes.ann_data import AnnData
from .readwrite import read, write, read_params
from .examples import show_exdata, show_examples, get_example
from . import preprocess
from .preprocess.advanced import subsample
from .tools.diffmap import diffmap, plot_diffmap
from .tools.tsne import tsne, plot_tsne
from .tools.dpt import dpt, plot_dpt
from .tools.pca import pca, plot_pca
from .tools.diffrank import diffrank, plot_diffrank
from .tools.sim import sim, plot_sim

# just an equivalent name
pp = preprocess

__all__ = [
    # example use cases
    'get_example', # call example
    'show_exdata', # show available example data
    'show_examples', # show available example use cases
    # help
    'help', # show help for a given tool
    # elementary operations
    'read',
    'write',
    # preprocessing
    'preprocess', 'pp',
    'subsample',
    # visualization
    'diffmap', 'plot_diffmap',
    'tsne', 'plot_tsne',
    'pca', 'plot_pca'
    # subgroup identification
    'dpt', 'plot_dpt'
    # differential expression testing
    'diffrank', 'plot_diffrank'
    # simulation
    'sim', 'plot_sim'
    # plotting
    'show',
    # classes
    'AnnData'
]

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
    from .compat.matplotlib import pyplot as pl
    pl.show()

def _run_command_line_args(toolkey, args):
    """
    Run specified tool, do preprocessing and read/write outfiles.

    Result files store the dictionary returned by the tool. File type is
    determined by variable sett.extd allowed are 'h5' (hdf5), 'xlsx' (Excel) or
    'csv' (comma separated value file).

    If called twice with the same settings the existing result file is used.

    Parameters
    ----------
    toolkey : str
        Name of the tool.
    args : dict containing
        exkey : str
            String that identifies the example use key.
    """
    # help on plot parameters
    if args['pparams']:
        if args['pparams'][0] == 'help':
            from sys import exit
            exit(get_tool(toolkey).plot.__doc__)

    # read parameters
    adata = None
    if os.path.exists(readwrite.get_filename_from_key(sett.basekey)) or toolkey != 'sim':
        adata, exmodule = get_example(args['exkey'], subsample=args['subsample'],
                                      return_module=True)
        oparams = {}
        # try to load tool parameters from dexamples
        try:
            did_not_find_params_in_exmodule = False
            dexample = exmodule.dexamples[args['exkey']]
            oparams = {}
            for key in dexample.keys():
                if toolkey in key:
                    oparams = dexample[key]
                    sett.m(0, '... appending "-o',
                           ' '.join([' '.join([k, str(v)]) for k, v in oparams.items()])
                           + '"',
                          'to call of', toolkey)
                    break
        except:
            did_not_find_params_in_exmodule = True
            pass
        # if optional parameters have been specified in a parameter file update
        # the current param dict with these
        if args['opfile'] != '':
            add_params = read_params(args['opfile'])
            oparams = utils.update_params(oparams, add_params)
        # same if optional parameters have been specified on the command line
        if args['oparams']:
            add_params = readwrite.get_params_from_list(args['oparams'])
            sett.m(0, '... overwriting optional params', '"' + 
                   ' '.join([' '.join([k, str(v)]) for k, v in add_params.items()])
                   + '"',
                  'to call of', toolkey)
            oparams = utils.update_params(oparams, add_params)
        elif did_not_find_params_in_exmodule and args['opfile'] != '':
            sett.m(0, 'using default parameters, change them using "--oparams"')
    elif toolkey == 'sim':
        if args['opfile'] != '':
            oparams = read_params(args['opfile'])
        else:
            from . import sim_models
            opfile_sim = os.path.dirname(sim_models.__file__) + '/' + args['exkey'] + '_oparams.txt'
            oparams = read_params(opfile_sim)
            sett.m(0,'--> you can specify your custom params file using the option\n'
                     '    "--opfile" or provide parameters directly via "--oparams"')
        if 'writedir' not in oparams:
            oparams['writedir'] = sett.writedir + sett.basekey + '_' + toolkey


    # read/write files
    writekey = sett.basekey + '_' + toolkey
    opfile = sett.writedir + writekey + '_oparams.txt'
    if args['logfile']:
        logfile = sett.writedir + writekey + '_log.txt'
        sett.logfile(logfile)
    
    # actual call of tool
    if (adata is None
        or toolkey not in adata['tools']
        or sett.recompute != 'none'):
        tool = get_tool(toolkey, func=True)
        if toolkey == 'sim':
            adata = tool(**oparams)
        else:
            adata = tool(adata, **oparams)
        # append toolkey to tools in adata
        if toolkey not in adata['tools']:
            adata['tools'] = np.append(adata['tools'], toolkey)
        write(sett.basekey, adata)
        sett.m(0, 'updated file',
               readwrite.get_filename_from_key(sett.basekey))
        # save a copy of the changed parameters
        readwrite.write_params(opfile, oparams)

    # plotting and postprocessing
    pparams = (readwrite.get_params_from_list(args['pparams']) 
               if args['pparams'] else {})
    # post-processing specific to example and tool
    # - only if we are not subsampling
    if toolkey != 'sim':
        postprocess = args['exkey'] + '_' + toolkey
        if postprocess in dir(exmodule) and args['subsample'] == 1:
            adata = getattr(exmodule, postprocess)(adata)
            write(sett.basekey, adata)
    getattr(get_tool(toolkey), 'plot_' + toolkey)(adata, **pparams)
   
def _read_command_line_args_run_single_tool(toolkey):
    """
    Read arguments and run tool specified by toolkey.
    """
    args = utils.read_args_tool(toolkey, examples.dexamples())
    run_args(toolkey, args)

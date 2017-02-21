# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Scanpy - Single-Cell Analysis in Python

Reference
---------
Wolf, Angerer & Theis, bioRxiv doi:... (2017)
"""

from . import settings as sett
from . import tools
from . import preprocess
from . import utils
from .tools import get_tool
from .ann_data import AnnData
from .readwrite import read, write, read_params
from .exs import exdata, examples, example
from . import preprocess
from .preprocess.advanced import subsample
from .tools.diffmap import diffmap
from .tools.tsne import tsne
from .tools.dpt import dpt
from .tools.pca import pca
from .tools.difftest import difftest
from .tools.sim import sim

# just an equivalent name
pp = preprocess

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
    'pca',
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
    # classes
    'AnnData'
]

def plot(toolkey, adata, **pparams):
    """
    Plot the result of a tool computation.
    """
    get_tool(toolkey).plot(adata, **pparams)

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

def run_args(toolkey, args):
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
    if toolkey == 'sim':
        if args['opfile'] != '':
            oparams = read_params(args['opfile'])
        else:
            opfile_sim = 'sim/' + args['exkey'] + '_oparams.txt'
            oparams = read_params(opfile_sim)
            sett.m(0,'--> you can specify your custom params file using the option\n'
                     '    "--opfile" or provide parameters directly via "--oparams"')
        if 'writedir' not in oparams:
            oparams['writedir'] = sett.writedir + sett.basekey + '_' + toolkey
    else:
        adata, exmodule = example(args['exkey'], subsample=args['subsample'],
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

    # read/write files
    writekey = sett.basekey + '_' + toolkey
    opfile = sett.writedir + writekey + '_oparams.txt'
    if args['logfile']:
        logfile = sett.writedir + writekey + '_log.txt'
        sett.logfile(logfile)
    
    # actual call of tool
    from os.path import exists
    if ((toolkey == 'sim'
         or toolkey not in adata['tools'])
        or sett.recompute != 'none'):
        tool = get_tool(toolkey, func=True)
        if toolkey == 'sim':
            adata = tool(**oparams)
        else:
            adata = tool(adata, **oparams)
        # append toolkey to tools in adata
        if (toolkey != 'sim' 
            and toolkey not in adata['tools']):
            import numpy as np
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
    plot(toolkey, adata, **pparams)
   
def read_args_run_tool(toolkey):
    """
    Read arguments and run tool specified by toolkey.
    """
    args = utils.read_args_tool(toolkey, exs.dexamples())
    run_args(toolkey, args)


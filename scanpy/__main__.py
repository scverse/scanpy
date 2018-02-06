# Author: F. Alex Wolf (http://falexwolf.de)
"""Scanpy - Single-Cell Analysis in Python

This is the command-line interface.
"""

import os
import sys
import argparse
from collections import OrderedDict as odict
from sys import argv, exit

# preprocessing
dpp = odict([
    ('pp', 'run a preprocessing function'),
])

# description of simple inquiries
dsimple = odict([
    ('exdata', 'show example data'),
])

# description of standard tools
dtools = odict([
    ('pca', 'visualize data using PCA'),
    ('diffmap', 'visualize data using Diffusion Map'''),
    ('tsne', 'visualize data using tSNE'),
    ('draw_graph', 'visualize data using force-directed graph drawing'),
    ('dpt', 'perform Diffusion Pseudotime analysis'),
    ('louvain', 'cluster cells using the louvain algorithm'),
    ('rank_genes_groups', ' test groups for differential expression'),
    ('score', ' calculates cell scores according to given genes'),
    ('sim', 'simulate stochastic gene expression models'),
])


# assemble main description
def main_descr():
    string = '\nsimple inquiries\n----------------'
    for key, descr in dpp.items():
        string += '\n{:12}'.format(key) + descr
    for key, descr in dsimple.items():
        string += '\n{:12}'.format(key) + descr
    for key, descr in dtools.items():
        if key == 'pca':
            string += '\n\nplotting tools\n--------------'
        if key == 'dpt':
            string += '\n\nother tools\n-----------'
        string += '\n{:12}'.format(key) + descr
    string += '\n\nNote: For convenience, you can swap the first argument (run_name) of any tool and the tool name.'
    return string


def add_args(p, dadd_args=None):
    """Add arguments to parser.

    Parameters
    -------
    dadd_args : dict
        Dictionary of additional arguments formatted as
            {'arg': {'type': int, 'default': 0, ... }}
    """
    aa = p.add_argument_group('Look up a preprocessing function').add_argument
    aa('run_name',
       type=str, default='', metavar='run_name',
       help='The "run_name" is used to read and preprocess the data and, if provided, default parameters of tools. '
            'It\'s the name of a function that returns an annotated data object and is stored in a user module "preprocessing_whatevername.py" in the current working directory. '
            'There are many builtin examples with such corresponding functions. '
            'All output files start with the run_name. '
            'Use the subcommands "exdata" and "exparams" to display the corresponding data and parameters, respectively.')
    aa = p.add_argument_group('Provide tool parameters').add_argument
    # example key default argument
    aa('-p', '--params',
       nargs='*', default=None, metavar='k=v',
       help='Optional tool parameters as list, e.g., `k=20 knn=True`, '
            'see detailed help below, also notice the option --opfile (default: "").')
    # arguments from dadd_args
    if dadd_args is not None:
        for key, val in dadd_args.items():
            if key != 'arg':
                aa(key, **val)
    # make sure there are is no conflict with dadd_args
    if dadd_args is None or '--pfile' not in dadd_args:
        aa('--pfile',
           type=str, default='', metavar='f',
           help='Path to file with parameters (default: "").')

    aa = p.add_argument_group('Provide plotting parameters').add_argument
    aa('-q', '--plot_params',
       nargs='*', default=None, metavar='k v',
       help='Plotting parameters as list, e.g., `comps=1,3 legendloc="upper left"`. '
            'Display possible paramaters via "-p help" (default: "").')

    # settings-related arguments
    from .settings import add_args
    p = add_args(p)

    return p


def init_main_parser():
    """Add sub parsers for each tool to main parser.
    """
    # the actual parser and parser container
    main_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)
    sub_parsers = main_parser.add_subparsers(metavar='',
                                             description=main_descr())

    from .api import tools
    for key in dtools:
        # use the doc string of the tool for outputting the description
        descr = 78*'-' + '\n' + getattr(tools, key).__doc__
        sub_p = sub_parsers.add_parser(
            key,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False,
            epilog=descr)
        try:
            sub_p = getattr(tools, key).add_args(sub_p)
        except AttributeError:
            sub_p = add_args(sub_p)
        sub_p.set_defaults(key=key)

    # preprocessing
    for key in dpp:
        # use the doc string of the tool for outputting the description
        descr = 78*'-' + '\n'
        sub_p = sub_parsers.add_parser(
            key,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False,
            epilog=descr)
        try:
            sub_p = getattr(tools, key).add_args(sub_p)
        except AttributeError:
            sub_p = add_args(sub_p)
        sub_p.set_defaults(key=key)

    return main_parser


def run_command_line_args(toolkey, args):
    """Run specified tool, do preprocessing and read/write outfiles.

    Result files store the dictionary returned by the tool. File type is
    determined by variable settings.file_format_data allowed are 'h5' (hdf5), 'xlsx' (Excel) or
    'csv' (comma separated value file).

    If called twice with the same settings the existing result file is used.

    Parameters
    ----------
    toolkey : str
        Name of the tool or 'pp' for simply calling preprocessing.
    args : dict containing at least
        run_name : str
            String that identifies the example use key.
    """
    from . import settings
    from . import logging as logg
    from . import readwrite

    # help on plot parameters
    if args['plot_params']:
        from . import plotting
        if args['plot_params'][0] == 'help':
            sys.exit(getattr(plotting, toolkey).__doc__)

    adata = None
    exmodule = None
    if toolkey != 'sim':
        from .examples import init_run
        adata, exmodule = init_run(
            args['run_name'],
            suffix=settings._run_suffix,
            return_module=True,
            recompute=settings.recompute == 'pp' or toolkey == 'pp',
            reread=settings.recompute == 'read')
        params = {}
        # try to load tool parameters from dexamples
        did_not_find_params_in_exmodule = False
        try:
            params = getattr(exmodule, args['run_name'] + '_' + toolkey + '_params')
            logg.m('... loading run params "-p {}"'.format(
                   ' '.join(['='.join([k, str(v)]) for k, v in params.items()])))
        except AttributeError:
            did_not_find_params_in_exmodule = True
        # if optional parameters have been specified in a parameter file update
        # the current param dict with these
        from . import utils
        if args['pfile'] != '':
            add_params = readwrite.read_params(args['pfile'])
            params = utils.update_params(params, add_params)
        # same if optional parameters have been specified on the command line
        if args['params']:
            add_params = readwrite.get_params_from_list(args['params'])
            logg.m('... overwriting params', '"' +
                   ' '.join(['='.join([k, str(v)]) for k, v in add_params.items()])
                   + '"', 'in call of', toolkey)
            params = utils.update_params(params, add_params)
        elif did_not_find_params_in_exmodule and args['pfile'] != '':
            logg.m('using default parameters, change them using "--params"')
    elif toolkey == 'sim':
        if args['pfile'] != '':
            params = readwrite.read_params(args['pfile'])
        else:
            from . import sim_models
            from . import readwrite
            pfile_sim = os.path.dirname(sim_models.__file__) + '/' + args['run_name'] + '_params.txt'
            params = readwrite.read_params(pfile_sim)
            logg.m('--> you can specify your custom params file using the option\n'
                   '    "--pfile" or provide parameters directly via "--params"')
        if 'writedir' not in params:
            params['writedir'] = settings.writedir + settings.run_name + '_' + toolkey

    if toolkey == 'pp': exit()

    if (adata is not None and settings.recompute == 'none'):
        try:
            run_plotting_and_postprocessing(args, adata, toolkey, exmodule)
            return
        except KeyError:
            pass

    # actual call of tool
    from .api import tools
    tool = getattr(tools, toolkey)
    if toolkey == 'sim': adata = tool(**params)
    else: tool(adata, **params)
    readwrite.write(settings.run_name, adata)
    if settings.file_format_data not in {'h5', 'npz'}:
        readwrite.write(settings.run_name, adata, ext='h5')
    # save a copy of the changed parameters
    pfile = settings.writedir + settings.run_name + '_params/' + toolkey + '.txt'
    readwrite.write_params(pfile, params)
    run_plotting_and_postprocessing(args, adata, toolkey, exmodule)


def run_plotting_and_postprocessing(args, adata, toolkey, exmodule):
    from . import readwrite
    from . import settings
    plot_params = (readwrite.get_params_from_list(args['plot_params'])
               if args['plot_params'] else {})
    # post-processing specific to example and tool
    # - only if we are not subsampling
    if exmodule is not None:
        postprocess = args['run_name'] + '_' + toolkey
        if postprocess in dir(exmodule) and args['subsample'] == 1:
            getattr(exmodule, postprocess)(adata)
            readwrite.write(settings.run_name, adata)
    from . import plotting
    getattr(plotting, toolkey)(adata, **plot_params)


def read_command_line_args_run_single_tool(toolkey):
    """Read arguments and run tool specified by toolkey.
    """
    sys.exit('`read_command_line_args_run_single_tool()` is currently not supported')
    args = utils.read_args_tool(toolkey, examples._example_parameters())
    run_command_line_args(toolkey, args)


def main():
    # check whether at least one subcommand has been supplied
    if len(argv) == 1 or argv[1] == 'h' or argv[1] == '--help':
        init_main_parser().print_help()
        exit(0)
    # simple inquiries
    if argv[1] in dsimple:
        print('See https://github.com/theislab/scanpy/blob/master/scanpy/examples/builtin.py')
        print('for builtin examples. Function names serve as "run_names" '
              '(keys for looking up preprocessing and default parameters).')
        exit(0)
    # init the parsers for each tool
    main_parser = init_main_parser()
    # test whether run_name is provided first
    valid_keys = set(list(dtools.keys()) + list(dpp.keys()))
    if argv[1] not in valid_keys:
        if len(argv) > 2 and argv[2] in valid_keys:
            run_name = argv[1]
            argv[1] = argv[2]
            argv[2] = run_name
        else:
            print('normal usage:    ' + argv[0] + ' tool run_name')
            print('efficient usage: ' + argv[0] + ' run_name tool')
            print('help:            ' + argv[0] + ' -h')
            exit(0)
    # transform to dictionary
    args = vars(main_parser.parse_args(argv[1:]))
    from .settings import process_args
    args = process_args(args)
    # run Scanpy
    run_command_line_args(args['key'], args)


# enable `python -m scanpy ...` call in addition to callin the
# entrypoint `main()` via `scanpy ...`
if __name__ == '__main__':
    main()

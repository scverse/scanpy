# Author: F. Alex Wolf (http://falexwolf.de)
"""
Scanpy - Single-Cell Analysis in Python

This is the general purpose command-line utility.
"""

# this is necessary to import scanpy from within package
from __future__ import absolute_import
import argparse
from collections import OrderedDict as odict
from sys import argv, exit
import scanpy as sc

# description of simple inquiries
dsimple = odict([
    ('exdata', 'show example data'),
    ('examples', 'show example use cases'),
])

# description of standard tools
dtools = odict([
    ('pca', 'visualize data using PCA'),
    ('diffmap', 'visualize data using Diffusion Map'''),
    ('tsne', 'visualize data using tSNE'),
    ('spring', 'visualize data using force-directed graph drawing'),
    ('dpt', 'perform Diffusion Pseudotime analysis'),
    ('diffrank', 'test for differential expression'),
    ('sim', 'simulate stochastic gene expression models'),
])

# assemble main description
def main_descr():
    descr = '\nsimple inquiries\n----------------'
    for key, help in dsimple.items():
        descr += '\n{:12}'.format(key) + help
    for key, help in dtools.items():
        if key == 'pca':
            descr += '\n\nplotting tools\n--------------'
        if key == 'dpt':
            descr += '\n\nother tools\n-----------'
        descr += '\n{:12}'.format(key) + help
    descr += '\n\nexkey tool\n----------'
    descr += ('\n{:12}'.format('exkey tool')
                   + 'shortcut for providing exkey argument to tool')
    return descr

def init_main_parser():
    """
    Init subparser for each tool.
    """
    # the actual parser and parser container
    main_parser = argparse.ArgumentParser(
                      description=__doc__,
                      formatter_class=argparse.RawDescriptionHelpFormatter,
                      add_help=False)
    sub_parsers = main_parser.add_subparsers(metavar='',
                                             description=main_descr())

    for key, help in dtools.items():
        descr = 78*'-' + '\n' + sc.help(key, string=True) + 78*'-'
        sub_p = sub_parsers.add_parser(
                    key,
                    description=descr,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    add_help=False)
        try:
            sub_p = sc.get_tool(key).add_args(sub_p)
        except:
            sub_p = sc.utils.add_args(sub_p)
        sub_p.set_defaults(toolkey=key)
    
    return main_parser

def main():
    # check whether at least one subcommand has been supplied
    if len(argv) == 1 or argv[1] == 'h' or argv[1] == '--help':
        init_main_parser().print_help()
        exit(0)
    # simple inquiries
    if argv[1] in dsimple:
        # same keys as in dsimple
        func = {
            'exdata': sc.show_exdata,
            'examples': sc.show_examples
        }
        if len(argv) > 2:
            func[argv[1]](argv[2])
        else:
            func[argv[1]]()
        exit(0)
    # init the parsers for each tool
    main_parser = init_main_parser()
    # test whether exkey is provided first
    if argv[1] not in dtools:
        if len(argv) > 2 and argv[2] in dtools:
            exkey = argv[1]
            argv[1] = argv[2]
            argv[2] = exkey
        else:
            print('normal usage:    ' + argv[0] + ' tool exkey')
            print('efficient usage: ' + argv[0] + ' exkey tool')
            print('help:            ' + argv[0] + ' -h')
            exit(0)
    # transform to dictionary
    args = vars(main_parser.parse_args(argv[1:]))
    args = sc.sett.process_args(args)
    # run Scanpy
    sc._run_command_line_args(args['toolkey'], args)

if __name__ == '__main__':
    main()

from sphinx.ext.napoleon import NumpyDocstring, GoogleDocstring
from sphinx.ext.napoleon import Config
import inspect
import scanpy as sc


def iterate_over_functions():
    functions_all = {}
    for mn in dir(sc):
        # skip privates and duplicates and stuff we haven't formatted in the righ way, yet
        if mn.startswith('_') or mn in {'tl', 'pl', 'pp', 'readwrite', 'utils', 'logging', 'neighbors'}:
            continue
        module = sc.__dict__[mn]
        if not inspect.ismodule(module):
            continue
        function_names = [fn for fn in dir(module) if inspect.isfunction(module.__dict__[fn])]
        functions = {mn + '.' + fn: module.__dict__[fn] for fn in function_names}
        functions_all.update(functions)
    return functions_all


def test_function_headers():
    for descr, f in iterate_over_functions().items():
        lines = f.__doc__.split('\n')
        if lines[0].startswith('    ') or lines[0] == '':
            raise Exception(
                'Header of docstring of function `{}` should start with one-line description:\n'
                '"""My one-line description.\n'
                'not\n{}'.format(descr, lines[:2]))

import inspect
import scanpy as sc


def iterate_over_functions():
    functions_all = {}
    for mn in dir(sc):
        # skip privates and duplicates and stuff we haven't formatted in the righ way, yet
        if mn.startswith('_') or mn in {'tl', 'pl', 'pp', 'cli', 'readwrite', 'utils', 'logging', 'neighbors'}:
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
        assert f.__doc__ is not None, '{} has no docstring'.format(f)
        lines = f.__doc__.split('\n')
        # assert lines[0] and not lines[0].startswith('    '), (
        #     'Header of docstring of function `{}` should start with one-line description:\n'
        #     '"""My one-line description.\n'
        #     'not\n{}'.format(descr, lines[:2])
        # )

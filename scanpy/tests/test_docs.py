import inspect
import traceback

import pytest

import scanpy as sc


def iterate_over_functions():
    functions_all = {}
    for mn in dir(sc):
        # skip privates and duplicates and stuff we haven’t formatted
        # in the right way, yet
        if mn.startswith('_') or mn in {
            'tl', 'pl', 'pp',
            'cli', 'readwrite', 'utils', 'logging', 'neighbors'
        }:
            continue
        module = sc.__dict__[mn]
        if not inspect.ismodule(module):
            continue
        function_names = [fn for fn in dir(module) if inspect.isfunction(module.__dict__[fn])]
        functions = {mn + '.' + fn: module.__dict__[fn] for fn in function_names}
        functions_all.update(functions)
    return functions_all


@pytest.mark.parametrize("descr,f", iterate_over_functions().items())
def test_function_headers(descr, f):
    assert f.__doc__ is not None, f"{f} has no docstring"
    lines = getattr(f, "__orig_doc__", f.__doc__).split("\n")
    assert lines[0], "Expected single-line summary"
    broken = [
        i for i, l in enumerate(lines)
        if l and not l.startswith("    ")
    ]
    if any(broken):
        msg = f'''\
Header of function `{descr}`’s docstring should start with one-line description:

␣␣␣␣"""\\
␣␣␣␣My one-line␣description.

␣␣␣␣…
␣␣␣␣"""
'''
        filename = inspect.getsourcefile(f)
        _, lineno = inspect.getsourcelines(f)
        text = f">{lines[broken[0]]}<"
        raise SyntaxError(msg, (filename, lineno, 2, text))

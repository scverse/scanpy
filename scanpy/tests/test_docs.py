import inspect
from types import FunctionType

import pytest
from scanpy._utils import descend_classes_and_funcs

# CLI is locally not imported by default but on travis it is?
import scanpy.cli


scanpy_functions = [
    c_or_f
    for c_or_f in descend_classes_and_funcs(scanpy, "scanpy")
    if isinstance(c_or_f, FunctionType)
]


@pytest.mark.parametrize("f", scanpy_functions)
def test_function_headers(f):
    name = f"{f.__module__}.{f.__qualname__}"
    assert f.__doc__ is not None, f"{name} has no docstring"
    lines = getattr(f, "__orig_doc__", f.__doc__).split("\n")
    assert lines[0], f"{name} needs a single-line summary"
    broken = [i for i, l in enumerate(lines) if l and not l.startswith("    ")]
    if any(broken):
        msg = f'''\
Header of function `{name}`’s docstring should start with one-line description:

␣␣␣␣"""\\
␣␣␣␣My one-line␣description.

␣␣␣␣…
␣␣␣␣"""
'''
        filename = inspect.getsourcefile(f)
        _, lineno = inspect.getsourcelines(f)
        text = f">{lines[broken[0]]}<"
        raise SyntaxError(msg, (filename, lineno, 2, text))


def test_plot_doc_signatures():
    pass

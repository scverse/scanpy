import email
import inspect
import os
from importlib.util import find_spec
from types import FunctionType
from pathlib import Path

import pytest
from scanpy._utils import descend_classes_and_funcs

# CLI is locally not imported by default but on travis it is?
import scanpy.cli


mod_dir = Path(scanpy.__file__).parent
proj_dir = mod_dir.parent

scanpy_functions = [
    c_or_f
    for c_or_f in descend_classes_and_funcs(scanpy, "scanpy")
    if isinstance(c_or_f, FunctionType)
]


@pytest.fixture
def in_project_dir():
    wd_orig = Path.cwd()
    os.chdir(proj_dir)
    try:
        yield proj_dir
    finally:
        os.chdir(wd_orig)


@pytest.mark.parametrize("f", scanpy_functions)
def test_function_headers(f):
    name = f"{f.__module__}.{f.__qualname__}"
    assert f.__doc__ is not None, f"{name} has no docstring"
    lines = getattr(f, "__orig_doc__", f.__doc__).split("\n")
    broken = [i for i, l in enumerate(lines) if l.strip() and not l.startswith("    ")]
    if any(broken):
        msg = f'''\
Header of function `{name}`’s docstring should start with one-line description
and be consistently indented like this:

␣␣␣␣"""\\
␣␣␣␣My one-line␣description.

␣␣␣␣…
␣␣␣␣"""

The displayed line is under-indented.
'''
        filename = inspect.getsourcefile(f)
        _, lineno = inspect.getsourcelines(f)
        text = f">{lines[broken[0]]}<"
        raise SyntaxError(msg, (filename, lineno, 2, text))

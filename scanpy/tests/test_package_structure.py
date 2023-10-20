import os
from inspect import Parameter, signature
from types import FunctionType, ModuleType
from pathlib import Path

import pytest
from scanpy._utils import descend_classes_and_funcs

# CLI is locally not imported by default but on travis it is?
import scanpy.cli


mod_dir = Path(scanpy.__file__).parent
proj_dir = mod_dir.parent


def get_name(func: FunctionType) -> str:
    # TODO: more declarative
    mod = func.__module__
    if mod.startswith("scanpy.readwrite"):
        mod = "sc"
    elif mod.startswith("scanpy.datasets"):
        mod = "sc.datasets"
    elif mod.startswith("scanpy.queries"):
        mod = "sc.queries"
    elif mod.startswith("scanpy.metrics"):
        mod = "sc.metrics"
    elif mod.startswith("scanpy.get"):
        mod = "sc.get"
    elif mod.startswith("scanpy.neighbors"):
        mod = "sc.neighbors"
    elif mod.startswith("scanpy.logging"):
        mod = "sc.logging"
    elif mod.startswith("scanpy.preprocessing"):
        mod = "sc.pp"
    elif mod.startswith("scanpy.tools"):
        mod = "sc.tl"
    elif mod.startswith("scanpy.plotting"):
        mod = "sc.pl"
    elif mod.startswith("scanpy.experimental.pp"):
        mod = "sc.experimental.pp"
    elif mod.startswith("scanpy.external.pp"):
        mod = "sc.external.pp"
    elif mod.startswith("scanpy.external.tl"):
        mod = "sc.external.tl"
    elif mod.startswith("scanpy.external.pl"):
        mod = "sc.external.pl"
    elif mod.startswith("scanpy.external.exporting"):
        mod = "sc.external.exporting"
    return f"{mod}.{func.__name__}"


scanpy_functions = [
    pytest.param(func, id=get_name(func))
    for func in sorted(descend_classes_and_funcs(scanpy, "scanpy"), key=get_name)
    if isinstance(func, FunctionType)
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
    filename = getsourcefile(f)
    lines, lineno = getsourcelines(f)
    if f.__doc__ is None:
        msg = f"Function `{name}` has no docstring"
        text = lines[0]
    else:
        lines = getattr(f, "__orig_doc__", f.__doc__).split("\n")
        broken = [
            i for i, l in enumerate(lines) if l.strip() and not l.startswith("    ")
        ]
        if not any(broken):
            return
        msg = f'''\
Header of function `{name}`’s docstring should start with one-line description
and be consistently indented like this:

␣␣␣␣"""\\
␣␣␣␣My one-line␣description.

␣␣␣␣…
␣␣␣␣"""

The displayed line is under-indented.
'''
        text = f">{lines[broken[0]]}<"
    raise SyntaxError(msg, (filename, lineno, 2, text))


@pytest.mark.parametrize("f", scanpy_functions)
def test_function_positional_args(f):
    """See https://github.com/astral-sh/ruff/issues/3269#issuecomment-1772632200"""
    sig = signature(f)
    pos_kinds = {
        Parameter.POSITIONAL_ONLY,
        Parameter.POSITIONAL_OR_KEYWORD,
    }
    n_pos = sum(1 for p in sig.parameters.values() if p.kind in pos_kinds)
    n_pos_max = 5
    if n_pos <= n_pos_max:
        return

    msg = f"Function `{f.__module__}.{f.__qualname__}` has too many positional arguments ({n_pos}>{n_pos_max})"
    filename = getsourcefile(f)
    lines, lineno = getsourcelines(f)
    text = lines[0]
    raise SyntaxError(msg, (filename, lineno, 1, text))


def getsourcefile(obj):
    """inspect.getsourcefile, but supports singledispatch"""
    from inspect import getsourcefile

    if wrapped := getattr(obj, "__wrapped__", None):
        return getsourcefile(wrapped)

    return getsourcefile(obj)


def getsourcelines(obj):
    """inspect.getsourcelines, but supports singledispatch"""
    from inspect import getsourcelines

    if wrapped := getattr(obj, "__wrapped__", None):
        return getsourcelines(wrapped)

    return getsourcelines(obj)

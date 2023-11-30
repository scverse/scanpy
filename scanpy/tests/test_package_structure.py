from __future__ import annotations

import os
from inspect import Parameter, signature
from pathlib import Path

import pytest

# CLI is locally not imported by default but on travis it is?
import scanpy.cli
from scanpy._utils import _import_name, descend_classes_and_funcs

mod_dir = Path(scanpy.__file__).parent
proj_dir = mod_dir.parent


api_module_names = [
    "sc",
    "sc.pp",
    "sc.tl",
    "sc.pl",
    "sc.experimental.pp",
    "sc.external.pp",
    "sc.external.tl",
    "sc.external.pl",
    "sc.external.exporting",
    "sc.get",
    "sc.logging",
    # "sc.neighbors",  # Not documented
    "sc.datasets",
    "sc.queries",
    "sc.metrics",
]
api_modules = {
    mod_name: _import_name(f"scanpy{mod_name.removeprefix('sc')}")
    for mod_name in api_module_names
}


api_functions = [
    pytest.param(func, id=f"{mod_name}.{name}")
    for mod_name, mod in api_modules.items()
    for name in mod.__all__
    if callable(func := getattr(mod, name))
]


@pytest.fixture
def in_project_dir():
    wd_orig = Path.cwd()
    os.chdir(proj_dir)
    try:
        yield proj_dir
    finally:
        os.chdir(wd_orig)


def test_descend_classes_and_funcs():
    # TODO: unclear if we want this to totally match, let’s see
    funcs = set(descend_classes_and_funcs(scanpy, "scanpy"))
    assert {p.values[0] for p in api_functions} == funcs


@pytest.mark.parametrize("f", api_functions)
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


@pytest.mark.parametrize("f", api_functions)
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

"""Extension to warn about numpydoc-style parameter types in docstrings."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from sphinx.ext.napoleon import NumpyDocstring

if TYPE_CHECKING:
    from sphinx.application import Sphinx

_format_docutils_params_orig = NumpyDocstring._format_docutils_params
param_warnings = {}


def scanpy_log_param_types(self, fields, field_role="param", type_role="type"):
    """Wrap ``NumpyDocstring._format_docutils_params``."""
    for _name, _type, _desc in fields:
        if not _type or not self._obj.__module__.startswith("scanpy"):
            continue
        w_list = param_warnings.setdefault((self._name, self._obj), [])
        if (_name, _type) not in w_list:
            w_list.append((_name, _type))
    return _format_docutils_params_orig(self, fields, field_role, type_role)


def show_param_warnings(app, exception):
    """Warn about numpydoc-style parameter types in docstring."""
    import inspect

    for (fname, fun), params in param_warnings.items():
        _, line = inspect.getsourcelines(fun)
        file_name = inspect.getsourcefile(fun)
        assert file_name is not None
        params_str = "\n".join(f"\t{n}: {t}" for n, t in params)
        warnings.warn_explicit(
            f"\nParameters in `{fname}` have types in docstring.\n"
            f"Replace them with type annotations.\n{params_str}",
            UserWarning,
            file_name,
            line,
        )
    if param_warnings:
        msg = "Encountered text parameter type. Use annotations."
        raise RuntimeError(msg)


def setup(app: Sphinx):
    """App setup hook."""
    NumpyDocstring._format_docutils_params = scanpy_log_param_types
    app.connect("build-finished", show_param_warnings)

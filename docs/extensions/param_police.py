import warnings

from sphinx.application import Sphinx
from sphinx.ext.napoleon import NumpyDocstring


_format_docutils_params_orig = NumpyDocstring._format_docutils_params
param_warnings = {}


def scanpy_log_param_types(self, fields, field_role='param', type_role='type'):
    for _name, _type, _desc in fields:
        if not _type or not self._obj.__module__.startswith('scanpy'):
            continue
        w_list = param_warnings.setdefault((self._name, self._obj), [])
        if (_name, _type) not in w_list:
            w_list.append((_name, _type))
    return _format_docutils_params_orig(self, fields, field_role, type_role)


def show_param_warnings(app, exception):
    import inspect

    for (fname, fun), params in param_warnings.items():
        _, line = inspect.getsourcelines(fun)
        file_name = inspect.getsourcefile(fun)
        params_str = '\n'.join(f'\t{n}: {t}' for n, t in params)
        warnings.warn_explicit(
            f'\nParameters in `{fname}` have types in docstring.\n'
            f'Replace them with type annotations.\n{params_str}',
            UserWarning,
            file_name,
            line,
        )
    if param_warnings:
        raise RuntimeError('Encountered text parameter type. Use annotations.')


def setup(app: Sphinx):
    NumpyDocstring._format_docutils_params = scanpy_log_param_types
    app.connect('build-finished', show_param_warnings)

"""Logging and Profiling
"""

import time as time_module
import datetime
from typing import Union

from anndata import logging
from ._settings import settings


_VERBOSITY_LEVELS_FROM_STRINGS = {
    'error': 0,
    'warn': 1,
    'info': 2,
    'hint': 3,
    'debug': 4,
}


def info(*args, **kwargs):
    return _msg(*args, v='info', **kwargs)


def error(*args, **kwargs):
    args = ('Error:',) + args
    return _msg(*args, v='error', **kwargs)


def warn(*args, **kwargs):
    args = ('WARNING:',) + args
    return _msg(*args, v='warn', **kwargs)


def hint(*args, **kwargs):
    return _msg(*args, v='hint', **kwargs)


def debug(*args, **kwargs):
    return _msg(*args, v='debug', **kwargs)


def _settings_verbosity_greater_or_equal_than(v):
    if isinstance(settings.verbosity, str):
        settings_v = _VERBOSITY_LEVELS_FROM_STRINGS[settings.verbosity]
    else:
        settings_v = settings.verbosity
    return settings_v >= v


def _msg(
    *msg,
    v: Union[str, int],
    time: bool = False,
    reset: bool = False,
    end: str = '\n',
    no_indent: bool = False,
    t: bool = None,
    r: bool = None,
):
    """Write message to logging output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.

    v : {'error', 'warn', 'info', 'hint'} or int
        0/'error', 1/'warn', 2/'info', 3/'hint', 4/'debug'
    time, t
        Print timing information; restart the clock.
    reset, r
        Reset timing and memory measurement. Is automatically reset
        when passing one of ``time`` or ``memory``.
    end
        Same meaning as in builtin ``print()`` function.
    no_indent
        Do not indent for ``v >= 4``.
    """
    # variable shortcuts
    if t is not None: time = t
    if r is not None: reset = r
    if isinstance(v, str):
        v = _VERBOSITY_LEVELS_FROM_STRINGS[v]
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v >= 4 and not no_indent:
        msg = ('   ',) + msg
    if _settings_verbosity_greater_or_equal_than(v):
        if not time and len(msg) > 0:
            _write_log(*msg, end=end)
        if reset:
            try:
                settings._previous_memory_usage, _ = get_memory_usage()
            except:
                pass
            settings._previous_time = time_module.time()
        if time:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(_sec_to_str(elapsed)),)
            _write_log(*msg, end=end)


def _write_log(*msg, end='\n'):
    """Write message to log output, ignoring the verbosity level.

    This is the most basic function.

    Parameters
    ----------
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    out = ' '.join(map(str, msg))
    if hasattr(settings.logfile, 'open'):
        with settings.logfile.open('a') as f:
            f.write(out + end)
    elif hasattr(settings.logfile, 'write'):
        settings.logfile.write(out + end)
    else:
        raise TypeError(f'settings.logfile of unknown type {type(settings.logfile)}')


def _sec_to_str(t):
    """Format time in seconds.

    Parameters
    ----------
    t : int
        Time in seconds.
    """
    from functools import reduce
    return "%d:%02d:%02d.%02d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*100,), 100, 60, 60])


print_memory_usage = logging.print_memory_usage


get_memory_usage = logging.get_memory_usage


def get_passed_time():
    now = time_module.time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def print_version_and_date():
    from . import __version__
    _write_log(f'Running Scanpy {__version__}, on {get_date_string()}.')


_DEPENDENCIES_NUMERICS = [
    'anndata',  # anndata actually shouldn't, but as long as it's in development
    'umap',
    'numpy',
    'scipy',
    'pandas',
    ('sklearn', 'scikit-learn'),
    'statsmodels',
    ('igraph', 'python-igraph'),
    'louvain',
]


_DEPENDENCIES_PLOTTING = ['matplotlib', 'seaborn']


def _print_versions_dependencies(dependencies):
    # this is not the same as the requirements!
    for mod in dependencies:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try:
            imp = __import__(mod_name)
            print('{}=={}'.format(mod_install, imp.__version__), end=' ')
        except (ImportError, AttributeError):
            pass
    print()


def print_versions():
    """Versions that might influence the numerical results.

    Matplotlib and Seaborn are excluded from this.
    """
    _print_versions_dependencies(['scanpy'] + _DEPENDENCIES_NUMERICS)


def get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

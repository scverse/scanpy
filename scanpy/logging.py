"""Logging and Profiling
"""

import time as time_module
import datetime
from anndata import logging
from . import settings


verbosity_levels_from_strings = {
    'error': 0,
    'warn': 1,
    'info': 2,
    'hint': 3,
}


def info(*args, **kwargs):
    return msg(*args, v='info', **kwargs)


def error(*args, **kwargs):
    args = ('Error:',) + args
    return msg(*args, v='error', **kwargs)


def warn(*args, **kwargs):
    args = ('WARNING:',) + args
    return msg(*args, v='warn', **kwargs)


def hint(*args, **kwargs):
    return msg(*args, v='hint', **kwargs)


def verbosity_greater_or_equal_than(v):
    if isinstance(v, str):
        v = verbosity_levels_from_strings[v]
    if isinstance(settings.verbosity, str):
        global_v = verbosity_levels_from_strings[settings.verbosity]
    else:
        global_v = settings.verbosity
    return global_v >= v


def msg(*msg, v=4, time=False, memory=False, reset=False, end='\n',
        no_indent=False, t=None, m=None, r=None):
    """Write message to logging output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.

    v : {'error', 'warn', 'info', 'hint'} or int, (default: 4)
        0/'error', 1/'warn', 2/'info', 3/'hint', 4, 5, 6...
    time, t : bool, optional (default: False)
        Print timing information; restart the clock.
    memory, m : bool, optional (default: Faulse)
        Print memory information.
    reset, r : bool, optional (default: False)
        Reset timing and memory measurement. Is automatically reset
        when passing one of ``time`` or ``memory``.
    end : str (default: '\n')
        Same meaning as in builtin ``print()`` function.
    no_indent : bool (default: False)
        Do not indent for ``v >= 4``.
    """
    # variable shortcuts
    if t is not None: time = t
    if m is not None: memory = m
    if r is not None: reset = r
    if isinstance(v, str):
        v = verbosity_levels_from_strings[v]
    if isinstance(settings.verbosity, str):
        global_verbosity = verbosity_levels_from_strings[settings.verbosity]
    else:
        global_verbosity = settings.verbosity
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v >= 4 and not no_indent:
        msg = ('   ',) + msg
    if global_verbosity >= v:
        if not time and not m and len(msg) > 0:
            settings.mi(*msg, end=end)
        if reset:
            try:
                settings._previous_memory_usage, _ = get_memory_usage()
            except:
                pass
            settings._previous_time = time_module.time()
        if time:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(settings._sec_to_str(elapsed)),)
            settings.mi(*msg, end=end)
        if memory:
            settings.mi(format_memory_usage(get_memory_usage()),
                        msg='' if time else msg, end=end)

m = msg
print_memory_usage = logging.print_memory_usage
get_memory_usage = logging.get_memory_usage

def get_passed_time():
    now = time_module.time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def print_version_and_date():
    from . import __version__
    settings.mi('Running Scanpy', __version__, 'on {}.'.format(get_date_string()))


_dependencies_numerics = [
    'anndata',  # anndata actually shouldn't, but as long as it's in development
    'numpy',
    'scipy',
    'pandas',
    ('sklearn', 'scikit-learn'),
    'statsmodels',
    ('igraph', 'python-igraph'),
    'louvain']


_dependencies_plotting = ['matplotlib', 'seaborn']


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
    _print_versions_dependencies(['scanpy'] + _dependencies_numerics)


def print_versions_dependencies_numerics():
    """Dependencies that might influence numerical results (computed data).
    """
    print('Dependencies:', end=' ')
    _print_versions_dependencies(_dependencies_numerics)


def print_versions_dependencies_plotting():
    """Dependencies that might influence plots (but not computed data).
    """
    print('Dependencies:', end=' ')
    _print_versions_dependencies(_dependencies_plotting)


def print_versions_dependencies_all():
    """All relevant dependencies.
    """
    _print_versions_dependencies(
        _dependencies_numerics + _dependencies_plotting)


def get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# Author: Alex Wolf (http://falexwolf.de)
"""Logging and Profiling
"""

import os
import time as time_module
import datetime
from . import settings


verbosity_levels_from_strings = {
    'error': 0,
    'warn': 1,
    'info': 2,
    'hint': 3,
}


def info(*msg, **kwargs):
    return m(*msg, v='info', **kwargs)


def error(*msg, **kwargs):
    msg = ('Error:',) + msg
    return m(*msg, v='error', **kwargs)


def warn(*msg, **kwargs):
    msg = ('WARNING:',) + msg
    return m(*msg, v='warn', **kwargs)


def hint(*msg, **kwargs):
    return m(*msg, v='hint', **kwargs)


def verbosity_greater_or_equal_than(v):
    if isinstance(v, str):
        v = verbosity_levels_from_strings[v]
    if isinstance(settings.verbosity, str):
        global_v = verbosity_levels_from_strings[settings.verbosity]
    else:
        global_v = settings.verbosity
    return global_v >= v


def msg(*msg, verbosity='info', time=False, memory=False, reset=False, end='\n',
        no_indent=False, v=None, t=None, m=None, r=None):
    """Write message to logging output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.

    verbosity : {'error', 'warn', 'info', 'hint'} or int, (default: 'info')
        0/'error', 1/'warn', 2/'info', 3/'hint' or integer 4.
    time : bool, optional (default: False)
        Print timing information; restart the clock.
    memory : bool, optional (default: Faulse)
        Print memory information.
    reset : bool, optional (default: False)
        Reset timing and memory measurement. Is automatically reset
        when passing one of ``time`` or ``memory``.
    end : str (default: '\n')
        Same meaning as in builtin ``print()`` function.
    no_indent : bool (default: False)
        Do not indent for ``verbosity >= 4``.
    """
    # all deprecated variable namings
    if v is not None: verbosity = v
    if t is not None: time = t
    if m is not None: memory = m
    if r is not None: reset = r
    if isinstance(verbosity, str):
        verbosity = verbosity_levels_from_strings[verbosity]
    if isinstance(settings.verbosity, str):
        global_verbosity = verbosity_levels_from_strings[settings.verbosity]
    else:
        global_verbosity = settings.verbosity
    if verbosity == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if verbosity >= 4 and not no_indent:
        msg = ('   ',) + msg
    if global_verbosity >= verbosity:
        if not time and not m and len(msg) > 0:
            settings.mi(*msg, end=end)
        if reset:
            settings._previous_memory_usage, _ = get_memory_usage()
            settings._previous_time = time_module.time()
        if time:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(settings._sec_to_str(elapsed)),)
            settings.mi(*msg, end=end)
        if memory:
            settings.mi(format_memory_usage(get_memory_usage()),
                        msg='' if time else msg, end=end)

m = msg


def get_passed_time():
    now = time_module.time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def get_memory_usage():
    import psutil
    process = psutil.Process(os.getpid())
    meminfo_attr = ('memory_info' if hasattr(process, 'memory_info')
                    else 'get_memory_info')
    mem = getattr(process, meminfo_attr)()[0] / 2**30  # output in GB
    mem_diff = mem
    if settings._previous_memory_usage != -1:
        mem_diff = mem - settings._previous_memory_usage
    settings._previous_memory_usage = mem
    return mem, mem_diff


def format_memory_usage(mem_usage, msg='', newline=False):
    mem, diff = mem_usage
    string = (('\n' if newline else '')
              + msg + (' \n... ' if msg != '' else '')
              + 'Memory usage: current {:.2f} GB, difference {:+.2f} GB'
              .format(mem, diff))
    return string


def print_memory_usage(msg='', newline=False):
    string = format_memory_usage(get_memory_usage(), msg, newline)
    settings.mi(string)


def print_version_and_date():
    from . import __version__
    settings.mi('Running Scanpy', __version__, 'on {}.'.format(get_date_string()))


_dependencies_numerics = [
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
    print('Dependencies:', end=' ')
    for mod in dependencies:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try:
            imp = __import__(mod_name)
            print('{}=={}'.format(mod_install, imp.__version__), end=' ')
        except (ImportError, AttributeError):
            pass
    print()

    
def print_versions_dependencies_numerics():
    """Dependencies that might influence numerical results (computed data).
    """
    _print_versions_dependencies(_dependencies_numerics)


def print_versions_dependencies_plotting():
    """Dependencies that might influence plots (but not computed data).
    """
    _print_versions_dependencies(_dependencies_plotting)


def print_versions_dependencies_all():
    """All relevant dependencies.
    """
    _print_versions_dependencies(
        _dependencies_numerics + _dependencies_plotting)


def get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

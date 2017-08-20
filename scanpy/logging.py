# Author: F. Alex Wolf (http://falexwolf.de)
"""Logging and Profiling
"""

import os
import time
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


def msg(*msg, verbosity='info', t=False, m=False, r=False, end='\n', no_indent=False,
        v=None):
    """Write message to logging output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.

    verbosity : {'error', 'warn', 'info', 'hint'} or int, (default: 'info')
        0/'error', 1/'warn', 2/'info', 3/'hint' or integer 4.
    t/timing : bool, optional (default: False)
        Print timing information.
    m/memory : bool, optional (default: Faulse)
        Print memory information.
    r/reset : bool, optional (default: False)
        Reset timing and memory measurement. Is automatically reset
        when passing one of t or m.
    end : str (default: '\n')
        As in python builtin print function.
    no_indent : bool (default: False)
        Indent depending on verbosity level.
    """
    if v is not None: verbosity = v
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
        if not t and not m and len(msg) > 0:
            settings.mi(*msg, end=end)
        if r:
            settings._previous_memory_usage, _ = get_memory_usage()
            settings._previous_time = time.time()
        if t:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(settings._sec_to_str(elapsed)),)
            settings.mi(*msg, end=end)
        if m:
            settings.mi(format_memory_usage(get_memory_usage()),
                    msg='' if t else msg, end=end)

m = msg

            
def get_passed_time():
    now = time.time()
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
    settings.mi('Running Scanpy version', __version__, 'on {}.'.format(get_date()))


def print_imported_modules():
    import types
    for _, var in sorted(globals().items()):
        if isinstance(var, types.ModuleType):
            print(var.__name__, var.__version__, end=' | ')


def get_date():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

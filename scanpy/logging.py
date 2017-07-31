# Author: F. Alex Wolf (http://falexwolf.de)
"""Logging and Profiling
"""

import os
import time
import datetime
from . import settings as sett


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
    if isinstance(sett.verbosity, str):
        global_v = verbosity_levels_from_strings[sett.verbosity]
    else:
        global_v = sett.verbosity
    return global_v >= v


def m(*msg, v='info', t=False, m=False, r=False, end='\n'):
    """Write message to log output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.sett.log_file = 'mylogfile.txt'`.

    v/verbosity : {'error', 'warn', 'info', 'hint'} or int, (default: 'info')
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
    """
    if isinstance(v, str):
        v = verbosity_levels_from_strings[v]
    if isinstance(sett.verbosity, str):
        global_v = verbosity_levels_from_strings[sett.verbosity]
    else:
        global_v = sett.verbosity
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v >= 4:
        msg = ('   ',) + msg
    if global_v >= v:
        if not t and not m and len(msg) > 0:
            sett.mi(*msg, end=end)
        if r:
            sett._previous_memory_usage, _ = get_memory_usage()
            sett._previous_time = time.time()
        if t:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(sett._sec_to_str(elapsed)),)
            sett.mi(*msg, end=end)
        if m:
            sett.mi(format_memory_usage(get_memory_usage()), msg='' if t else msg, end=end)


def get_passed_time():
    now = time.time()
    elapsed = now - sett._previous_time
    sett._previous_time = now
    return elapsed


def get_memory_usage():
    import psutil
    process = psutil.Process(os.getpid())
    meminfo_attr = ('memory_info' if hasattr(process, 'memory_info')
                    else 'get_memory_info')
    mem = getattr(process, meminfo_attr)()[0] / 2**30  # output in GB
    mem_diff = mem
    if sett._previous_memory_usage != -1:
        mem_diff = mem - sett._previous_memory_usage
    sett._previous_memory_usage = mem
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
    sett.mi(string)


def print_version_and_date():
    from . import __version__
    sett.mi('Running Scanpy version', __version__, 'on {}.'.format(get_date()))


def get_date():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

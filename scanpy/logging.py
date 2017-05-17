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


def m(*msg, v='info', t=False, m=False, r=False):
    """Write message to log output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.sett.log_file = 'mylogfile.txt'`.

    v/verbosity : {'error', 'warn', 'info', 'hint'} or int, (default: 'info')
         0/'error', 1/'warn', 2/'info', 3/'hint' or integer 4.
    t/timing : bool, optional (default: False)
         Prepend timing information.
    m/memory : bool, optional (default: Faulse)
         Append memory information.
    r/reset : bool, optional (default: False)
        Reset timing and memory measurement. Is automatically reset
        when passing one of t, m or p.
    """
    if isinstance(v, str):
        v = verbosity_levels_from_strings[v]
    if isinstance(sett.verbosity, str):
        global_v = verbosity_levels_from_strings[sett.verbosity]
    else:
        global_v = sett.verbosity
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v <= global_v:
        if not t and not m and not r:
            sett.mi(*msg)
        if r:
            sett._previous_memory_usage, _ = get_memory_usage()
            sett._previous_time = time.time()
        if t or r:
            elapsed = get_passed_time()
            sett.mi(sett._sec_to_str(elapsed), '-', *msg)
        if m:
            sett.mi(format_memory_usage(get_memory_usage()), msg='' if t else msg)


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
    print(string)


def get_date():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

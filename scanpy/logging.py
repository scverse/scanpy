"""Logging and Profiling
"""
<<<<<<< HEAD

import time as time_module
import datetime
from anndata import logging
from . import settings


_VERBOSITY_LEVELS_FROM_STRINGS = {
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


def _settings_verbosity_greater_or_equal_than(v):
    if isinstance(settings.verbosity, str):
        settings_v = _VERBOSITY_LEVELS_FROM_STRINGS[settings.verbosity]
    else:
        settings_v = settings.verbosity
    return settings_v >= v


def msg(*msg, v=4, time=False, memory=False, reset=False, end='\n',
        no_indent=False, t=None, m=None, r=None):
    """Write message to logging output.
=======
import logging
from functools import update_wrapper, partial
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET
from datetime import datetime, timedelta, timezone
from typing import Optional

import anndata.logging


HINT = (INFO + DEBUG) // 2
logging.addLevelName(HINT, 'HINT')


class _RootLogger(logging.RootLogger):
    def __init__(self, level):
        super().__init__(level)
        self.propagate = False
        _RootLogger.manager = logging.Manager(self)

    def log(
        self,
        level: int,
        msg: str,
        *,
        extra: Optional[dict] = None,
        time: datetime = None,
        deep: Optional[str] = None,
    ) -> datetime:
        from . import settings
        now = datetime.now(timezone.utc)
        time_passed: timedelta = None if time is None else now - time
        extra = {
            **(extra or {}),
            'deep': deep if settings.verbosity.level < level else None,
            'time_passed': time_passed
        }
        super().log(level, msg, extra=extra)
        return now

    def critical(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(CRITICAL, msg, time=time, deep=deep, extra=extra)

    def error(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(ERROR, msg, time=time, deep=deep, extra=extra)

    def warning(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(WARNING, msg, time=time, deep=deep, extra=extra)

    def info(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(INFO, msg, time=time, deep=deep, extra=extra)

    def hint(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(HINT, msg, time=time, deep=deep, extra=extra)

    def debug(self, msg, *, time=None, deep=None, extra=None) -> datetime:
        return self.log(DEBUG, msg, time=time, deep=deep, extra=extra)


def _set_log_file(settings):
    file = settings.logfile
    name = settings.logpath
    root = settings._root_logger
    h = logging.StreamHandler(file) if name is None else logging.FileHandler(name)
    h.setFormatter(_LogFormatter())
    h.setLevel(root.level)
    if len(root.handlers) == 1:
        root.removeHandler(root.handlers[0])
    elif len(root.handlers) > 1:
        raise RuntimeError('Scanpy’s root logger somehow got more than one handler')
    root.addHandler(h)


def _set_log_level(settings, level: int):
    root = settings._root_logger
    root.setLevel(level)
    h, = root.handlers  # may only be 1
    h.setLevel(level)


class _LogFormatter(logging.Formatter):
    def __init__(self, fmt='{levelname}: {message}', datefmt='%Y-%m-%d %H:%M', style='{'):
        super().__init__(fmt, datefmt, style)

    def format(self, record: logging.LogRecord):
        format_orig = self._style._fmt
        if record.levelno == INFO:
            self._style._fmt = '{message}'
        elif record.levelno == HINT:
            self._style._fmt = '--> {message}'
        elif record.levelno == DEBUG:
            self._style._fmt = '    {message}'
        if record.time_passed:
            # strip microseconds
            if record.time_passed.microseconds:
                record.time_passed = timedelta(seconds=int(record.time_passed.total_seconds()))
            if '{time_passed}' in record.msg:
                record.msg = record.msg.replace('{time_passed}', str(record.time_passed))
            else:
                self._style._fmt += ' ({time_passed})'
        if record.deep:
            record.msg = f'{record.msg}: {record.deep}'
        result = logging.Formatter.format(self, record)
        self._style._fmt = format_orig
        return result


print_memory_usage = anndata.logging.print_memory_usage
get_memory_usage = anndata.logging.get_memory_usage


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


def _versions_dependencies(dependencies):
    # this is not the same as the requirements!
    for mod in dependencies:
        mod_name, dist_name = mod if isinstance(mod, tuple) else (mod, mod)
        try:
            imp = __import__(mod_name)
            yield dist_name, imp.__version__
        except (ImportError, AttributeError):
            pass


def print_versions():
    """\
    Versions that might influence the numerical results.

    Matplotlib and Seaborn are excluded from this.
    """
    from ._settings import settings
    modules = ['scanpy'] + _DEPENDENCIES_NUMERICS
    print(' '.join(
        f'{mod}=={ver}'
        for mod, ver in _versions_dependencies(modules)
    ), file=settings.logfile)


def print_version_and_date():
    """\
    Useful for starting a notebook so you see when you started working.
    """
    from . import __version__
    from ._settings import settings
    print(
        f'Running Scanpy {__version__}, '
        f'on {datetime.now():%Y-%m-%d %H:%M}.',
        file=settings.logfile,
    )


def _copy_docs_and_signature(fn):
    return partial(update_wrapper, wrapped=fn, assigned=['__doc__', '__annotations__'])


def error(
    msg: str,
    *,
    time: datetime = None,
    deep: Optional[str] = None,
    extra: Optional[dict] = None,
) -> datetime:
    """\
    Log message with specific level and return current time.

    Parameters
    ==========
    msg
        Message to display.
    time
        A time in the past. If this is passed, the time difference from then
        to now is appended to `msg` as ` (HH:MM:SS)`.
        If `msg` contains `{time_passed}`, the time difference is instead
        inserted at that position.
    deep
        If the current verbosity is higher than the log function’s level,
        this gets displayed as well
    extra
        Additional values you can specify in `msg` like `{time_passed}`.
    """
    from ._settings import settings
    return settings._root_logger.error(msg, time=time, deep=deep, extra=extra)
>>>>>>> upstream/master


<<<<<<< HEAD
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
        v = _VERBOSITY_LEVELS_FROM_STRINGS[v]
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v >= 4 and not no_indent:
        msg = ('   ',) + msg
    if _settings_verbosity_greater_or_equal_than(v):
        if not time and not memory and len(msg) > 0:
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
        if memory:
            _write_log(format_memory_usage(get_memory_usage()),
                        msg='' if time else msg, end=end)

m = msg  # backwards compat


def _write_log(*msg, end='\n'):
    """Write message to log output, ignoring the verbosity level.

    This is the most basic function.

    Parameters
    ----------
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    from .settings import logfile
    if logfile == '':
        print(*msg, end=end)
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with open(logfile, 'a') as f:
            f.write(out + end)


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
    _write_log('Running Scanpy', __version__, 'on {}.'.format(get_date_string()))


_DEPENDENCIES_NUMERICS = [
    'anndata',  # anndata actually shouldn't, but as long as it's in development
    'numpy',
    'scipy',
    'pandas',
    ('sklearn', 'scikit-learn'),
    'statsmodels',
    ('igraph', 'python-igraph'),
    'louvain']


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
=======
@_copy_docs_and_signature(error)
def warning(msg, *, time=None, deep=None, extra=None) -> datetime:
    from ._settings import settings
    return settings._root_logger.warning(msg, time=time, deep=deep, extra=extra)
>>>>>>> upstream/master

def print_versions():
    """Versions that might influence the numerical results.

<<<<<<< HEAD
    Matplotlib and Seaborn are excluded from this.
    """
    _print_versions_dependencies(['scanpy'] + _DEPENDENCIES_NUMERICS)


def print_versions_dependencies_numerics():
    """Dependencies that might influence numerical results (computed data).
    """
    print('Dependencies:', end=' ')
    _print_versions_dependencies(_DEPENDENCIES_NUMERICS)


def print_versions_dependencies_plotting():
    """Dependencies that might influence plots (but not computed data).
    """
    print('Dependencies:', end=' ')
    _print_versions_dependencies(_DEPENDENCIES_PLOTTING)


def print_versions_dependencies_all():
    """All relevant dependencies.
    """
    _print_versions_dependencies(
        _DEPENDENCIES_NUMERICS + _DEPENDENCIES_PLOTTING)


def get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
=======
@_copy_docs_and_signature(error)
def info(msg, *, time=None, deep=None, extra=None) -> datetime:
    from ._settings import settings
    return settings._root_logger.info(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def hint(msg, *, time=None, deep=None, extra=None) -> datetime:
    from ._settings import settings
    return settings._root_logger.hint(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def debug(msg, *, time=None, deep=None, extra=None) -> datetime:
    from ._settings import settings
    return settings._root_logger.debug(msg, time=time, deep=deep, extra=extra)
>>>>>>> upstream/master

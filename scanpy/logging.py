"""Logging and Profiling
"""
import logging
import time as time_
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET
from datetime import datetime

import anndata.logging


HINT = (INFO + DEBUG) // 2
logging.addLevelName(HINT, 'HINT')


class RootLogger(logging.RootLogger):
    def __init__(self, level):
        super().__init__(level)
        self.propagate = False
        self.last_time = time_.time()

    def log(self, level, msg, *, extra=None, deep=None, time_passed=None):
        extra = {
            **(extra or {}),
            'deep': deep,
            'time_passed': time_passed,
        }
        super().log(level, msg, extra=extra)

    def error(self, msg, *, deep=None, extra=None):
        self.log(ERROR, msg, deep=deep, extra=extra)

    def warning(self, msg, *, deep=None, extra=None):
        self.log(WARNING, msg, deep=deep, extra=extra)

    def info(self, msg, *, time=False, deep=None, extra=None):
        current_time = time_.time()
        time_passed = None
        if time:
            time_passed = time_.strftime('%H:%M:%S', time_.gmtime(current_time - self.last_time))
        self.log(INFO, msg, time_passed=time_passed, deep=deep, extra=extra)
        self.last_time = time_.time()

    def hint(self, msg, *, deep=None, extra=None):
        self.log(HINT, msg, deep=deep, extra=extra)

    def debug(self, msg, *, deep=None, extra=None):
        self.log(DEBUG, msg, deep=deep, extra=extra)


def _set_log_file(settings):
    file = settings.logfile
    name = settings.logpath
    root = settings._root_logger
    h = logging.StreamHandler(file) if name is None else logging.FileHandler(name)
    h.setFormatter(LogFormatter())
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


class LogFormatter(logging.Formatter):
    def __init__(self, fmt='{levelname}: {message}', datefmt='%Y-%m-%d %H:%M', style='{'):
        super().__init__(fmt, datefmt, style)

    def format(self, record: logging.LogRecord):
        format_orig = self._style._fmt
        if record.levelno == INFO:
            self._style._fmt = '{asctime} | {message}'
            if '{time_passed}' not in record.msg and record.time_passed:
                self._style._fmt += ' ({time_passed})'
        elif record.levelno == HINT:
            self._style._fmt = '--> {message}'
        elif record.levelno == DEBUG:
            self._style._fmt = '    {message}'
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
    """Versions that might influence the numerical results.

    Matplotlib and Seaborn are excluded from this.
    """
    from ._settings import settings
    modules = ['scanpy'] + _DEPENDENCIES_NUMERICS
    print(' '.join(
        f'{mod}=={ver}'
        for mod, ver in _versions_dependencies(modules)
    ), file=settings.logfile)


def print_version_and_date():
    from . import __version__
    from ._settings import settings
    print(
        f'Running Scanpy {__version__}, '
        f'on {datetime.now():%Y-%m-%d %H:%M}.',
        file=settings.logfile,
    )


# will be replaced in settings
def error(msg, *, deep=None): pass
def warn(msg, *, deep=None): pass
def info(msg, *, deep=None, time=False): pass
def hint(msg, *, deep=None): pass
def debug(msg, *, deep=None): pass

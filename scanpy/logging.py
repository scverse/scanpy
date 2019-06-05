"""Logging and Profiling
"""
import logging
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET
from datetime import datetime, timedelta, timezone
from typing import Optional

import anndata.logging


HINT = (INFO + DEBUG) // 2
logging.addLevelName(HINT, 'HINT')


class RootLogger(logging.RootLogger):
    def __init__(self, level):
        super().__init__(level)
        self.propagate = False

    def log(
        self,
        level: int,
        msg: str,
        *,
        extra: Optional[dict] = None,
        time: datetime = None,
        deep: Optional[str] = None,
    ) -> datetime:
        now = datetime.now(timezone.utc)
        time_passed: timedelta = None if time is None else now - time
        extra = {
            **(extra or {}),
            'deep': deep,
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
        elif record.levelno == HINT:
            self._style._fmt = '--> {message}'
        elif record.levelno == DEBUG:
            self._style._fmt = '    {message}'
        if '{time_passed}' not in record.msg and record.time_passed:
            self._style._fmt += ' ({time_passed})'
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
def error(msg, *, time=None, deep=None, extra=None) -> datetime: return datetime.now()
def warning(msg, *, time=None, deep=None, extra=None) -> datetime: return datetime.now()
def info(msg, *, time=None, deep=None, extra=None) -> datetime: return datetime.now()
def hint(msg, *, time=None, deep=None, extra=None) -> datetime: return datetime.now()
def debug(msg, *, time=None, deep=None, extra=None) -> datetime: return datetime.now()

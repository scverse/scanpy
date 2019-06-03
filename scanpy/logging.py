"""Logging and Profiling
"""
import logging
import time
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET
from datetime import datetime

import anndata.logging


HINT = (INFO + DEBUG) / 2
logging.addLevelName(HINT, 'HINT')


class RootLogger(logging.RootLogger):
    def __init__(self, level):
        super().__init__(level)
        self.propagate = False

    def hint(self, msg, *args, **kwargs):
        return self.log(HINT, msg, *args, **kwargs)


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
    def __init__(self, fmt='%(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M', style='%', passed_time=False):
        super().__init__(fmt, datefmt, style)
        self.passed_time = passed_time
        self.last_time = time.time()

    def format(self, record: logging.LogRecord):
        format_orig = self._style._fmt
        if record.levelno == logging.INFO:
            current_time = time.time()
            passed_time_str = time.strftime('%H:%M:%S', time.gmtime(current_time - self.last_time))
            if self.passed_time:
                self._style._fmt = passed_time_str + ' - %(message)s'
            else:
                self._style._fmt = f'%(asctime)s | {passed_time_str} - %(message)s'
            self.last_time = time.time()
        if record.levelno == logging.DEBUG:
            self._style._fmt = '%(message)s'
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
def error(msg, *args, **kwargs): pass
def warn(msg, *args, **kwargs): pass
def info(msg, *args, **kwargs): pass
def hint(msg, *args, **kwargs): pass
def debug(msg, *args, **kwargs): pass

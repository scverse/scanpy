from packaging import version

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    cache = lru_cache(maxsize=None)

try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        pass


__all__ = ['cache', 'DaskArray', 'fullname', 'pkg_metadata', 'pkg_version']


def fullname(typ: type) -> str:
    module = typ.__module__
    name = typ.__qualname__
    if module == 'builtins' or module is None:
        return name
    return f'{module}.{name}'


def pkg_metadata(package):
    from importlib.metadata import metadata as m

    return m(package)


@cache
def pkg_version(package):
    from importlib.metadata import version as v

    return version.parse(v(package))

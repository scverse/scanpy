from packaging import version

try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        pass


__all__ = ['DaskArray', 'DaskTask', 'fullname', 'pkg_metadata', 'pkg_version']


def fullname(typ: type) -> str:
    module = typ.__module__
    name = typ.__qualname__
    if module == 'builtins' or module is None:
        return name
    return f'{module}.{name}'


def pkg_metadata(package):
    from importlib.metadata import metadata as m

    return m(package)


def pkg_version(package):
    from importlib.metadata import version as v

    return version.parse(v(package))

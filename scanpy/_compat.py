from packaging import version

try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type('Literal_', (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


def pkg_metadata(package):
    try:
        from importlib.metadata import metadata as m
    except ImportError:  # < Python 3.8: Use backport module
        from importlib_metadata import metadata as m
    return m(package)


def pkg_version(package):
    try:
        from importlib.metadata import version as v
    except ImportError:  # < Python 3.8: Use backport module
        from importlib_metadata import version as v
    return version.parse(v(package))

"""External exporting functions."""

from __future__ import annotations

from scverse_misc import Deprecation, deprecated


@deprecated(
    Deprecation("1.13.0", "Import `write_cellbrowser` from `scanpy.io` instead")
)
def cellbrowser(*args, **kwargs):
    from ..io import write_cellbrowser

    return write_cellbrowser(*args, **kwargs)


@deprecated(
    Deprecation("1.13.0", "Import `write_spring_project` from `scanpy.io` instead")
)
def spring_project(*args, **kwargs):
    from ..io import write_spring_project

    return write_spring_project(*args, **kwargs)

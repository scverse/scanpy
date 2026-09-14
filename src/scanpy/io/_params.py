"""Reading and writing parameter files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._utils.random import _LegacyRng

if TYPE_CHECKING:
    from os import PathLike


def read_params(
    filename: PathLike[str] | str, *, as_header: bool = False
) -> dict[str, int | float | bool | str | None]:
    """Read parameter dictionary from text file.

    Assumes that parameters are specified in the format::

        par1 = value1
        par2 = value2

    Comments that start with '#' are allowed.

    Parameters
    ----------
    filename
        Filename of data file.
    asheader
        Read the dictionary from the header (comment section) of a file.

    Returns
    -------
    Dictionary that stores parameters.

    """
    filename = Path(filename)  # allow passing str objects
    from collections import OrderedDict

    params = OrderedDict([])
    with filename.open() as f:
        for line_raw in f:
            if "=" not in line_raw or (as_header and not line_raw.startswith("#")):
                continue
            line = line_raw[1:] if line_raw.startswith("#") else line_raw
            key, val = line.split("=")
            key = key.strip()
            val = val.strip()
            params[key] = convert_string(val)
    return params


def write_params(path: PathLike[str] | str, *args, **maps):
    """Write parameters to file, so that it's readable by read_params.

    Uses INI file format.
    """
    path = Path(path)
    if not path.parent.is_dir():
        path.parent.mkdir(parents=True)
    if len(args) == 1:
        maps[None] = args[0]
    with path.open("w") as f:
        for header, map in maps.items():
            if header is not None:
                f.write(f"[{header}]\n")
            for key, val in map.items():
                if key == "rng":
                    if isinstance(val, _LegacyRng) and val.arg is not None:
                        f.write(f"seed = {val.arg}\n")
                else:
                    f.write(f"{key} = {val}\n")


def is_float(string: str) -> float:
    """Check whether string is float.

    See Also
    --------
    https://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python

    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string: str) -> bool:
    """Check whether string is integer."""
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string: str) -> tuple[bool, bool]:
    """Check whether string is boolean."""
    if string == "True":
        return True, True
    elif string == "False":
        return True, False
    else:
        return False, False


def convert_string(string: str) -> int | float | bool | str | None:
    """Convert string to int, float or bool."""
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == "None":
        return None
    else:
        return string

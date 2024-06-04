from __future__ import annotations

import sys
from enum import Enum, auto
from importlib.util import find_spec

import pytest


def _next_val(name: str, start: int, count: int, last_values: list[str]) -> str:
    """Distribution name for matching modules"""
    return name.replace("_", "-")


class QuietMarkDecorator(pytest.MarkDecorator):
    def __init__(self, mark: pytest.Mark) -> None:
        super().__init__(mark, _ispytest=True)


class needs(QuietMarkDecorator, Enum):
    """
    Pytest skip marker evaluated at module import.

    This allows us to see the amount of skipped tests at the start of a test run.
    :func:`pytest.importorskip` skips tests after they started running.
    """

    # _generate_next_value_ needs to come before members, also it’s finnicky:
    # https://github.com/python/mypy/issues/7591#issuecomment-652800625
    _generate_next_value_ = (
        staticmethod(_next_val) if sys.version_info >= (3, 10) else _next_val
    )

    dask = auto()
    dask_ml = auto()
    fa2 = auto()
    gprofiler = "gprofiler-official"
    leidenalg = auto()
    louvain = auto()
    openpyxl = auto()
    igraph = auto()
    pybiomart = auto()
    skimage = "scikit-image"
    skmisc = "scikit-misc"
    zarr = auto()
    zappy = auto()
    # external
    bbknn = auto()
    harmony = "harmonyTS"
    harmonypy = auto()
    magic = "magic-impute"
    palantir = auto()
    phate = auto()
    phenograph = auto()
    pypairs = auto()
    samalg = "sam-algorithm"
    scanorama = auto()
    trimap = auto()
    wishbone = "wishbone-dev"

    def __init__(self, mod: str) -> None:
        reason = f"needs module `{self._name_}`"
        if self._name_.casefold() != mod.casefold().replace("-", "_"):
            reason = f"{reason} (`pip install {mod}`)"
        dec = pytest.mark.skipif(not find_spec(self._name_), reason=reason)
        super().__init__(dec.mark)

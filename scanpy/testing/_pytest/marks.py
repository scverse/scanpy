from __future__ import annotations
from enum import Enum, auto
from importlib.util import find_spec
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pytest


class needs(pytest.MarkDecorator, Enum):
    """
    Pytest skip marker evaluated at module import.

    This allows us to see the amount of skipped tests at the start of a test run.
    :func:`pytest.importorskip` skips tests after they started running.
    """

    # _generate_next_value_ needs to come before members
    @staticmethod
    def _generate_next_value_(
        name: str, start: int, count: int, last_values: list[str]
    ) -> str:
        """Distribution name for matching modules"""
        return name.replace("_", "-")

    dask = auto()
    dask_ml = auto()
    fa2 = auto()
    gprofiler = "gprofiler-official"
    leidenalg = auto()
    louvain = auto()
    igraph = auto()
    pybiomart = auto()
    skimage = "scikit-image"
    skmisc = "scikit-misc"
    zarr = auto()
    zappy = auto()
    # external
    harmony = "harmonyTS"
    harmonypy = auto()
    magic = "magic-impute"
    palantir = auto()
    phate = auto()
    phenograph = "PhenoGraph"
    pypairs = auto()
    samalg = "sam-algorithm"
    scanorama = auto()
    scrublet = auto()
    trimap = auto()
    wishbone = "wishbone-dev"

    def __init__(self, mod: str) -> None:
        reason = f"needs module `{self._name_}`"
        if self._name_.casefold() != mod.casefold().replace("-", "_"):
            reason = f"{reason} (`pip install {mod}`)"
        dec = pytest.mark.skipif(not find_spec(self._name_), reason=reason)
        super().__init__(dec.mark)

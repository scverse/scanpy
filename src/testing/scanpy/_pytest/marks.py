from __future__ import annotations

from enum import Enum, auto
from importlib.util import find_spec
from typing import TYPE_CHECKING

import pytest
from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Callable


SKIP_EXTRA: dict[str, Callable[[], str | None]] = {}


def _skip_if_skmisc_too_old() -> str | None:
    import numpy as np
    import skmisc

    if Version(skmisc.__version__) <= Version("0.3.1") and Version(
        np.__version__
    ) >= Version("2"):
        return "scikit-misc≤0.3.1 requires numpy<2"
    return None


SKIP_EXTRA["skmisc"] = _skip_if_skmisc_too_old


class QuietMarkDecorator(pytest.MarkDecorator):
    def __init__(self, mark: pytest.Mark) -> None:
        super().__init__(mark, _ispytest=True)


class needs(QuietMarkDecorator, Enum):
    """Pytest skip marker evaluated at module import.

    This allows us to see the amount of skipped tests at the start of a test run.
    :func:`pytest.importorskip` skips tests after they started running.
    """

    # _generate_next_value_ needs to come before members
    @staticmethod
    def _generate_next_value_(
        name: str, start: int, count: int, last_values: list[str]
    ) -> str:
        """Distribution name for matching modules."""
        return name.replace("_", "-")

    mod: str

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
        self.mod = mod
        reason = self.skip_reason
        dec = pytest.mark.skipif(bool(reason), reason=reason or "")
        super().__init__(dec.mark)

    @property
    def skip_reason(self) -> str | None:
        if find_spec(self._name_):
            if skip_extra := SKIP_EXTRA.get(self._name_):
                return skip_extra()
            return None
        reason = f"needs module `{self._name_}`"
        if self._name_.casefold() != self.mod.casefold().replace("-", "_"):
            reason = f"{reason} (`pip install {self.mod}`)"
        return reason

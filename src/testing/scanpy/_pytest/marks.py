from __future__ import annotations

from enum import Enum, auto
from importlib.metadata import distributions, requires
from importlib.util import find_spec

import pytest
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name


def _missing_scanpy2_deps() -> list[Requirement]:
    dist_names = {canonicalize_name(d.name) for d in distributions()}
    return [
        r
        for r in map(Requirement, requires("scanpy") or ())
        if r.marker
        and r.marker.evaluate({"extra": "scanpy2"}, "requirement")
        and canonicalize_name(r.name) not in dist_names
    ]


class QuietMarkDecorator(pytest.MarkDecorator):
    def __init__(self, mark: pytest.Mark) -> None:
        super().__init__(mark, _ispytest=True)


class needs(QuietMarkDecorator, Enum):  # noqa: N801
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

    scanpy2 = "scanpy[scanpy2]"

    colour = "colour-science"
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
    illico = auto()
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
        if self._name_ == "scanpy2":
            if not (missing := _missing_scanpy2_deps()):
                return None
            return f"scanpy 2 deps missing: {', '.join(m.name for m in missing)}"

        if find_spec(self._name_):
            return None
        reason = f"needs module `{self._name_}`"
        if self._name_.casefold() != self.mod.casefold().replace("-", "_"):
            reason = f"{reason} (`pip install {self.mod}`)"
        return reason

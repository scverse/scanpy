from __future__ import annotations

from typing import TypedDict


class UmapMethodKwds(TypedDict, total=False):
    dens_lambda: float
    dens_frac: float
    dens_var_shift: float

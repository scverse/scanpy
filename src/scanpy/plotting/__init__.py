"""Plotting functions and classes."""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

from . import legacy
from ._common import dot_area

if TYPE_CHECKING:
    from types import ModuleType
    from typing import Any

__all__ = ["dot_area", "legacy"]


def _v2_module() -> ModuleType:
    return importlib.import_module("scanpy.plotting._v2")


def __dir__() -> list[str]:
    from scanpy._settings import Preset, settings

    if settings.preset is Preset.ScanpyV2Preview:
        return sorted(set(_v2_module().__all__) | set(__all__))
    return sorted(set(legacy.__all__) | set(__all__))


def __getattr__(name: str) -> Any:
    from scanpy._settings import Preset, settings

    _backend = _v2_module() if settings.preset is Preset.ScanpyV2Preview else legacy
    try:
        return getattr(_backend, name)
    except AttributeError:
        msg = f"module 'scanpy.plotting' has no attribute {name!r}"
        raise AttributeError(msg) from None

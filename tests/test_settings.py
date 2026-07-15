from __future__ import annotations

import inspect

import pytest

import scanpy as sc
from scanpy._settings import presets
from testing.scanpy._pytest import marks


# TODO: reset everything
@pytest.mark.parametrize("_attempt", range(3))
def test_resets(_attempt: int) -> None:
    """Test that changes made reset."""
    assert sc.settings.autoshow
    sc.settings.autoshow = False
    assert not sc.settings.autoshow


def test_set_figure_params_warns() -> None:
    with pytest.warns(FutureWarning, match=r"scanpy\.set_figure_params"):
        sc.settings.set_figure_params()


def test_preset_scanpy_v2_preview_checks_deps() -> None:
    if presets._missing_scanpy2_deps():
        with pytest.raises(ImportError, match=r"scanpy\[scanpy2\]"):
            sc.settings.preset = sc.Preset.ScanpyV2Preview
    else:
        sc.settings.preset = sc.Preset.ScanpyV2Preview
        assert sc.settings.preset is sc.Preset.ScanpyV2Preview
        sc.settings.preset = sc.Preset.ScanpyV1
    assert sc.settings.preset is sc.Preset.ScanpyV1


@pytest.mark.parametrize(
    "func", ["_missing_scanpy2_deps", "_req_satisfied", "dist_names"]
)
def test_no_divergence(func: str) -> None:
    """Unfortunately this function has to be duplicated.

    - we can’t import `scanpy` too early for coverage
    - we can’t import `testing.scanpy` in `scanpy`
    """
    a, b = (inspect.getsource(getattr(mod, func)) for mod in [presets, marks])
    assert a == b

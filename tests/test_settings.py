from __future__ import annotations

import pytest

import scanpy as sc
from scanpy._settings.presets import _missing_scanpy2_deps


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
    if _missing_scanpy2_deps():
        with pytest.raises(ImportError, match=r"scanpy\[scanpy2\]"):
            sc.settings.preset = sc.Preset.ScanpyV2Preview
    else:
        sc.settings.preset = sc.Preset.ScanpyV2Preview
        assert sc.settings.preset is sc.Preset.ScanpyV2Preview
        sc.settings.preset = sc.Preset.ScanpyV1
    assert sc.settings.preset == sc.Preset.ScanpyV1

from __future__ import annotations

from importlib.metadata import packages_distributions, requires

import pytest
from packaging.requirements import Requirement

import scanpy as sc


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
    dists = {d for m, ds in packages_distributions().items() for d in ds}
    scanpy2_deps_missing = any(
        r.name in dists
        for r in map(Requirement, requires("scanpy"))
        if r.marker and r.marker.evaluate({"extra": "scanpy2"})
    )

    if scanpy2_deps_missing:
        with pytest.raises(ImportError, match=r"scanpy\[scanpy2\]"):
            sc.settings.preset = sc.Preset.ScanpyV2Preview
    else:
        sc.settings.preset = sc.Preset.ScanpyV2Preview
        assert sc.settings.preset is sc.Preset.ScanpyV2Preview
        sc.settings.preset = sc.Preset.ScanpyV1
    assert sc.settings.preset == sc.Preset.ScanpyV1

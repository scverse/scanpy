from __future__ import annotations

import pytest

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

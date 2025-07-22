"""Some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from contextlib import chdir
from typing import TYPE_CHECKING

import numpy as np
import pytest

from .data import (
    backed_adata,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
    pbmc3ks_parametrized_session,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

__all__ = [
    "_doctest_env",
    "backed_adata",
    "float_dtype",
    "pbmc3k_parametrized",
    "pbmc3k_parametrized_small",
    "pbmc3ks_parametrized_session",
]


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param


@pytest.fixture
def _doctest_env(cache: pytest.Cache, tmp_path: Path) -> Generator[None, None, None]:
    showwarning_orig = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None) -> None:  # noqa: PLR0917
        if file is None:
            if line is None:
                import linecache

                line = linecache.getline(filename, lineno)
            line = line.strip()
            print(f"{category.__name__}: {message}\n    {line}")
        else:
            showwarning_orig(message, category, filename, lineno, file, line)

    # ignore plt.show() warning only in doctests.
    warnings.filterwarnings(
        "ignore",
        r"FigureCanvasAgg is non-interactive, and thus cannot be shown",
        UserWarning,
    )
    # make errors visible and the rest ignored
    action_map = defaultdict(lambda: "ignore", error="default")
    warnings.filters = [
        (action_map[action], *rest) for action, *rest in warnings.filters
    ]

    warnings.showwarning = showwarning
    with chdir(tmp_path):
        yield
    warnings.showwarning = showwarning_orig

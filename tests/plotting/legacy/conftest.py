from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from packaging.version import Version

from scanpy._compat import pkg_version

if TYPE_CHECKING:
    import os


HERE = Path(__file__).parent


@pytest.fixture(scope="session")
def image_dir() -> Path:
    return HERE / "_images"


@pytest.fixture
def plot_cmp(check_same_image, image_dir: Path):
    from matplotlib import pyplot as plt

    def save_and_compare(
        path_str: Path | os.PathLike, *, root: Path | None = None, tol: int = 15
    ):
        __tracebackhide__ = True
        if root is None:
            root = image_dir
        if tol > 25:
            msg = "with this loose of a tolerance, you might as well delete the test"
            raise ValueError(msg)

        mpl_dir = "3.11" if pkg_version("matplotlib") >= Version("3.11") else "3.10"
        base_pth = root / mpl_dir / path_str

        if not base_pth.is_dir():
            base_pth.mkdir()
        expected_pth = base_pth / "expected.png"
        actual_pth = base_pth / "actual.png"
        plt.savefig(actual_pth, dpi=40)
        plt.close()
        if not expected_pth.is_file():
            msg = f"No expected output found at {expected_pth}."
            raise OSError(msg)
        check_same_image(expected_pth, actual_pth, tol=tol, root=root)

    return save_and_compare

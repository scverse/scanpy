from __future__ import annotations

from pathlib import PurePosixPath, PureWindowsPath

import pytest

from scanpy.readwrite import _slugify


@pytest.mark.parametrize(
    "path",
    [
        PureWindowsPath(r"C:\foo\bar"),
        PureWindowsPath(r".\C\foo\bar"),
        PureWindowsPath(r"C\foo\bar"),
        PurePosixPath("/C/foo/bar"),
        PurePosixPath("./C/foo/bar"),
        PurePosixPath("C/foo/bar"),
    ],
)
def test_slugify(path):
    assert _slugify(path) == "C-foo-bar"

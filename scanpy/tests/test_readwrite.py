from pathlib import PureWindowsPath, PurePosixPath

import pytest

from scanpy.readwrite import _slugify


@pytest.mark.parametrize('path', [
    PureWindowsPath(r'C:\foo\bar'),
    PureWindowsPath(r'.\C\foo\bar'),
    PureWindowsPath(r'C\foo\bar'),
    PurePosixPath('/C/foo/bar'),
    PurePosixPath('./C/foo/bar'),
    PurePosixPath('C/foo/bar'),
])
def test_slugify(path):
    assert _slugify(path) == 'C-foo-bar'

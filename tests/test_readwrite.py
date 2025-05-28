from __future__ import annotations

from pathlib import PurePosixPath, PureWindowsPath

import numpy as np
import pytest
from anndata import AnnData

import scanpy as sc
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


def test_read_ext_match(tmp_path):
    adata_path = tmp_path / "foo.bar.anndata.h5ad"
    AnnData(np.array([[1, 2], [3, 4]])).write_h5ad(adata_path)
    with pytest.raises(ValueError, match="does not end in expected extension"):
        sc.read(adata_path, ext="zarr")
    # should not warn: https://github.com/scverse/scanpy/issues/2288
    sc.read(adata_path, ext="h5ad")

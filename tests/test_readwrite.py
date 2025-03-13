from __future__ import annotations

from contextlib import nullcontext
from pathlib import PurePosixPath, PureWindowsPath
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pytest
import zarr
from anndata import AnnData
from anndata.tests.helpers import assert_equal

import scanpy as sc
from scanpy.readwrite import _slugify

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


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


@pytest.mark.parametrize("ext", ["h5ad", "zarr", "csv"])
@pytest.mark.parametrize("style", ["path", "ext", "default"])
def test_write(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    ext: Literal["h5ad", "zarr", "csv"],
    style: Literal["path", "ext", "default"],
) -> None:
    monkeypatch.chdir(tmp_path)
    adata = AnnData(np.array([[1, 2], [3, 4]]))

    # test that writing works (except style="default" and ext="csv")
    ctx = (
        pytest.warns(FutureWarning, match=r"removed from this function")
        if ext == "csv"
        else nullcontext()
    )
    match style, ext:
        case "path", _:
            with ctx:
                sc.write(f"test.{ext}", adata)
            d = tmp_path
        case "ext", _:
            with ctx:
                sc.write("test", adata, ext=ext)
            d = sc.settings.writedir
        case "default", "csv":
            # check that it throws an error instead
            ff = sc.settings.file_format_data
            with pytest.raises(ValueError, match=r"Cannot set file_format_data to csv"):
                sc.settings.file_format_data = ext  # type: ignore[assignment]
            assert sc.settings.file_format_data == ff
            return  # return early
        case "default", _:
            monkeypatch.setattr(sc.settings, "file_format_data", ext)
            with ctx:
                sc.write("test", adata)
            d = sc.settings.writedir
        case _:
            pytest.fail("add branch for new style")

    path = d / ("test" if ext == "csv" else f"test.{ext}")
    assert tuple(d.iterdir()) == (path,)
    assert path.is_file() if ext == "h5ad" else path.is_dir()

    # test that roundtripping works
    if ext != "csv":  # no reader for this
        adata_read = sc.read(path)
        assert_equal(adata_read, adata)


@pytest.mark.parametrize("fmt", ["h5ad", "zarr"])
@pytest.mark.parametrize("s2c", [True, False], ids=["s2c", "no_s2c"])
def test_write_strings_to_cats(fmt: Literal["h5ad", "zarr"], *, s2c: bool) -> None:
    adata = AnnData(np.array([[1, 2], [3, 4]]), obs=dict(a=["a", "b"]))

    sc.write("test", adata.copy(), convert_strings_to_categoricals=s2c, ext=fmt)
    p = sc.settings.writedir / f"test.{fmt}"

    if fmt == "h5ad":
        with h5py.File(p, "r") as f:
            et = f["obs"]["a"].attrs["encoding-type"]
    elif fmt == "zarr":
        et = zarr.open_group(p, "r")["obs"]["a"].attrs["encoding-type"]
    else:
        pytest.fail("add branch for new fmt")
    assert et == ("categorical" if s2c else "string-array")

    adata_read = sc.read(p)
    assert_equal(adata_read, adata)
    assert adata_read.obs["a"].dtype == "object"

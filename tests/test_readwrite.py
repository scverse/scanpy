from __future__ import annotations

from contextlib import nullcontext, suppress
from pathlib import PurePosixPath, PureWindowsPath
from typing import TYPE_CHECKING

import anndata.io
import numpy as np
import pytest
from anndata import AnnData
from anndata.tests.helpers import assert_equal

import scanpy as sc
from scanpy.io._download import slugify
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable
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
    assert slugify(path) == "C-foo-bar"


def test_read_ext_match(tmp_path):
    adata_path = tmp_path / "foo.bar.anndata.h5ad"
    AnnData(np.array([[1, 2], [3, 4]])).write_h5ad(adata_path)
    with pytest.raises(ValueError, match="does not end in expected extension"):
        sc.io.read(adata_path, ext="zarr")
    # should not warn: https://github.com/scverse/scanpy/issues/2288
    sc.io.read(adata_path, ext="h5ad")


@pytest.mark.parametrize("ext", ["h5ad", pytest.param("zarr", marks=needs.zarr), "csv"])
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
                sc.io.write(f"test.{ext}", adata)
            d = tmp_path
        case "ext", _:
            with ctx:
                sc.io.write("test", adata, ext=ext)
            d = sc.settings.writedir
        case "default", "csv":
            # check that it throws an error instead
            ff = sc.settings.file_format_data
            with pytest.raises(ValueError, match=r"should be 'h5ad' or 'zarr'.*'csv'"):
                sc.settings.file_format_data = ext  # type: ignore[assignment]
            assert sc.settings.file_format_data == ff
            return  # return early
        case "default", _:
            sc.settings.file_format_data, old = ext, sc.settings.file_format_data
            try:
                with ctx:
                    sc.io.write("test", adata)
            finally:
                sc.settings.file_format_data = old
            d = sc.settings.writedir
        case _:
            pytest.fail("add branch for new style")

    path = d / ("test" if ext == "csv" else f"test.{ext}")
    assert tuple(d.iterdir()) == (path,)
    assert path.is_file() if ext == "h5ad" else path.is_dir()

    # test that roundtripping works
    if ext != "csv":  # no reader for this
        adata_read = sc.io.read(path)
        assert_equal(adata_read, adata)


@pytest.mark.parametrize("fmt", ["h5ad", pytest.param("zarr", marks=needs.zarr)])
@pytest.mark.parametrize("s2c", [True, False], ids=["s2c", "no_s2c"])
def test_write_strings_to_cats(fmt: Literal["h5ad", "zarr"], *, s2c: bool) -> None:
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]), obs=dict(a=["a", "b", "a"]))

    sc.io.write("test", adata, convert_strings_to_categoricals=s2c, ext=fmt)
    p = sc.settings.writedir / f"test.{fmt}"
    adata_read = sc.io.read(p)

    assert_equal(adata_read, adata)
    assert adata_read.obs["a"].dtype == adata.obs["a"].dtype
    assert adata_read.obs["a"].dtype in (("category",) if s2c else ("object", "string"))


@pytest.mark.parametrize(
    "access",
    [
        pytest.param(lambda: sc.read, id="read"),
        pytest.param(lambda: sc.read_10x_h5, id="read_10x_h5"),
        pytest.param(lambda: sc.read_10x_mtx, id="read_10x_mtx"),
        pytest.param(lambda: sc.read_csv, id="read_csv"),
        pytest.param(lambda: sc.read_excel, id="read_excel"),
        pytest.param(lambda: sc.read_h5ad, id="read_h5ad"),
        pytest.param(lambda: sc.read_hdf, id="read_hdf"),
        pytest.param(lambda: sc.read_mtx, id="read_mtx"),
        pytest.param(lambda: sc.read_text, id="read_text"),
        pytest.param(lambda: sc.read_umi_tools, id="read_umi_tools"),
        pytest.param(lambda: sc.write, id="write"),
        pytest.param(sc.external.exporting.cellbrowser, id="exporting.cellbrowser"),
        pytest.param(
            sc.external.exporting.spring_project, id="exporting.spring_project"
        ),
    ],
)
def test_moved_to_io(access: Callable[[], object]) -> None:
    with (
        pytest.warns(FutureWarning, match=r"from `scanpy\.io`"),
        suppress(TypeError),  # the `sc.external.exporting` shims get no arguments
    ):
        access()


def test_read_loom_deprecated() -> None:
    with pytest.warns(FutureWarning, match=r"`read_loom` is deprecated"):
        assert sc.read_loom is anndata.io.read_loom

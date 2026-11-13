from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scanpy.external as sce
from testing.scanpy._pytest.marks import needs


@pytest.fixture
def cb_stub(monkeypatch):
    cb = Mock()
    package = ModuleType("cellbrowser")
    package.cellbrowser = cb
    monkeypatch.setitem(sys.modules, "cellbrowser", package)
    monkeypatch.setitem(sys.modules, "cellbrowser.cellbrowser", cb)
    return cb


@pytest.mark.parametrize(
    "path_kind",
    ["relative-str", "relative-path", "absolute-str", "absolute-path", "parent"],
)
def test_cellbrowser_paths(cb_stub, tmp_path, path_kind):
    data_dir = tmp_path / "données with spaces"
    if path_kind == "parent":
        data_dir = data_dir / ".." / "project"
    if not path_kind.startswith("absolute"):
        # Keep '..' components: resolving them can change symlink semantics.
        data_dir = Path(os.path.relpath(tmp_path)) / data_dir.relative_to(tmp_path)
    expected = str(data_dir if data_dir.is_absolute() else Path.cwd() / data_dir)
    if path_kind.endswith("str"):
        data_dir = str(data_dir)
    adata = object()
    coords = {"X_umap": "UMAP"}
    annotations = {"sample": "Sample"}
    with pytest.warns(FutureWarning, match="cellbrowser"):
        result = sce.exporting.cellbrowser(
            adata,
            data_dir,
            "tiny",
            embedding_keys=coords,
            annot_keys=annotations,
            cluster_field="sample",
            nb_marker=7,
            skip_matrix=True,
            html_dir=Path("relative html"),
            port=8765,
            do_debug=True,
        )
    assert result is None
    cb_stub.setDebug.assert_called_once_with(True)  # noqa: FBT003
    cb_stub.scanpyToCellbrowser.assert_called_once_with(
        adata,
        expected,
        "tiny",
        coordFields=coords,
        metaFields=annotations,
        clusterField="sample",
        nb_marker=7,
        skipMatrix=True,
        doDebug=None,
    )
    cb_stub.build.assert_called_once_with(expected, "relative html", doDebug=None)
    cb_stub.serve.assert_called_once_with("relative html", 8765)


def test_cellbrowser_symlink(cb_stub, tmp_path):
    target = tmp_path / "target" / "nested"
    target.mkdir(parents=True)
    link = tmp_path / "link"
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError as exc:
        pytest.skip(f"Symlinks unavailable: {exc}")
    data_dir = Path(os.path.relpath(link)) / ".." / "project"
    expected = str(Path.cwd() / data_dir)
    assert expected != str(data_dir.resolve())
    with pytest.warns(FutureWarning, match="cellbrowser"):
        sce.exporting.cellbrowser(object(), data_dir, "tiny")
    assert cb_stub.scanpyToCellbrowser.call_args.args[1] == expected


@pytest.mark.parametrize(
    ("html_dir", "port"), [(None, None), (None, 8765), ("relative html", None)]
)
def test_cellbrowser_defaults(cb_stub, tmp_path, html_dir, port):
    adata = object()
    with pytest.warns(FutureWarning, match="cellbrowser"):
        sce.exporting.cellbrowser(adata, tmp_path, "tiny", html_dir=html_dir, port=port)
    cb_stub.setDebug.assert_called_once_with(False)  # noqa: FBT003
    cb_stub.scanpyToCellbrowser.assert_called_once_with(
        adata,
        str(tmp_path),
        "tiny",
        coordFields=None,
        metaFields=("louvain", "percent_mito", "n_genes", "n_counts"),
        clusterField="louvain",
        nb_marker=50,
        skipMatrix=False,
        doDebug=None,
    )
    if html_dir is None:
        cb_stub.build.assert_not_called()
    else:
        cb_stub.build.assert_called_once_with(str(tmp_path), html_dir, doDebug=None)
    cb_stub.serve.assert_not_called()


@pytest.fixture
def adata() -> AnnData:
    adata = AnnData(
        np.array([[1, 0, 3], [2, 1, 0], [0, 4, 1], [3, 2, 1]], dtype=np.float32),
        obs=pd.DataFrame(
            {
                "louvain": pd.Categorical(["a", "a", "b", "b"]),
                "sample": pd.Categorical(["one", "two", "one", "two"]),
            },
            index=["cell1", "cell2", "cell3", "cell4"],
        ),
        var=pd.DataFrame(index=["gene1", "gene2", "gene3"]),
        obsm={"X_umap": np.array([[0, 1], [1, 0], [2, 3], [3, 2]])},
        uns={
            "sample_colors": ["#ff0000", "#0000ff"],
            "rank_genes_groups": {
                "names": np.array([("gene1", "gene2")], dtype=[("a", "O"), ("b", "O")]),
                "scores": np.array([(2.0, 3.0)], dtype=[("a", "f4"), ("b", "f4")]),
            },
        },
    )
    return adata


@pytest.mark.filterwarnings(
    "ignore:The 'warn' function is deprecated:DeprecationWarning"
)
@pytest.mark.filterwarnings(
    "ignore:The inplace keyword in DataFrame.rename is deprecated"
)
@pytest.mark.filterwarnings("ignore:unclosed file:ResourceWarning")
@pytest.mark.parametrize("relative", [True, False], ids=["relative", "absolute"])
@needs.cellbrowser
def test_cellbrowser_export_build(adata, tmp_path, monkeypatch, relative):
    cb = pytest.importorskip("cellbrowser.cellbrowser")  # noqa: TID251

    def nonspatial(adata, out_dir, config, coords):
        # CellBrowser's spatial helper imports removed private Scanpy APIs,
        # even when it would return immediately for this non-spatial fixture.
        assert "spatial" not in adata.obsm
        return config, coords

    monkeypatch.setattr(cb, "exportScanpySpatial", nonspatial)

    def forbidden(*args, **kwargs):
        pytest.fail("CellBrowser must not download data or start a server in this test")

    monkeypatch.setattr(cb, "downloadUrlBinary", forbidden)
    monkeypatch.setattr(cb, "serve", forbidden)
    monkeypatch.setattr(cb, "startHttpServer", forbidden)
    # Optional reference marker databases require downloads; JS/CSS still copy normally.
    monkeypatch.setattr(cb, "copyMarkers", Mock())
    # Use CellBrowser's Python gzip/md5 fallbacks without external processes.
    monkeypatch.setattr(cb, "which", lambda command: None)
    monkeypatch.setattr(cb, "setDebug", Mock())
    monkeypatch.setattr(cb, "getConfig", lambda tag, defValue=None: defValue)  # noqa: N803

    data_dir = tmp_path / "données"
    html_dir = tmp_path / "html"
    data_arg = os.path.relpath(data_dir) if relative else data_dir
    html_arg = os.path.relpath(html_dir)
    kwargs = dict(
        embedding_keys=["X_umap"],
        annot_keys=["louvain", "sample"],
        cluster_field=None,  # Use the supplied markers, avoiding automatic DE.
    )
    with pytest.warns(FutureWarning, match="cellbrowser"):
        sce.exporting.cellbrowser(adata, data_arg, "tiny", html_dir=html_arg, **kwargs)

    conf = cb.loadConfig(str(data_dir / "cellbrowser.conf"))
    assert Path(conf["colors"]["sample"]).is_absolute()
    assert Path(conf["colors"]["sample"]).samefile(data_dir / "sample_colors.tsv")
    colors = pd.read_csv(data_dir / "sample_colors.tsv", sep="\t")
    assert colors.to_numpy().tolist() == [["one", "#ff0000"], ["two", "#0000ff"]]
    coords = pd.read_csv(data_dir / "umap_coords.tsv", sep="\t", index_col=0)
    assert coords.index.tolist() == adata.obs_names.tolist()
    np.testing.assert_array_equal(coords.to_numpy(), adata.obsm["X_umap"])

    dataset_dir = html_dir / "tiny"
    dataset = json.loads((dataset_dir / "dataset.json").read_text())
    assert dataset["sampleCount"] == adata.n_obs
    sample = next(field for field in dataset["metaFields"] if field["name"] == "sample")
    assert dict(
        zip((row[0] for row in sample["valCounts"]), sample["colors"], strict=True)
    ) == {
        "one": "ff0000",
        "two": "0000ff",
    }
    assert (html_dir / "index.html").is_file()
    for directory in (data_dir, dataset_dir):
        matrix = pd.read_csv(directory / "exprMatrix.tsv.gz", sep="\t", index_col=0)
        assert matrix.columns.tolist() == adata.obs_names.tolist()
        assert matrix.index.tolist() == adata.var_names.tolist()
        np.testing.assert_array_equal(matrix.to_numpy(), adata.X.T)
        meta = pd.read_csv(directory / "meta.tsv", sep="\t", index_col=0)
        assert meta.index.tolist() == adata.obs_names.tolist()
        assert meta["sample"].tolist() == adata.obs["sample"].tolist()

    original_conf = (data_dir / "cellbrowser.conf").read_bytes()
    original_matrix = (data_dir / "exprMatrix.tsv.gz").read_bytes()
    adata.X = adata.X + 10
    with pytest.warns(FutureWarning, match="cellbrowser"):
        sce.exporting.cellbrowser(adata, data_arg, "tiny", skip_matrix=True, **kwargs)
    assert (data_dir / "cellbrowser.conf").read_bytes() == original_conf
    assert (data_dir / "exprMatrix.tsv.gz").read_bytes() == original_matrix

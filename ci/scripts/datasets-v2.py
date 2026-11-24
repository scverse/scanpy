#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "anndata==0.13.4",  # needs the read-time compat code removed in scverse/anndata#2673
#   "scanpy==1.13.0",
#   "igraph==1.0.0",
#   "numpy==2.5.5",
#   "scipy==1.18.1",
#   "pandas==3.1.1",
#   "h5py==3.16.0",
#   "zarr==3.4.2",
# ]
# [tool.uv]
# exclude-newer = "2026-11-16T00:00:00Z"  # pins the rest of the dependency tree
# ///
"""Convert v1 example datasets into v2s that don’t rely on anndata read-time compat.

For all datasets:
- zipped zarr v3 instead of h5ad
- v2 embedding keys (`scanpy.datasets` restores v1’s for the v1 preset): `obsm` keys without
  `X_` prefix, `varm["PCs"]` → `varm["pca"]`, `draw_graph_fr`/`uns["draw_graph"]` → `graph_fr`
- `uns["neighbors"]` names its graph keys
- new `obs["leiden"]` clustering on the existing neighbors graph

`pbmc3k-processed` (from S3 key `scanpy/pbmc3k.h5ad`):
- current on-disk format: categorical `obs` columns, adjacency graphs in `obsp`
- `obs["louvain"]` → `obs["cell_type"]` (it holds cell type annotations, not clusters)
- `.raw` (log1p of the counts, all genes) → `X`, with `var["highly_variable"]` marking v1’s genes;
  v1’s scaled `X` is dropped, `varm["PCs"]` is zero outside the mask like `pp.pca(mask=…)`
- `uns["rank_genes_groups"]` recomputed (v1’s was stale: groups named "0"…"7")

`pbmc68k-reduced` (from `git show 43072780^:src/scanpy/datasets/10x_pbmc68k_reduced.h5ad`):
- `.raw` (log-normalized, rounded to 3 decimals) → `layers["counts"]`

Run on Linux x86_64 for the canonical sha256s: other platforms differ in the last digits of
p-values. From the repo root (paths must be inside it):

    docker run --rm --platform=linux/amd64 -v .:/work -w /work \
      ghcr.io/astral-sh/uv:python3.13-bookworm-slim@sha256:531f855bda2c73cd6ef67d56b733b357cea384185b3022bd09f05e002cd144ca \
      uv run --script ci/scripts/datasets-v2.py pbmc3k-processed data/pbmc3k_processed.h5ad pbmc3k_processed.zarr.zip
"""

from __future__ import annotations

import argparse
import hashlib
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZipInfo

import anndata as ad
import numpy as np
import pandas as pd
import zarr.storage

import scanpy as sc
from scanpy._compat import CSRBase


def pbmc3k_processed(adata: ad.AnnData) -> ad.AnnData:
    """Apply changes specific to `pbmc3k_processed`."""
    assert isinstance(adata.obs, pd.DataFrame)
    assert isinstance(adata.obs["louvain"].dtype, pd.CategoricalDtype)
    adata.obs = adata.obs.rename(columns={"louvain": "cell_type"})
    adata.uns["cell_type_colors"] = adata.uns.pop("louvain_colors")
    del adata.uns["louvain"]  # clustering params don’t describe a manual annotation

    # `.raw` (log1p of the counts, all genes) becomes `X`; the scaled HVG `X` is dropped.
    assert adata.raw is not None
    idx = adata.raw.var_names.get_indexer(adata.var_names)
    assert (idx >= 0).all()
    hvg = np.zeros(adata.raw.n_vars, dtype=bool)
    hvg[idx] = True
    # like `pp.pca(mask=…)`: zero loadings for genes outside the mask
    pcs = np.zeros(
        (adata.raw.n_vars, adata.varm["PCs"].shape[1]), adata.varm["PCs"].dtype
    )
    pcs[idx] = adata.varm["PCs"]
    adata = ad.AnnData(
        adata.raw.X,
        obs=adata.obs,
        var=adata.raw.var.assign(highly_variable=hvg),
        uns=adata.uns,
        obsm=dict(adata.obsm),
        varm=dict(PCs=pcs),
        obsp=dict(adata.obsp),
    )
    adata.uns["pca"]["params"] = dict(zero_center=True, mask_var="highly_variable")

    # same params as v1 (which used `.raw`), see `uns["rank_genes_groups"]["params"]` there
    del adata.uns["rank_genes_groups"]
    sc.tl.rank_genes_groups(
        adata,
        "cell_type",
        method="t-test_overestim_var",
        reference="rest",
        n_genes=100,
    )
    return adata


def pbmc68k_reduced(adata: ad.AnnData) -> ad.AnnData:
    """Apply changes specific to `pbmc68k_reduced`."""
    assert adata.raw is not None
    assert isinstance(adata.raw.X, CSRBase)
    assert list(adata.raw.var_names) == list(adata.var_names)
    # `.raw` is `log1p(normalize_total(counts, target_sum=1e4))` rounded to 3 decimals,
    # which is precise enough to restore the integer counts exactly.
    counts = adata.raw.X.copy()
    size_factors = adata.obs["n_counts"].to_numpy() / 1e4
    counts.data = np.expm1(counts.data) * np.repeat(
        size_factors, np.diff(counts.indptr)
    )
    assert np.abs(counts.data - np.rint(counts.data)).max() < 0.25
    adata.layers["counts"] = counts.rint().astype(np.int32)
    del adata.raw
    return adata


DATASETS = {"pbmc3k-processed": pbmc3k_processed, "pbmc68k-reduced": pbmc68k_reduced}


def convert(name: str, src: Path, dst: Path) -> None:
    """Read v1 from `src`, write v2 to `dst`."""
    with warnings.catch_warnings():
        # the compat code we rely on here warns about moving graphs to `obsp`
        warnings.simplefilter("ignore", (FutureWarning, ad.OldFormatWarning))
        adata = ad.read_h5ad(src)
    assert {"distances", "connectivities"} <= adata.obsp.keys()

    adata = DATASETS[name](adata)
    assert adata.raw is None

    # v2 key defaults, see `scanpy._keys._embedding_keys`
    for k in list(adata.obsm):
        adata.obsm[k.removeprefix("X_")] = adata.obsm.pop(k)
    adata.varm["pca"] = adata.varm.pop("PCs")
    if "draw_graph_fr" in adata.obsm:
        adata.obsm["graph_fr"] = adata.obsm.pop("draw_graph_fr")
        adata.uns["graph_fr"] = adata.uns.pop("draw_graph")

    adata.uns["neighbors"] |= dict(
        connectivities_key="connectivities", distances_key="distances"
    )

    # all params explicit so changing defaults (scanpy 2.0) don’t change the result
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False, rng=0)

    # writing into a `ZipStore` directly leaves duplicate `zarr.json` entries
    with TemporaryDirectory() as tmp:
        store = Path(tmp) / "adata.zarr"
        adata.write_zarr(store)
        zip_dir(store, dst)


def zip_dir(src: Path, dst: Path) -> None:
    """Zip the contents of `src` with fixed order and timestamps, so the sha256 is reproducible."""
    with ZipFile(dst, "w") as zf:  # chunks are already compressed by zarr
        for p in sorted(src.rglob("*")):
            if p.is_file():
                info = ZipInfo(p.relative_to(src).as_posix(), (2026, 11, 16, 0, 0, 0))
                zf.writestr(info, p.read_bytes())


def main() -> None:
    """Convert, check the result, and print the sha256 for `registry.yaml`."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dataset", choices=DATASETS)
    parser.add_argument("src", type=Path, help="v1 `.h5ad` file")
    parser.add_argument("dst", type=Path, help="output path for v2 (`….zarr.zip`)")
    args = parser.parse_args()

    # write strings as `string-array`: older anndata refuses to write back `nullable-string-array`
    pd.set_option("future.infer_string", False)  # noqa: FBT003
    convert(args.dataset, args.src, args.dst)
    # read back the same way `scanpy.datasets` reads `10x_pbmc68k_reduced.zarr.zip`
    with zarr.storage.ZipStore(args.dst, mode="r") as store:
        v2 = ad.read_zarr(zarr.open_group(store=store, mode="r"))
    assert isinstance(v2.obs["leiden"].dtype, pd.CategoricalDtype)
    assert not any(k.startswith("X_") for k in v2.obsm)
    print(v2)
    print("sha256:", hashlib.sha256(args.dst.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()

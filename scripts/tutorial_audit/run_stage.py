"""Run one stage of the published clustering notebook and save its evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import pooch

import scanpy as sc

mpl.use("Agg")
from matplotlib import pyplot as plt

STAGES = {
    1: [1, 2, 4, 5],
    2: [9, 10, 12, 14, 16],
    3: [18],
    4: [22, 23],
    5: [25, 26],
    6: [28, 30, 32],
    7: [34, 36, 38],
    8: [41, 42, 45, 46, 52, 54, 55],
    9: [59, 60, 62, 63],
    10: [67, 69, 72, 73],
}


def main():
    """Execute original cells for one stage and record their output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", type=int, choices=STAGES)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    notebook = root / "docs/tutorials/basics/clustering.ipynb"
    content = notebook.read_bytes()
    cells = json.loads(content)["cells"]
    for directory in ["checkpoints", "plots", "metrics"]:
        (args.output / directory).mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 8
    sc.settings.autoshow = False
    sc.settings.set_figure_params(dpi=50, facecolor="white")
    namespace = {"ad": ad, "sc": sc, "pooch": pooch}
    if args.stage > 1:
        namespace["adata"] = ad.read_h5ad(
            args.output / "checkpoints" / f"{args.stage - 1:02}.h5ad"
        )
    timings = {}
    for index in STAGES[args.stage]:
        source = "".join(cells[index]["source"])
        print(f"START cell {index}\n{source}", flush=True)
        started = time.perf_counter()
        exec(compile(source, f"{notebook}:cell-{index}", "exec"), namespace)
        timings[str(index)] = time.perf_counter() - started
        for figure in plt.get_fignums():
            plt.figure(figure).savefig(
                args.output / "plots" / f"{args.stage:02}-cell-{index}-{figure}.png",
                bbox_inches="tight",
            )
        plt.close("all")
        print(f"DONE cell {index}: {timings[str(index)]:.3f}s", flush=True)
    adata = namespace["adata"]
    adata.write_h5ad(args.output / "checkpoints" / f"{args.stage:02}.h5ad")
    evidence = {
        "stage": args.stage,
        "cells": STAGES[args.stage],
        "cell_seconds": timings,
        "shape": list(adata.shape),
        "notebook_sha256": hashlib.sha256(content).hexdigest(),
        "scanpy_version": sc.__version__,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
    }
    (args.output / "metrics" / f"{args.stage:02}-execution.json").write_text(
        json.dumps(evidence, indent=2) + "\n"
    )
    print(json.dumps(evidence, indent=2), flush=True)


if __name__ == "__main__":
    main()

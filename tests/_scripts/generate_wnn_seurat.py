"""Regenerate the small, synthetic Seurat WNN regression fixture.

Requires NumPy and R with Seurat 5.5.1 and RANN installed::

    python tests/_scripts/generate_wnn_seurat.py --rscript /path/to/Rscript

Inputs are rounded to float32 before R reads them, matching the test input.
The accompanying R script runs Seurat itself with exact RANN neighbor search;
it does not use the Python WNN implementation or its NumPy test reference.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    script_dir = Path(__file__).resolve().parent

    rng = np.random.default_rng(0)
    n = 200
    labels = rng.integers(0, 4, n)
    reps = {}
    for name, dims, groups in (("pca", 12, 4), ("apca", 8, 2)):
        centers = rng.normal(size=(4, dims)) * 3
        reps[name] = (centers[labels % groups] + rng.normal(size=(n, dims))).astype(
            np.float32
        )

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        for name, values in reps.items():
            # Seventeen digits preserve the exact float32 values as R doubles.
            np.savetxt(directory / f"{name}.csv", values, delimiter=",", fmt="%.17g")
        subprocess.run(
            [args.rscript, str(script_dir / "generate_wnn_seurat.R"), tmp],
            check=True,
        )
        outputs = {
            name: np.loadtxt(directory / f"{name}.csv", delimiter=",")
            for name in ("weights", "indices", "distances", "snn")
        }
        outputs["indices"] = outputs["indices"].astype(np.int32)
        provenance = (
            f"NumPy {np.__version__}\n" + (directory / "provenance.txt").read_text()
        )

    destination = script_dir.parent / "_data" / "wnn_seurat_5.5.1.npz"
    np.savez_compressed(destination, **reps, **outputs, provenance=provenance)
    print(f"Wrote {destination} ({destination.stat().st_size:,} bytes)")
    print(provenance)


if __name__ == "__main__":
    main()

"""Independently audit full tutorial HVG metrics with NumPy, SciPy, and pandas."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import scanpy as sc

COLUMNS = ["means", "dispersions", "dispersions_norm"]


def per_batch(x, names):
    """Compute sample variance and binned log dispersion without Scanpy helpers."""
    x = x.astype(np.float64).tocsr()
    detected = np.asarray((x > 0).sum(axis=0)).ravel() > 0
    x = x[:, detected]
    n = x.shape[0]
    mean = np.asarray(x.sum(axis=0)).ravel() / n
    # A centered sum of squares avoids Scanpy's E[X²] - E[X]² formula.
    csc = x.tocsc()
    repeated_mean = np.repeat(mean, np.diff(csc.indptr))
    squared = csc.copy()
    squared.data = (csc.data - repeated_mean) ** 2
    ss = np.asarray(squared.sum(axis=0)).ravel()
    ss += (n - np.diff(csc.indptr)) * mean**2
    variance = ss / (n - 1)
    frame = pd.DataFrame(
        {"means": np.log1p(mean), "dispersions": np.log(variance / mean)},
        index=names[detected],
    )
    bins, edges = pd.cut(frame.means, bins=20, retbins=True)
    grouped = frame.dispersions.groupby(bins, observed=True)
    avg = grouped.transform("mean")
    std = grouped.transform("std")
    singleton = std.isna()
    std.loc[singleton] = avg.loc[singleton]
    avg.loc[singleton] = 0
    frame["dispersions_norm"] = (frame.dispersions - avg) / std
    z = frame.dispersions_norm.to_numpy()
    cutoff = np.sort(z[~np.isnan(z)])[-2000]
    frame["highly_variable"] = z >= cutoff
    frame["rank_min"] = frame.dispersions_norm.rank(ascending=False, method="min")
    info = {
        "cells": n,
        "expressed_genes": int(detected.sum()),
        "absent_genes": int((~detected).sum()),
        "selected": int(frame.highly_variable.sum()),
        "cutoff": float(cutoff),
        "cutoff_ties": int((z == cutoff).sum()),
        "nonfinite_normalized_dispersions": int((~np.isfinite(z)).sum()),
        "singleton_genes": frame.index[singleton].tolist(),
        "bin_counts": {str(k): int(v) for k, v in grouped.size().items()},
        "bin_edges": edges.tolist(),
    }
    frame = frame.reindex(names, fill_value=0)
    frame.loc[names[~detected], "rank_min"] = np.nan
    frame["highly_variable"] = frame.highly_variable.astype(bool)
    return frame, info


def combine(frames):
    """Combine equal-weight batch metrics and select exactly 2,000 genes."""
    joined = pd.concat(frames)
    result = joined.groupby(level=0, sort=True).agg(
        means=("means", "mean"),
        dispersions=("dispersions", "mean"),
        dispersions_norm=("dispersions_norm", "mean"),
        highly_variable_nbatches=("highly_variable", "sum"),
    )
    result = result.sort_values(
        ["highly_variable_nbatches", "dispersions_norm"], ascending=False
    )
    result["highly_variable"] = np.arange(len(result)) < 2000
    return result


def differences(actual, expected):
    """Measure the largest absolute error in each reported metric."""
    return {
        column: float(np.max(np.abs(actual[column] - expected[column])))
        for column in COLUMNS
    }


def main():  # noqa: PLR0915
    """Run the numerical audit and write its evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = ad.read_h5ad(args.output / "checkpoints/04.h5ad")
    after = ad.read_h5ad(args.output / "checkpoints/05.h5ad")
    assert before.shape == after.shape
    assert before.obs_names.equals(after.obs_names)
    assert before.var_names.equals(after.var_names)
    assert (before.X != after.X).nnz == 0
    assert (before.layers["counts"] != after.layers["counts"]).nnz == 0
    pd.testing.assert_frame_equal(before.obs, after.obs)
    pd.testing.assert_frame_equal(before.var, after.var[before.var.columns])
    assert after.uns["hvg"]["flavor"] == "seurat"
    assert before.uns["log1p"]["base"] is None
    counts = before.layers["counts"].astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    target = np.median(totals[totals > 0])
    normalized = counts.multiply((target / totals)[:, None]).tocsr()
    recovered = before.X.copy()
    recovered.data = np.expm1(recovered.data)
    delta = recovered.astype(np.float64) - normalized
    batches = before.obs["sample"].cat.categories
    frames, ideal_frames, details = [], [], {}
    for batch in batches:
        rows = (before.obs["sample"] == batch).to_numpy()
        frame, info = per_batch(recovered[rows], before.var_names)
        ideal, _ = per_batch(normalized[rows], before.var_names)
        observed = sc.pp.highly_variable_genes(
            before[rows].copy(),
            n_top_genes=2000,
            flavor="seurat",
            filter_unexpressed_genes=True,
            inplace=False,
        )
        for column in COLUMNS:
            np.testing.assert_allclose(
                observed[column], frame[column], rtol=1e-5, atol=1e-6
            )
        np.testing.assert_array_equal(observed.highly_variable, frame.highly_variable)
        info["scanpy_metric_max_absolute_errors"] = differences(observed, frame)
        info["scanpy_selection_disagreements"] = 0
        info["singleton_final_selected"] = {
            gene: bool(after.var.loc[gene, "highly_variable"])
            for gene in info["singleton_genes"]
        }
        info["count_float64_selection_disagreements"] = int(
            (frame.highly_variable != ideal.highly_variable).sum()
        )
        info["count_float64_metric_max_absolute_errors"] = differences(frame, ideal)
        frames.append(frame)
        ideal_frames.append(ideal)
        details[str(batch)] = info
    ranked = combine(frames)
    expected = ranked.loc[before.var_names]
    ideal = combine(ideal_frames).loc[before.var_names]
    actual = after.var
    for column in COLUMNS:
        np.testing.assert_allclose(
            actual[column], expected[column], rtol=1e-5, atol=1e-6
        )
    np.testing.assert_array_equal(actual.highly_variable, expected.highly_variable)
    np.testing.assert_array_equal(
        actual.highly_variable_nbatches, expected.highly_variable_nbatches
    )
    np.testing.assert_array_equal(
        actual.highly_variable_intersection,
        expected.highly_variable_nbatches == len(batches),
    )
    selected = actual.index[actual.highly_variable]
    assert len(selected) == 2000
    boundary = ranked.iloc[1999:2001]
    cutoff = ranked.iloc[1999]
    tied = ranked[
        (ranked.highly_variable_nbatches == cutoff.highly_variable_nbatches)
        & (ranked.dispersions_norm == cutoff.dispersions_norm)
    ]
    pooled, _ = per_batch(recovered, before.var_names)
    pure_dispersion = expected.sort_values("dispersions_norm", ascending=False).index[
        :2000
    ]
    weighted_means = (
        sum(
            frame.means * info["cells"]
            for frame, info in zip(frames, details.values(), strict=True)
        )
        / before.n_obs
    )
    selected_frame = expected.loc[selected]
    records = []
    for gene, row in boundary.iterrows():
        records.append({
            "gene": gene,
            "nbatches": int(row.highly_variable_nbatches),
            "dispersion_norm": float(row.dispersions_norm),
            "selected": bool(row.highly_variable),
        })
    external = args.output / "metrics"
    ranked.to_csv(external / "05-independent-gene-metrics.csv")
    for batch, frame in zip(batches, frames, strict=True):
        frame.to_csv(external / f"05-{batch}-gene-metrics.csv")
    gene_bytes = ("\n".join(selected) + "\n").encode()
    (external / "05-selected-genes.txt").write_bytes(gene_bytes)
    detected_cells = np.asarray((counts > 0).sum(axis=0)).ravel()
    selected_mask = actual.highly_variable.to_numpy()
    evidence = {
        "shape": list(before.shape),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["scanpy", "numpy", "scipy", "pandas", "fast-array-utils"]
        },
        "input_dtype": str(before.X.dtype),
        "stored_metric_dtypes": {
            column: str(actual[column].dtype) for column in COLUMNS
        },
        "normalization_target": float(target),
        "count_normalization_vs_expm1_max_absolute_error": float(abs(delta).max()),
        "input_matrix_counts_annotations_unchanged": True,
        "batches": details,
        "metric_max_absolute_errors": differences(actual, expected),
        "metric_tolerance": {"rtol": 1e-5, "atol": 1e-6},
        "selected_genes": len(selected),
        "selection_agreement": 1.0,
        "nbatches_and_intersection_exact": True,
        "selected_gene_list_sha256": hashlib.sha256(gene_bytes).hexdigest(),
        "selected_nbatches_counts": {
            str(k): int(v)
            for k, v in selected_frame.highly_variable_nbatches.value_counts().items()
        },
        "all_nbatches_counts": {
            str(k): int(v)
            for k, v in expected.highly_variable_nbatches.value_counts().items()
        },
        "final_boundary": records,
        "final_boundary_exact_ties": len(tied),
        "selected_low_mean_count": int((selected_frame.means <= 0.0125).sum()),
        "selected_high_mean_count": int((selected_frame.means >= 3).sum()),
        "selected_detected_in_at_most_10_cells": int(
            (detected_cells[selected_mask] <= 10).sum()
        ),
        "selected_detection_count_range": [
            int(detected_cells[selected_mask].min()),
            int(detected_cells[selected_mask].max()),
        ],
        "selected_absent_in_one_batch": int(
            sum(((frame.means == 0) & actual.highly_variable).sum() for frame in frames)
        ),
        "selected_qc_gene_groups": {
            column: int(actual.loc[selected, column].sum())
            for column in ["mt", "ribo", "hb"]
        },
        "selected_ignored_mean_cutoff_count": int(
            ((selected_frame.means <= 0.0125) | (selected_frame.means >= 3)).sum()
        ),
        "selected_ignored_dispersion_cutoff_count": int(
            (selected_frame.dispersions_norm <= 0.5).sum()
        ),
        "equal_batch_vs_cell_weighted_log_mean_max_difference": float(
            abs(weighted_means - expected.means).max()
        ),
        "count_float64_final_selection_disagreements": int(
            (ideal.highly_variable != actual.highly_variable).sum()
        ),
        "count_float64_metric_max_absolute_errors": differences(actual, ideal),
        "pooled_hvg_overlap": len(
            selected.intersection(pooled.index[pooled.highly_variable])
        ),
        "mean_dispersion_only_top2000_overlap": len(
            selected.intersection(pure_dispersion)
        ),
        "external_evidence": [
            str(external / "05-independent-gene-metrics.csv"),
            str(external / "05-selected-genes.txt"),
        ],
    }
    dest = Path(__file__).resolve().parents[2] / "audit/tutorial/evidence/05.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()

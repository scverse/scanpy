"""Audit tutorial normalization against independent float64 arithmetic."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

import scanpy as sc


def summary(values):
    """Return distribution summaries."""
    return dict(
        zip(
            ["min", "p05", "median", "p95", "max"],
            np.quantile(values, [0, 0.05, 0.5, 0.95, 1]).tolist(),
            strict=True,
        )
    )


def association(x, y):
    """Describe associations without causal claims."""
    return {
        "pearson": float(stats.pearsonr(x, y).statistic),
        "spearman": float(stats.spearmanr(x, y).statistic),
    }


def error_summary(actual, expected):
    """Compare positive values with float64 and rounded float32 references."""
    delta = actual.astype(np.float64) - expected
    ulps = abs(
        actual.view(np.int32).astype(np.int64)
        - expected.astype(np.float32).view(np.int32).astype(np.int64)
    )
    return {
        "max_absolute": float(np.max(abs(delta))),
        "max_relative": float(np.max(abs(delta) / expected)),
        "max_float32_ulps": int(ulps.max()),
        "sum_squared_error": float(delta @ delta),
    }


def main():  # noqa: PLR0915
    """Check all counts and report numerical and scientific effects."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = ad.read_h5ad(args.output / "checkpoints/03.h5ad")
    after = ad.read_h5ad(args.output / "checkpoints/04.h5ad")
    pd.testing.assert_frame_equal(before.obs, after.obs)
    pd.testing.assert_frame_equal(before.var, after.var)
    counts, logged = before.X, after.X
    saved = after.layers["counts"]
    for other in [saved, logged]:
        np.testing.assert_array_equal(counts.indptr, other.indptr)
        np.testing.assert_array_equal(counts.indices, other.indices)
    np.testing.assert_array_equal(counts.data, saved.data)
    assert counts.dtype == saved.dtype == logged.dtype == np.float32
    assert np.all(counts.data > 0)
    np.testing.assert_array_equal(counts.data, np.floor(counts.data))
    totals = np.asarray(counts.astype(np.float64).sum(axis=1)).ravel()
    target = float(np.median(totals[totals > 0]))
    assert np.all(totals > 0)
    assert after.uns["log1p"] == {"base": None}
    normalized = sc.pp.normalize_total(before, inplace=False)
    x = normalized["X"]
    np.testing.assert_array_equal(
        normalized["norm_factor"], totals.astype(np.float32) / np.float32(target)
    )
    prelog_totals = np.asarray(x.astype(np.float64).sum(axis=1)).ravel()
    log_totals = np.asarray(logged.astype(np.float64).sum(axis=1)).ravel()
    detected = np.diff(counts.indptr)
    samples = before.obs["sample"].astype(str).to_numpy()
    batches = sorted(set(samples))
    sample_targets = {s: float(np.median(totals[samples == s])) for s in batches}
    row_targets = np.array([sample_targets[s] for s in samples])
    errors = {kind: [] for kind in ["normalized", "log1p"]}
    highly_expressed = np.zeros(before.n_vars, dtype=bool)
    alternative = {
        s: {"n": 0, "sum": 0.0, "sum_squared": 0.0, "max": 0.0} for s in batches
    }
    for start in range(0, before.n_obs, 512):
        end = min(start + 512, before.n_obs)
        lo, hi = counts.indptr[[start, end]]
        rows = np.repeat(np.arange(start, end), detected[start:end])
        raw = counts.data[lo:hi].astype(np.float64)
        expected = raw * target / totals[rows]
        log_expected = np.log1p(expected)
        for kind, actual, reference in [
            ("normalized", x.data[lo:hi], expected),
            ("log1p", logged.data[lo:hi], log_expected),
        ]:
            errors[kind].append(error_summary(actual, reference))
            np.testing.assert_allclose(actual, reference, rtol=3e-7, atol=0)
        highly_expressed[counts.indices[lo:hi][raw > 0.05 * totals[rows]]] = True
        delta = np.log1p(raw * row_targets[rows] / totals[rows]) - log_expected
        for sample in batches:
            d = delta[samples[rows] == sample]
            result = alternative[sample]
            result["n"] += len(d)
            result["sum"] += float(d.sum())
            result["sum_squared"] += float(d @ d)
            result["max"] = max(result["max"], float(np.max(abs(d), initial=0)))
    for kind, chunks in errors.items():
        errors[kind] = {
            key: max(c[key] for c in chunks)
            for key in ["max_absolute", "max_relative", "max_float32_ulps"]
        } | {
            "rmse": float(
                np.sqrt(sum(c["sum_squared_error"] for c in chunks) / counts.nnz)
            )
        }
    sc.pp.log1p(x)
    np.testing.assert_array_equal(x.data, logged.data)
    np.testing.assert_array_equal(counts.data, saved.data)
    excluded_totals = np.asarray(
        counts[:, highly_expressed].astype(np.float64).sum(axis=1)
    ).ravel()
    denominator_without_high = totals - excluded_totals
    optional_target = float(
        np.median(denominator_without_high[denominator_without_high > 0])
    )
    optional = sc.pp.normalize_total(
        before, inplace=False, exclude_highly_expressed=True
    )
    np.testing.assert_array_equal(
        optional["norm_factor"],
        denominator_without_high.astype(np.float32) / np.float32(optional_target),
    )
    optional_max_relative_error = 0.0
    for start in range(0, before.n_obs, 512):
        end = min(start + 512, before.n_obs)
        lo, hi = counts.indptr[[start, end]]
        rows = np.repeat(np.arange(start, end), detected[start:end])
        expected = (
            counts.data[lo:hi].astype(np.float64)
            * optional_target
            / denominator_without_high[rows]
        )
        actual = optional["X"].data[lo:hi]
        np.testing.assert_allclose(actual, expected, rtol=3e-7, atol=0)
        optional_max_relative_error = max(
            optional_max_relative_error,
            float(np.max(abs(actual - expected) / expected)),
        )
    remaining_sums = np.asarray(
        optional["X"][:, ~highly_expressed].astype(np.float64).sum(axis=1)
    ).ravel()
    by_sample = {}
    for sample in batches:
        mask = samples == sample
        alt = alternative[sample]
        by_sample[sample] = {
            "cells": int(mask.sum()),
            "library_size": summary(totals[mask]),
            "scaling_multiplier": summary(target / totals[mask]),
            "one_count_log_value": summary(np.log1p(target / totals[mask])),
            "detected_genes": summary(detected[mask]),
            "sum_log1p": summary(log_totals[mask]),
            "library_vs_sum_log1p": association(totals[mask], log_totals[mask]),
            "library_vs_detected_genes": association(totals[mask], detected[mask]),
            "sample_specific_target": sample_targets[sample],
            "sample_specific_over_common_scale": sample_targets[sample] / target,
            "sample_specific_minus_common_log_nonzero_mean": alt["sum"] / alt["n"],
            "sample_specific_minus_common_log_nonzero_rms": float(
                np.sqrt(alt["sum_squared"] / alt["n"])
            ),
            "sample_specific_minus_common_log_max_absolute": alt["max"],
            "highly_expressed_gene_count_fraction": summary(
                excluded_totals[mask] / totals[mask]
            ),
        }
    evidence = {
        "baseline": "a6f1a2d2",
        "notebook_cells": [22, 23],
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["scanpy", "anndata", "numpy", "scipy"]
        },
        "shape": list(before.shape),
        "nonzero_values_checked": int(counts.nnz),
        "implicit_zeros_preserved": int(before.n_obs * before.n_vars - counts.nnz),
        "zero_library_cells": int(np.count_nonzero(totals == 0)),
        "dtype": str(counts.dtype),
        "counts_layer_exact": True,
        "obs_var_unchanged": True,
        "full_replay_bitwise_equal": True,
        "default_target": target,
        "library_size": summary(totals),
        "stored_qc_total_minus_current_total": summary(
            before.obs["total_counts"].to_numpy() - totals
        ),
        "prelog_row_totals": summary(prelog_totals),
        "prelog_max_absolute_target_error": float(np.max(abs(prelog_totals - target))),
        "prelog_max_relative_target_error": float(
            np.max(abs(prelog_totals / target - 1))
        ),
        "errors": errors,
        "pooled_library_vs_sum_log1p": association(totals, log_totals),
        "pooled_library_vs_detected_genes": association(totals, detected),
        "highly_expressed_exclusion_enabled": False,
        "optional_exclusion": {
            "threshold": 0.05,
            "all_nonzero_values_replayed": True,
            "max_relative_error": optional_max_relative_error,
            "remaining_prelog_sum_max_absolute_target_error": float(
                np.max(abs(remaining_sums - optional_target))
            ),
            "gene_count": int(highly_expressed.sum()),
            "genes": before.var_names[highly_expressed].tolist(),
            "remaining_denominator": summary(denominator_without_high),
            "remaining_zero_denominator_cells": int(
                np.count_nonzero(denominator_without_high == 0)
            ),
            "remaining_nonzero_median_target": float(
                np.median(denominator_without_high[denominator_without_high > 0])
            ),
        },
        "by_sample": by_sample,
    }
    destination = Path("audit/tutorial/evidence/04.json")
    destination.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()

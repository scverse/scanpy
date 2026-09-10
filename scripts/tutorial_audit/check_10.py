"""Audit every tutorial DGE result against independent NumPy/SciPy calculations."""

from __future__ import annotations

import argparse
import json
import subprocess
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy import stats

mpl.use("Agg")
from matplotlib import pyplot as plt

import scanpy as sc


def bh(p):
    """Adjust the complete gene family with the BH step-up formula."""
    order = np.argsort(p)
    q = np.empty_like(p)
    q[order] = np.minimum(
        1,
        np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[
            ::-1
        ],
    )
    return q


def aligned(result, field, group, genes):
    """Restore the matrix gene order for an independently checked field."""
    return (
        pd
        .Series(result[field][group], index=result["names"][group])
        .reindex(genes)
        .to_numpy()
    )


def compare(actual, expected, *, rtol=2e-6, atol=2e-6):
    """Record numerical errors without hiding small p-value differences."""
    return {
        "max_abs_error": float(np.max(np.abs(actual - expected))),
        "mismatch_count": int(
            (~np.isclose(actual, expected, rtol=rtol, atol=atol)).sum()
        ),
        "rtol": rtol,
        "atol": atol,
    }


def main():  # noqa: PLR0915
    """Run the complete dataset audit and write reviewable evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("/home/fdr/scanpy-tutorial-audit")
    )
    parser.add_argument("--chunk-size", type=int, default=128)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    evidence_dir = root / "audit/tutorial/evidence"
    sc.settings.n_jobs = 8
    sc.settings.autoshow = False
    a = ad.read_h5ad(args.output / "checkpoints/10.h5ad")
    result = a.uns["rank_genes_groups"]
    params = result["params"]
    x = a.raw.X if params["use_raw"] else a.X
    genes = a.raw.var_names if params["use_raw"] else a.var_names
    groups = list(result["names"].dtype.names)
    labels = a.obs["leiden_res_0.50"].astype(str).to_numpy()
    masks = [labels == group for group in groups]
    n, m = x.shape
    k = len(groups)
    z = np.empty((k, m))
    lfc = np.empty_like(z)
    arithmetic = np.empty_like(z)
    fractions = np.empty_like(z)
    means = np.empty_like(z)
    tc = np.empty(m)
    for left in range(0, m, args.chunk_size):
        right = min(m, left + args.chunk_size)
        values = x[:, left:right].toarray().astype(np.float64)
        ranks = stats.rankdata(values, axis=0, method="average")
        for j in range(right - left):
            _, counts = np.unique(values[:, j], return_counts=True)
            counts = counts.astype(np.float64)
            tc[left + j] = 1 - np.sum(counts**3 - counts) / (n**3 - n)
        linear = np.expm1(values)
        for g, mask in enumerate(masks):
            ng = int(mask.sum())
            z[g, left:right] = (ranks[mask].sum(axis=0) - ng * (n + 1) / 2) / np.sqrt(
                ng * (n - ng) * (n + 1) / 12
            )
            mean_g, mean_r = values[mask].mean(axis=0), values[~mask].mean(axis=0)
            lfc[g, left:right] = np.log2(
                (np.expm1(mean_g) + 1e-9) / (np.expm1(mean_r) + 1e-9)
            )
            means[g, left:right] = linear[mask].mean(axis=0)
            arithmetic[g, left:right] = np.log2(
                (means[g, left:right] + 1e-9) / (linear[~mask].mean(axis=0) + 1e-9)
            )
            fractions[g, left:right] = (values[mask] > 0).mean(axis=0)
        print(f"Independent genes {right}/{m}", flush=True)
    p = 2 * stats.norm.sf(np.abs(z))
    q = np.array([bh(row) for row in p])
    # Constant genes have zero variance in the tie-adjusted statistic.
    z_tie = np.divide(z, np.sqrt(tc), out=np.zeros_like(z), where=tc > 0)
    p_tie = 2 * stats.norm.sf(np.abs(z_tie))
    q_tie = np.array([bh(row) for row in p_tie])
    print("Run full tie-correction sensitivity", flush=True)
    sc.tl.rank_genes_groups(
        a,
        groupby="leiden_res_0.50",
        method="wilcoxon",
        tie_correct=True,
        key_added="audit_tie",
    )
    print("Run full arithmetic-fold-change sensitivity", flush=True)
    sc.tl.rank_genes_groups(
        a,
        groupby="leiden_res_0.50",
        method="wilcoxon",
        mean_in_log_space=False,
        key_added="audit_arithmetic",
    )
    print("Run top-five BH check", flush=True)
    sc.tl.rank_genes_groups(
        a,
        groupby="leiden_res_0.50",
        method="wilcoxon",
        n_genes=5,
        key_added="audit_top5",
    )
    checks = {}
    sensitivity = []
    marker_rows = []
    reversal_examples = []
    marker_names = [
        "NKG7",
        "KLRF1",
        "TRDC",
        "SH2D1B",
        "GNLY",
        "KLRD1",
        "PRF1",
        "FCGR3A",
        "CD3D",
        "CD3E",
        "TRAC",
        "CD8A",
        "CD8B",
        "MS4A1",
        "CD79A",
        "CD14",
        "LYZ",
        "FCER1A",
        "CLEC10A",
    ]
    for g, group in enumerate(groups):
        fields = {
            "scores": z[g],
            "pvals": p[g],
            "pvals_adj": q[g],
            "logfoldchanges": lfc[g],
        }
        checks[group] = {
            field: compare(
                aligned(result, field, group, genes),
                expected,
                atol=0 if field.startswith("pvals") else 2e-6,
            )
            for field, expected in fields.items()
        }
        checks[group]["names_unique_complete"] = bool(
            len(set(result["names"][group])) == m
            and set(result["names"][group]) == set(genes)
        )
        checks[group]["scores_descending"] = bool(
            np.all(np.diff(result["scores"][group]) <= 0)
        )
        df = sc.get.rank_genes_groups_df(a, group=group)
        checks[group]["get_alignment"] = bool(
            all(
                np.array_equal(df[field], result[field][group])
                for field in ["names", *fields]
            )
        )
        checks[group]["tie_scores"] = compare(
            aligned(a.uns["audit_tie"], "scores", group, genes), z_tie[g]
        )
        checks[group]["tie_pvals"] = compare(
            aligned(a.uns["audit_tie"], "pvals", group, genes), p_tie[g], atol=0
        )
        checks[group]["tie_bh"] = compare(
            aligned(a.uns["audit_tie"], "pvals_adj", group, genes), q_tie[g], atol=0
        )
        checks[group]["arithmetic_lfc"] = compare(
            aligned(a.uns["audit_arithmetic"], "logfoldchanges", group, genes),
            arithmetic[g],
        )
        checks[group]["arithmetic_scores_unchanged"] = bool(
            np.array_equal(
                aligned(a.uns["audit_arithmetic"], "scores", group, genes),
                aligned(result, "scores", group, genes),
            )
        )
        top = a.uns["audit_top5"]
        indices = genes.get_indexer(top["names"][group])
        checks[group]["top5_bh_full_family"] = compare(
            top["pvals_adj"][group], q[g, indices], atol=0
        )
        top_default = set(result["names"][group][:20])
        top_tie = set(a.uns["audit_tie"]["names"][group][:20])
        delta = arithmetic[g] - lfc[g]
        reversal = arithmetic[g] * lfc[g] < 0
        material = reversal & (abs(arithmetic[g]) > 1) & (abs(lfc[g]) > 1)
        material_significant = material & (q[g] < 0.05)
        # Include each material, significant pair with both expression fractions.
        reversal_examples.extend(
            {
                "group": group,
                "gene": genes[j],
                "q_default": float(q[g, j]),
                "log2fc_default": float(lfc[g, j]),
                "log2fc_arithmetic": float(arithmetic[g, j]),
                "fraction_group": float(fractions[g, j]),
                "fraction_rest": float((x[~masks[g], j] > 0).mean()),
                "mean_normalized_group": float(means[g, j]),
                "mean_normalized_rest": float(
                    (means[g, j] + 1e-9) / np.exp2(arithmetic[g, j]) - 1e-9
                ),
            }
            for j in np.flatnonzero(material_significant)
        )
        sensitivity.append({
            "group": group,
            "cells": int(masks[g].sum()),
            "rest_cells": int((~masks[g]).sum()),
            "significant_default": int((q[g] < 0.05).sum()),
            "significant_tie": int((q_tie[g] < 0.05).sum()),
            "new_significant_tie": int(((q_tie[g] < 0.05) & (q[g] >= 0.05)).sum()),
            "top20_overlap_tie": len(top_default & top_tie),
            "p_zero_default": int((p[g] == 0).sum()),
            "p_zero_tie": int((p_tie[g] == 0).sum()),
            "q_zero_default": int((q[g] == 0).sum()),
            "q_zero_tie": int((q_tie[g] == 0).sum()),
            "lfc_abs_delta_median": float(np.median(abs(delta))),
            "lfc_abs_delta_p95": float(np.quantile(abs(delta), 0.95)),
            "lfc_abs_delta_max": float(abs(delta).max()),
            "lfc_sign_reversals": int(reversal.sum()),
            "lfc_reversals_at_least_one_abs_lte1": int((reversal & ~material).sum()),
            "lfc_reversals_both_abs_gt1": int(material.sum()),
            "lfc_reversals_both_abs_gt1_default_q_lt005": int(
                material_significant.sum()
            ),
            "positive_lfc_gt1_default": int(((q[g] < 0.05) & (lfc[g] > 1)).sum()),
            "positive_lfc_gt1_arithmetic": int(
                ((q[g] < 0.05) & (arithmetic[g] > 1)).sum()
            ),
            "top5_default": result["names"][group][:5].tolist(),
            "top5_tie": a.uns["audit_tie"]["names"][group][:5].tolist(),
        })
        for marker in marker_names:
            j = genes.get_loc(marker)
            marker_rows.append({
                "group": group,
                "gene": marker,
                "rank": int(np.flatnonzero(result["names"][group] == marker)[0] + 1),
                "score": float(z[g, j]),
                "q": float(q[g, j]),
                "log2fc_default": float(lfc[g, j]),
                "log2fc_arithmetic": float(arithmetic[g, j]),
                "fraction_expressing": float(fractions[g, j]),
                "mean_normalized_count": float(means[g, j]),
            })
    print("Check dotplot means and fractions", flush=True)
    plot = sc.pl.rank_genes_groups_dotplot(
        a, groupby="leiden_res_0.50", standard_scale="var", n_genes=5, return_fig=True
    )
    plot_genes = plot.dot_color_df.columns
    pi = genes.get_indexer(plot_genes)
    logmeans = np.empty((k, len(pi)))
    for g, mask in enumerate(masks):
        logmeans[g] = np.asarray(x[mask][:, pi].astype(np.float64).mean(axis=0)).ravel()
    plot_groups = plot.dot_color_df.index.astype(str)
    gi = [groups.index(group) for group in plot_groups]
    scaled = logmeans - logmeans.min(axis=0)
    scaled /= scaled.max(axis=0)
    plot_checks = {
        "mean_colors": compare(plot.dot_color_df.to_numpy(), scaled[gi]),
        "expression_fractions": compare(
            plot.dot_size_df.to_numpy(), fractions[np.ix_(gi, pi)]
        ),
    }
    plt.close("all")
    pd.DataFrame(sensitivity).to_csv(evidence_dir / "10-sensitivity.csv", index=False)
    pd.DataFrame(reversal_examples).to_csv(
        evidence_dir / "10-foldchange-reversals.csv", index=False
    )
    pd.DataFrame(marker_rows).to_csv(evidence_dir / "10-markers.csv", index=False)
    sample_counts = pd.crosstab(a.obs["leiden_res_0.50"], a.obs["sample"])
    sample_counts.to_csv(evidence_dir / "10-sample-counts.csv")
    evidence = {
        "shape": [n, m],
        "comparisons": k * m,
        "groups": groups,
        "params": params,
        "head_at_audit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "scanpy_version": sc.__version__,
        "raw_present": a.raw is not None,
        "layers": list(a.layers),
        "highly_variable_count": int(a.var["highly_variable"].sum()),
        "matrix_dtype": str(x.dtype),
        "log1p": a.uns["log1p"],
        "chunk_size": args.chunk_size,
        "tie_coefficient_quantiles": np.quantile(tc, [0, 0.25, 0.5, 0.75, 1]).tolist(),
        "constant_genes": genes[tc == 0].tolist(),
        "checks": checks,
        "plot_checks": plot_checks,
        "numeric_counts": {
            name: {
                "nan": int(np.isnan(arr).sum()),
                "inf": int(np.isinf(arr).sum()),
                "zero": int((arr == 0).sum()),
            }
            for name, arr in {
                "z": z,
                "p": p,
                "q": q,
                "lfc": lfc,
                "p_tie": p_tie,
                "q_tie": q_tie,
                "lfc_arithmetic": arithmetic,
            }.items()
        },
        "zero_p_log10_tail_range": (
            stats.norm.logsf(abs(z[p == 0])) / np.log(10) + np.log10(2)
        ).tolist(),
        "sensitivity": sensitivity,
        "sample_counts": sample_counts.to_dict(),
    }
    evidence["resolved_defaults"] = {
        "preset": "ScanpyV1",
        "mask_var": None,
        "mean_in_log_space": True,
        "tie_correct": False,
        "rankby_abs": False,
        "n_genes": m,
        "continuity_correction": False,
        "alternative": "two-sided",
    }
    evidence["versions"] = {
        name: version(name)
        for name in ["numpy", "scipy", "pandas", "anndata", "numba", "statsmodels"]
    }
    evidence["execution"] = {
        str(stage): json.loads(
            (args.output / f"metrics/{stage:02}-execution.json").read_text()
        )
        for stage in range(6, 11)
    }
    evidence["runtime_failures"] = []
    evidence["actual_numeric_counts"] = {
        key: {
            field: {
                "nan": sum(
                    int(np.isnan(a.uns[key][field][group]).sum()) for group in groups
                ),
                "inf": sum(
                    int(np.isinf(a.uns[key][field][group]).sum()) for group in groups
                ),
                "zero": sum(
                    int((a.uns[key][field][group] == 0).sum()) for group in groups
                ),
            }
            for field in ["scores", "pvals", "pvals_adj", "logfoldchanges"]
        }
        for key in ["rank_genes_groups", "audit_tie", "audit_arithmetic"]
    }
    underflows = []
    for method, scores, pvalues in [("default", z, p), ("tie_correct", z_tie, p_tie)]:
        for g, j in zip(*np.where(pvalues == 0), strict=True):
            underflows.append({
                "method": method,
                "group": groups[g],
                "gene": genes[j],
                "score": scores[g, j],
                "log10_two_sided_p": float(
                    (stats.norm.logsf(abs(scores[g, j])) + np.log(2)) / np.log(10)
                ),
            })
    pd.DataFrame(underflows).to_csv(evidence_dir / "10-underflow.csv", index=False)
    # Store compact extrema in JSON and all affected gene/group pairs in CSV.
    tails = evidence.pop("zero_p_log10_tail_range")
    evidence["zero_p_log10_tail_min_max"] = [min(tails), max(tails)] if tails else []
    (evidence_dir / "10.json").write_text(json.dumps(evidence, indent=2) + "\n")
    for group_checks in [*checks.values(), plot_checks]:
        for value in group_checks.values():
            assert value["mismatch_count"] == 0 if isinstance(value, dict) else value
    print(
        json.dumps(
            {
                "comparisons": k * m,
                "numeric_counts": evidence["numeric_counts"],
                "plot_checks": plot_checks,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

"""Audit PCA arithmetic and associations on every tutorial cell."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import linalg, stats

import scanpy as sc
from scanpy.tools._utils import _choose_representation_compat


def correlations(scores, values):
    """Return Pearson and Spearman correlations for all components."""
    return {
        "pearson": [
            float(stats.pearsonr(scores[:, j], values).statistic) for j in range(50)
        ],
        "spearman": [
            float(stats.spearmanr(scores[:, j], values).statistic) for j in range(50)
        ],
    }


def group_metrics(mask, covariance, loadings, variance):
    """Attribute variance through covariance, retaining cross-gene terms."""
    # Cov(X_g v_g, X v) / Var(X v) sums to one across gene groups.
    contributions = loadings * (covariance @ loadings) / variance
    fractions = contributions[mask].sum(axis=0)
    return {
        "genes": int(mask.sum()),
        "input_variance_fraction": float(
            np.diag(covariance)[mask].sum() / np.trace(covariance)
        ),
        "squared_loading_mass_by_pc": (loadings[mask] ** 2).sum(axis=0).tolist(),
        "score_covariance_fraction_by_pc": fractions.tolist(),
        "retained_variance_fraction": float(fractions @ variance / variance.sum()),
    }


def main():  # noqa: PLR0915
    """Compare saved PCA with independent float64 calculations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = ad.read_h5ad(args.output / "checkpoints/05.h5ad")
    after = ad.read_h5ad(args.output / "checkpoints/06.h5ad")
    assert before.obs_names.equals(after.obs_names)
    assert before.var_names.equals(after.var_names)
    assert (before.X != after.X).nnz == 0
    assert (before.layers["counts"] != after.layers["counts"]).nnz == 0
    pd.testing.assert_frame_equal(before.obs, after.obs)
    pd.testing.assert_frame_equal(before.var, after.var)
    mask = before.var.highly_variable.to_numpy()
    x = before.X[:, mask].astype(np.float64)
    n, p = x.shape
    assert (n, p) == (17041, 2000)
    v = after.varm["PCs"][mask].astype(np.float64)
    scores = after.obsm["X_pca"].astype(np.float64)
    variance = after.uns["pca"]["variance"].astype(np.float64)
    ratios = after.uns["pca"]["variance_ratio"].astype(np.float64)
    assert scores.shape == (n, 50)
    assert after.varm["PCs"].shape == (before.n_vars, 50)
    assert np.count_nonzero(after.varm["PCs"][~mask]) == 0
    mean = np.asarray(x.mean(axis=0)).ravel()
    projected = x @ v - mean @ v
    covariance = ((x.T @ x).toarray() - n * np.outer(mean, mean)) / (n - 1)
    score_cov = np.cov(scores, rowvar=False)
    residual = covariance @ v - v * variance
    eigenvalues, eigenvectors = linalg.eigh(
        covariance, subset_by_index=(p - 51, p - 1), driver="evr"
    )
    eigenvalues, eigenvectors = eigenvalues[::-1], eigenvectors[:, ::-1]
    independent = eigenvectors[:, :50]
    signs = np.sign(np.sum(independent * v, axis=0))
    aligned = independent * signs
    angles = linalg.subspace_angles(independent, v)
    projection_error = float(np.max(abs(projected - scores)))
    np.testing.assert_allclose(projected, scores, atol=2e-4, rtol=1e-4)
    np.testing.assert_allclose(v.T @ v, np.eye(50), atol=2e-6)
    np.testing.assert_allclose(np.diag(score_cov), variance, rtol=1e-5)
    np.testing.assert_allclose(ratios, variance / np.trace(covariance), rtol=1e-5)
    np.testing.assert_allclose(variance, eigenvalues[:50], rtol=1e-5)
    assert float(np.max(angles)) < 0.001
    replay = before.copy()
    sc.tl.pca(replay)
    repeat_error = float(np.max(abs(replay.obsm["X_pca"] - after.obsm["X_pca"])))
    np.testing.assert_allclose(
        replay.obsm["X_pca"], after.obsm["X_pca"], atol=2e-4, rtol=1e-4
    )
    chosen = _choose_representation_compat(after, use_rep=None, n_pcs=None)
    np.testing.assert_array_equal(chosen, after.obsm["X_pca"])
    counts = before.layers["counts"][:, mask].astype(np.float64)
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    samples = before.obs["sample"].astype(str).to_numpy()
    batches = sorted(set(samples))
    absent = np.zeros(p, dtype=bool)
    for batch in batches:
        absent |= np.asarray(counts[samples == batch].sum(axis=0)).ravel() == 0
    selected = before.var.loc[mask]
    groups = {
        "detected_at_most_10_cells": detected <= 10,
        "absent_in_one_sample": absent,
        "mt": selected.mt.to_numpy(),
        "hb": selected.hb.to_numpy(),
        "high_mean_hvg_means_over_3": selected.means.to_numpy() > 3,
    }
    group_results = {
        name: group_metrics(m, covariance, v, variance) for name, m in groups.items()
    }
    associations = {
        name: correlations(scores, before.obs[name].to_numpy())
        for name in [
            "total_counts",
            "log1p_total_counts",
            "pct_counts_mt",
            "pct_counts_hb",
            "n_genes_by_counts",
        ]
    }
    associations["sample_s1d3_indicator"] = correlations(
        scores, (samples == "s1d3").astype(float)
    )
    associations["sample_s1d3_indicator"]["r_squared_by_pc"] = (
        np.asarray(associations["sample_s1d3_indicator"]["pearson"]) ** 2
    ).tolist()
    gene_results = {}
    for gene in selected.index[groups["high_mean_hvg_means_over_3"]]:
        j = selected.index.get_loc(gene)
        gene_results[gene] = {
            "log_expression_correlations": correlations(
                scores, x[:, j].toarray().ravel()
            ),
            "loading_by_pc": v[j].tolist(),
            "retained_variance_fraction": float(
                (v[j] ** 2) @ variance / variance.sum()
            ),
        }
    top_loadings = {}
    for pc in range(6):
        top = np.argsort(abs(v[:, pc]))[-10:][::-1]
        top_loadings[str(pc + 1)] = [
            {"gene": selected.index[j], "loading": float(v[j, pc])} for j in top
        ]
    centered_scores = scores - scores.mean(axis=0)
    ss_total = float(np.sum(centered_scores**2))
    sample_ss = sum(
        int((samples == batch).sum())
        * float(np.sum(centered_scores[samples == batch].mean(axis=0) ** 2))
        for batch in batches
    )
    evidence = {
        "baseline": "a6f1a2d2",
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["scanpy", "numpy", "scipy", "scikit-learn", "anndata"]
        },
        "shape": list(before.shape),
        "selected_genes": p,
        "selected_gene_list_sha256": hashlib.sha256(
            ("\n".join(selected.index) + "\n").encode()
        ).hexdigest(),
        "input_and_annotations_unchanged": True,
        "outside_hvg_loadings_nonzero": 0,
        "params": after.uns["pca"]["params"],
        "dtypes": {
            "input": str(before.X.dtype),
            "scores": str(after.obsm["X_pca"].dtype),
            "loadings": str(after.varm["PCs"].dtype),
            "variance": str(after.uns["pca"]["variance"].dtype),
        },
        "total_hvg_variance": float(np.trace(covariance)),
        "variance": variance.tolist(),
        "variance_ratio": ratios.tolist(),
        "cumulative_variance_ratio": np.cumsum(ratios).tolist(),
        "score_projection_max_absolute_error": projection_error,
        "score_projection_relative_frobenius_error": float(
            linalg.norm(projected - scores) / linalg.norm(scores)
        ),
        "loading_orthonormality_max_absolute_error": float(
            np.max(abs(v.T @ v - np.eye(50)))
        ),
        "score_mean_max_absolute": float(np.max(abs(scores.mean(axis=0)))),
        "score_variance_max_relative_error": float(
            np.max(abs(np.diag(score_cov) / variance - 1))
        ),
        "score_covariance_max_off_diagonal": float(
            np.max(abs(score_cov - np.diag(np.diag(score_cov))))
        ),
        "covariance_eigenvector_max_relative_residual": float(
            np.max(linalg.norm(residual, axis=0) / variance)
        ),
        "variance_ratio_max_absolute_error": float(
            np.max(abs(ratios - variance / np.trace(covariance)))
        ),
        "independent_eigh": {
            "max_eigenvalue_relative_error": float(
                np.max(abs(variance / eigenvalues[:50] - 1))
            ),
            "max_principal_angle_radians": float(np.max(angles)),
            "sign_aligned_loading_max_absolute_error": float(np.max(abs(aligned - v))),
            "minimum_relative_adjacent_gap_first_51": float(
                np.min((eigenvalues[:-1] - eigenvalues[1:]) / eigenvalues[:-1])
            ),
            "pc50_to_pc51_relative_gap": float(
                (eigenvalues[49] - eigenvalues[50]) / eigenvalues[49]
            ),
        },
        "default_replay_score_max_absolute_error": repeat_error,
        "default_neighbors_representation_shape": list(chosen.shape),
        "sample_between_group_fraction_of_retained_variance": sample_ss / ss_total,
        "associations": associations,
        "gene_groups": group_results,
        "high_mean_genes": gene_results,
        "top_absolute_loadings": top_loadings,
        "assertion_tolerances": {
            "projection_atol": 2e-4,
            "projection_rtol": 1e-4,
            "orthonormality_atol": 2e-6,
            "variance_rtol": 1e-5,
            "max_subspace_angle_radians": 0.001,
        },
    }
    dest = Path(__file__).resolve().parents[2] / "audit/tutorial/evidence/06.json"
    dest.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()

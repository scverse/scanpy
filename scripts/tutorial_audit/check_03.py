"""Audit real tutorial Scrublet results and captured neighbor score calculations."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from importlib.metadata import version
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from skimage.filters import threshold_minimum

import scanpy as sc
from scanpy.preprocessing import _scrublet as wrapper
from scanpy.preprocessing._scrublet import core, pipeline


def probability(nd, n, rho, ratio):
    """Calculate probability through smoothed, prior-adjusted doublet odds."""
    odds = (nd + 1) / (n - nd + 1) * rho / (ratio * (1 - rho))
    return odds / (1 + odds)


def main():  # noqa: PLR0915
    """Audit the saved tutorial and a graph-preserving sensitivity calculation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = ad.read_h5ad(args.output / "checkpoints/02.h5ad")
    after = ad.read_h5ad(args.output / "checkpoints/03.h5ad")
    assert (before.X != after.X).nnz == 0
    assert before.obs_names.equals(after.obs_names)
    assert before.var_names.equals(after.var_names)
    pd.testing.assert_frame_equal(before.var, after.var)
    # Scrublet refreshes these filter statistics inside each sample.
    old_columns = before.obs.columns.difference(["n_genes", "n_counts"])
    pd.testing.assert_frame_equal(before.obs[old_columns], after.obs[old_columns])
    evidence = {
        "counts_unchanged": True,
        "cell_gene_order_unchanged": True,
        "original_obs_unchanged_except_filter_statistics": True,
        "n_genes_changed_cells": int((before.obs.n_genes != after.obs.n_genes).sum()),
        "batches": {},
    }
    source_diff = subprocess.check_output(
        [
            "git",
            "diff",
            "a6f1a2d2",
            "--",
            "src/scanpy/preprocessing/_scrublet",
            "src/scanpy/neighbors",
            "src/scanpy/preprocessing/_utils.py",
            "src/scanpy/_utils/random.py",
        ],
        text=True,
    )
    assert source_diff == ""
    evidence["relevant_source_matches_baseline"] = True
    evidence["upstream_reference_commit"] = "67f8ecbad14e8e1aa9c89b43dac6638cebe38640"
    assert np.isfinite(before.X.data).all()
    assert (before.X.data >= 0).all()
    assert np.array_equal(before.X.data, np.rint(before.X.data))
    evidence.update({
        "baseline": "a6f1a2d2",
        "rerun_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "shape": list(before.shape),
        "raw_nonnegative_integer_counts": True,
        "versions": {
            name: version(name)
            for name in [
                "scanpy",
                "numpy",
                "scipy",
                "scikit-learn",
                "pynndescent",
                "scikit-image",
            ]
        },
        "checkpoint_sha256": {
            str(i): hashlib.file_digest(
                (args.output / f"checkpoints/{i:02}.h5ad").open("rb"), "sha256"
            ).hexdigest()
            for i in [2, 3]
        },
    })
    captured = []
    preprocessing = []
    simulate = wrapper.scrublet_simulate_doublets
    zscore = pipeline.zscore
    pca = pipeline.pca

    def capture_simulation(adata, **kwargs):
        result = simulate(adata, **kwargs)
        parents = result.obsm["doublet_parents"]
        raw = adata.layers["raw"]
        assert (raw[parents[:, 0]] + raw[parents[:, 1]] != result.X).nnz == 0
        batch = str(adata.obs["sample"].iloc[0])
        original = before[before.obs["sample"] == batch]
        gene_mask = np.asarray((original.X > 0).sum(axis=0)).ravel() >= 3
        cell_mask = np.asarray((original.X[:, gene_mask] > 0).sum(axis=1)).ravel() >= 3
        assert adata.n_obs == cell_mask.sum()
        assert np.array_equal(
            adata.obs_names.astype(int),
            np.flatnonzero(before.obs["sample"] == batch)[cell_mask],
        )
        assert (raw != original[cell_mask, adata.var_names].X).nnz == 0
        preprocessing.append({
            "batch": batch,
            "genes_passing_min_cells": int(gene_mask.sum()),
            "cells_passing_min_genes": int(cell_mask.sum()),
            "hvg_count": adata.n_vars,
            "raw_parent_sum_exact": True,
            "raw_hvg_input_exact": True,
            "same_parent_pairs": int((parents[:, 0] == parents[:, 1]).sum()),
            "unique_ordered_parent_pairs": len(np.unique(parents, axis=0)),
        })
        return result

    def capture_zscore(self):
        obs = self._counts_obs_norm.toarray().astype(np.float64)
        sim = self._counts_sim_norm.toarray().astype(np.float64)
        np.testing.assert_allclose(obs.sum(axis=1), 1e6, rtol=1e-6)
        np.testing.assert_allclose(sim.sum(axis=1), 1e6, rtol=1e-6)
        mean = obs.mean(axis=0)
        std = obs.std(axis=0, ddof=1)
        zscore(self)
        observed_error = np.max(
            np.abs(self._counts_obs_norm.toarray() - (obs - mean) / std)
        )
        simulated_error = np.max(
            np.abs(self._counts_sim_norm.toarray() - (sim - mean) / std)
        )
        assert observed_error < 1e-10
        assert simulated_error < 1e-10
        preprocessing[-1].update({
            "normalization_target_checked": 1e6,
            "observed_zscore_max_error": float(observed_error),
            "simulated_zscore_max_error": float(simulated_error),
        })

    def capture_pca(self, **kwargs):
        from sklearn.decomposition import PCA

        fit = PCA.fit
        fitted = []

        def capture_fit(estimator, x, *args, **kwargs):
            result = fit(estimator, x, *args, **kwargs)
            fitted.append(estimator)
            assert x.shape[0] == self._counts_obs_norm.shape[0]
            return result

        with patch.object(PCA, "fit", capture_fit):
            pca(self, **kwargs)
        estimator = fitted[0]
        errors = []
        for label in ["obs", "sim"]:
            x = getattr(self, f"_counts_{label}_norm").toarray()
            expected = (x - estimator.mean_) @ estimator.components_.T
            errors.append(
                float(np.max(np.abs(expected - getattr(self, f"manifold_{label}_"))))
            )
        assert max(errors) < 1e-8
        preprocessing[-1].update({
            "pca_fit_observed_only": True,
            "pca_components": estimator.n_components_,
            "pca_projection_max_errors": errors,
        })

    extract = core._get_indices_distances_from_sparse_matrix
    calculate = core.Scrublet.calculate_doublet_scores

    def capture_graph(matrix, k):
        indices, distances = extract(matrix, k)
        # Inspect stored rows directly, without Scanpy extraction for the correction.
        widths = np.diff(matrix.indptr)
        assert np.all(widths == k + 1)
        stored_indices = matrix.indices.reshape(matrix.shape[0], k + 1).copy()
        stored_distances = matrix.data.reshape(matrix.shape[0], k + 1).copy()
        assert np.all(
            (stored_indices == np.arange(matrix.shape[0])[:, None]).sum(axis=1) == 1
        )
        assert np.all(stored_distances[:, 0] == 0)
        assert np.all(stored_distances >= 0)
        assert np.all(np.diff(stored_distances, axis=1) >= 0)
        assert np.array_equal(indices, stored_indices[:, :k])
        captured.append({
            "indices": indices.copy(),
            "stored_indices": stored_indices,
            "stored_distances": stored_distances,
        })
        return indices, distances

    def capture_scores(self, **kwargs):
        result = calculate(self, **kwargs)
        captured[-1]["obs"] = self.manifold_obs_.copy()
        captured[-1]["sim"] = self.manifold_sim_.copy()
        captured[-1]["scores"] = np.r_[
            self.doublet_scores_obs_, self.doublet_scores_sim_
        ]
        captured[-1]["errors"] = np.r_[
            self.doublet_errors_obs_, self.doublet_errors_sim_
        ]
        return result

    sc.settings.n_jobs = 8
    with (
        patch.object(core, "_get_indices_distances_from_sparse_matrix", capture_graph),
        patch.object(core.Scrublet, "calculate_doublet_scores", capture_scores),
        patch.object(wrapper, "scrublet_simulate_doublets", capture_simulation),
        patch.object(pipeline, "zscore", capture_zscore),
        patch.object(pipeline, "pca", capture_pca),
    ):
        sc.pp.scrublet(before, batch_key="sample")
    np.testing.assert_allclose(
        before.obs.doublet_score, after.obs.doublet_score, rtol=0, atol=0
    )
    assert np.array_equal(before.obs.predicted_doublet, after.obs.predicted_doublet)
    evidence["instrumented_rerun_exactly_matches"] = True
    for (batch, meta), graph in zip(
        after.uns["scrublet"]["batches"].items(), captured, strict=True
    ):
        selected = after.obs["sample"] == batch
        obs = after.obs.loc[selected, "doublet_score"].to_numpy()
        sim = meta["doublet_scores_sim"]
        np.testing.assert_array_equal(sim, graph["scores"][len(obs) :])
        np.testing.assert_array_equal(
            meta["doublet_parents"],
            before.uns["scrublet"]["batches"][batch]["doublet_parents"],
        )
        assert preprocessing[len(evidence["batches"])]["batch"] == batch
        threshold = float(meta["threshold"])
        assert np.isfinite(obs).all()
        assert np.isfinite(sim).all()
        assert ((obs > 0) & (obs < 1)).all()
        assert ((sim > 0) & (sim < 1)).all()
        assert np.array_equal(
            obs > threshold, after.obs.loc[selected, "predicted_doublet"]
        )
        assert threshold == threshold_minimum(sim)
        assert meta["doublet_parents"].shape == (2 * len(obs), 2)
        assert meta["doublet_parents"].min() >= 0
        assert meta["doublet_parents"].max() < len(obs)
        indices = graph["indices"]
        n = indices.shape[1]
        nd = (indices >= len(obs)).sum(axis=1)
        rho = meta["parameters"]["expected_doublet_rate"]
        ratio = len(sim) / len(obs)
        independent = probability(nd, n, rho, ratio)
        error = float(np.max(np.abs(independent - graph["scores"])))
        assert error < 1e-13
        self_count = int(
            (indices == np.arange(len(indices))[:, None]).any(axis=1).sum()
        )
        assert self_count == len(indices)
        # Independent delta-method derivatives of the odds-based probability.
        q = (nd + 1) / (n + 2)
        se_q = np.sqrt(q * (1 - q) / (n + 3))
        se = (
            independent
            * (1 - independent)
            * np.sqrt((se_q / (q * (1 - q))) ** 2 + (0.02 / (rho * (1 - rho))) ** 2)
        )
        se_error = float(np.max(np.abs(se - graph["errors"])))
        assert se_error < 1e-13
        # Controlled correction: keep k and graph fixed, replace self with the next stored neighbor.
        stored = graph["stored_indices"]
        corrected_indices = stored[stored != np.arange(len(stored))[:, None]].reshape(
            len(stored), n
        )
        assert corrected_indices.shape == indices.shape
        corrected = probability(
            (corrected_indices >= len(obs)).sum(axis=1), n, rho, ratio
        )
        corrected_threshold = float(threshold_minimum(corrected[len(obs) :]))
        corrected_delta = corrected[: len(obs)] - obs
        assert corrected_delta.min() > -1e-13
        explicit_true = probability(
            (indices[:, 1:] >= len(obs)).sum(axis=1), n, rho, ratio
        )
        graph["corrected_scores"] = corrected
        # Exercise both real classifier branches on the captured graph, without a new search.
        replay = core.Scrublet(
            np.zeros((len(obs), 1)), n_neighbors=round(n / 3), expected_doublet_rate=rho
        )
        replay.set_manifold(graph["obs"], graph["sim"])
        with (
            patch.object(core.Neighbors, "compute_neighbors"),
            patch.object(
                core,
                "_get_indices_distances_from_sparse_matrix",
                return_value=(indices, None),
            ),
        ):
            replay.calculate_doublet_scores(use_approx_neighbors=True)
        np.testing.assert_allclose(
            np.r_[replay.doublet_scores_obs_, replay.doublet_scores_sim_],
            explicit_true,
            atol=1e-13,
            rtol=0,
        )
        with (
            patch.object(core.Neighbors, "compute_neighbors"),
            patch.object(
                core,
                "_get_indices_distances_from_sparse_matrix",
                return_value=(corrected_indices, None),
            ),
        ):
            replay.calculate_doublet_scores()
        np.testing.assert_allclose(
            np.r_[replay.doublet_scores_obs_, replay.doublet_scores_sim_],
            corrected,
            atol=1e-13,
            rtol=0,
        )
        # Hold the other neighbors fixed. This measures self-vote sensitivity,
        # not the full upstream k-other-neighbor correction.
        nd_without = nd - (np.arange(len(indices)) >= len(obs))
        without = probability(nd_without, n - 1, rho, ratio)
        threshold_without = float(threshold_minimum(without[len(obs) :]))
        changed_fixed = np.flatnonzero(
            (without[: len(obs)] > threshold) != (obs > threshold)
        )
        changed_rethreshold = np.flatnonzero(
            (without[: len(obs)] > threshold_without) != (obs > threshold)
        )
        # Independent exact Euclidean search for 128 real observed cells.
        rows = np.unique(np.linspace(0, len(obs) - 1, 128, dtype=int))
        manifold = np.vstack([graph["obs"], graph["sim"]])
        distances = cdist(manifold[rows], manifold)
        distances[np.arange(len(rows)), rows] = np.inf
        exact_indices = np.argpartition(distances, n - 1, axis=1)[:, :n]
        exact_scores = probability(
            (exact_indices >= len(obs)).sum(axis=1), n, rho, ratio
        )
        entry = {
            "preprocessing": preprocessing[len(evidence["batches"])],
            "score_error_formula_max_absolute_error": se_error,
            "stored_neighbors_including_self": graph["stored_indices"].shape[1],
            "expected_doublet_rate": float(rho),
            "stdev_doublet_rate": 0.02,
            "detectable_fraction": float((sim > threshold).mean()),
            "inferred_overall_rate": float(
                (obs > threshold).mean() / (sim > threshold).mean()
            ),
            "corrected_nonself_neighbors": n,
            "corrected_observed_mean_score_increase": float(corrected_delta.mean()),
            "corrected_observed_max_score_increase": float(corrected_delta.max()),
            "corrected_observed_increased_cells": int((corrected_delta > 1e-13).sum()),
            "corrected_calls_changed_fixed_threshold": int(
                ((corrected[: len(obs)] > threshold) != (obs > threshold)).sum()
            ),
            "corrected_threshold": corrected_threshold,
            "corrected_predicted_doublets": int(
                (corrected[: len(obs)] > corrected_threshold).sum()
            ),
            "corrected_calls_changed_recomputed_threshold": int(
                (
                    (corrected[: len(obs)] > corrected_threshold) != (obs > threshold)
                ).sum()
            ),
            "corrected_gained_cell_names": after.obs_names[selected][
                (corrected[: len(obs)] > corrected_threshold) & ~(obs > threshold)
            ].tolist(),
            "corrected_lost_cell_names": after.obs_names[selected][
                ~(corrected[: len(obs)] > corrected_threshold) & (obs > threshold)
            ].tolist(),
            "explicit_true_same_graph_observed_max_difference": float(
                np.max(np.abs(explicit_true[: len(obs)] - obs))
            ),
            "explicit_true_neighbor_count": n - 1,
            "explicit_true_formula_denominator_count": n,
            "observed_cells": len(obs),
            "simulated_cells": len(sim),
            "predicted_doublets": int((obs > threshold).sum()),
            "threshold": threshold,
            "observed_score_range": [float(obs.min()), float(obs.max())],
            "simulated_score_range": [float(sim.min()), float(sim.max())],
            "score_formula_max_absolute_error": error,
            "neighbors_including_self": n,
            "rows_including_self": self_count,
            "rows_with_self_outside_first_column": int(
                (indices[:, 0] != np.arange(len(indices))).sum()
            ),
            "self_vote_sensitivity_max_score_change": float(
                np.max(np.abs(without[: len(obs)] - obs))
            ),
            "self_vote_sensitivity_calls_changed_fixed_threshold": len(changed_fixed),
            "self_vote_sensitivity_calls_changed_recomputed_threshold": len(
                changed_rethreshold
            ),
            "self_vote_sensitivity_recomputed_threshold": threshold_without,
            "exact_other_neighbor_subset_cells": len(rows),
            "exact_other_neighbor_subset_max_score_change": float(
                np.max(np.abs(exact_scores - obs[rows]))
            ),
            "exact_other_neighbor_subset_calls_changed_fixed_threshold": int(
                ((exact_scores > threshold) != (obs[rows] > threshold)).sum()
            ),
        }
        entry["histogram_bin_sensitivity"] = {
            str(bins): {
                "threshold": float(cut),
                "predicted_doublets": int((obs > cut).sum()),
            }
            for bins in [64, 128, 256, 512]
            for cut in [threshold_minimum(sim, nbins=bins)]
        }
        entry["classifier_replays_match_independent_formulas"] = True
        entry["explicit_true_rows_still_containing_self"] = int(
            (indices[:, 1:] == np.arange(len(indices))[:, None]).any(axis=1).sum()
        )
        entry["exact_subset_corrected_neighbor_mean_recall"] = float(
            np.mean([
                len(np.intersect1d(exact_indices[j], corrected_indices[row])) / n
                for j, row in enumerate(rows)
            ])
        )
        fig, axes = plt.subplots(1, 2, figsize=(10, 3), constrained_layout=True)
        for axis, observed_scores, simulated_scores, cut, title in zip(
            axes,
            [obs, corrected[: len(obs)]],
            [sim, corrected[len(obs) :]],
            [threshold, corrected_threshold],
            ["Tutorial checkpoint", "Self-exclusion sensitivity"],
            strict=True,
        ):
            axis.hist(
                observed_scores,
                bins=60,
                density=True,
                histtype="step",
                label="Observed",
            )
            axis.hist(
                simulated_scores,
                bins=60,
                density=True,
                histtype="step",
                label="Simulated",
            )
            axis.axvline(cut, color="black", linestyle="--")
            axis.set(
                title=f"{batch}: {title}", xlabel="Doublet score", ylabel="Density"
            )
            axis.legend()
        fig.savefig(args.output / "plots" / f"03-{batch}-sensitivity.png", dpi=150)
        plt.close(fig)
        evidence["batches"][batch] = entry
        np.savez_compressed(
            args.output / "metrics" / f"03-{batch}-neighbor-audit.npz", **graph
        )
    dest = Path(__file__).resolve().parents[2] / "audit/tutorial/evidence/03.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()

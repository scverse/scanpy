"""Audit saved tutorial neighbors and UMAP against independent calculations."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist


def summary(values):
    """Describe the distribution without a parametric model."""
    values = np.asarray(values)
    return dict(
        zip(
            ["min", "p05", "median", "p95", "max"],
            np.quantile(values, [0, 0.05, 0.5, 0.95, 1]).tolist(),
            strict=True,
        )
    ) | {"mean": float(values.mean())}


def main():  # noqa: PLR0915
    """Measure graph accuracy and embedding distortion on saved cells."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifacts", type=Path, default=Path("/home/fdr/scanpy-tutorial-audit")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("audit/tutorial/evidence/07.json")
    )
    args = parser.parse_args()
    before = ad.read_h5ad(args.artifacts / "checkpoints/06.h5ad", backed="r")
    data = ad.read_h5ad(args.artifacts / "checkpoints/07.h5ad", backed="r")
    x = data.obsm["X_pca"].astype(np.float64)
    y = data.obsm["X_umap"].astype(np.float64)
    n = len(x)
    d = data.obsp["distances"].tocsr()
    c = data.obsp["connectivities"].tocsr()
    labels = data.obs["sample"].astype(str).to_numpy()
    groups = np.unique(labels)
    # Independently sort the saved distances, excluding self by row identity.
    indices = []
    distances = []
    for i in range(n):
        start, end = d.indptr[i : i + 2]
        js, ds = d.indices[start:end], d.data[start:end]
        order = np.argsort(ds[js != i], kind="stable")[:14]
        indices.append(js[js != i][order])
        distances.append(ds[js != i][order])
    indices, distances = np.array(indices), np.array(distances)
    rows = np.repeat(np.arange(n), np.diff(d.indptr))
    actual_distances = np.linalg.norm(x[rows] - x[d.indices], axis=1)
    nonself = rows != d.indices
    components, component_labels = connected_components(c, directed=False)
    cr = np.repeat(np.arange(n), np.diff(c.indptr))
    cross = labels[cr] != labels[c.indices]
    cross_weight = np.bincount(cr, weights=c.data * cross, minlength=n)
    total_weight = np.asarray(c.astype(np.float64).sum(axis=1)).ravel()
    mixing = cross_weight / total_weight

    # Independent float64 root solution of the UMAP local membership equation.
    rho = distances[:, 0].astype(np.float64)
    delta = np.maximum(0, distances.astype(np.float64) - rho[:, None])
    low, high = np.zeros(n), np.maximum(distances[:, -1], 1).astype(np.float64) * 100
    for _ in range(80):
        sigma = (low + high) / 2
        mass = np.exp(-delta / sigma[:, None]).sum(axis=1)
        high = np.where(mass > np.log2(15), sigma, high)
        low = np.where(mass <= np.log2(15), sigma, low)
    sigma = np.maximum((low + high) / 2, 0.001 * distances.sum(axis=1) / 15)
    membership = np.exp(-delta / sigma[:, None])
    directed = sparse.csr_matrix(  # noqa: TID251
        (membership.ravel(), (np.repeat(np.arange(n), 14), indices.ravel())),
        shape=(n, n),
    )
    rebuilt = directed + directed.T - directed.multiply(directed.T)
    formula_error = rebuilt - c
    embedding_graph = c.copy()
    pruned = embedding_graph.data < embedding_graph.data.max() / 200
    pruned_count = int(pruned.sum())
    embedding_graph.data[pruned] = 0
    embedding_graph.eliminate_zeros()
    embedding_components, embedding_labels = connected_components(
        embedding_graph, directed=False
    )

    # Equal sample strata, deterministic queries, all cells as reference candidates.
    rng = np.random.default_rng(703)
    query = np.sort(
        np.concatenate([
            rng.choice(np.flatnonzero(labels == group), 512, replace=False)
            for group in groups
        ])
    )
    recall, overlap, trust, continuity, inflation, missed_ranks = [], [], [], [], [], []
    for start in range(0, len(query), 128):
        q = query[start : start + 128]
        high_d = cdist(x[q], x)
        low_d = cdist(y[q], y)
        high_d[np.arange(len(q)), q] = np.inf
        low_d[np.arange(len(q)), q] = np.inf
        high_order = np.argsort(high_d, axis=1, kind="stable")
        low_order = np.argsort(low_d, axis=1, kind="stable")
        high_ranks = np.empty_like(high_order)
        low_ranks = np.empty_like(low_order)
        np.put_along_axis(high_ranks, high_order, np.arange(1, n + 1)[None, :], axis=1)
        np.put_along_axis(low_ranks, low_order, np.arange(1, n + 1)[None, :], axis=1)
        for j, row in enumerate(q):
            exact, embedded = high_order[j, :14], low_order[j, :14]
            approximate = indices[row]
            recall.append(len(np.intersect1d(exact, approximate)) / 14)
            overlap.append(len(np.intersect1d(exact, embedded)) / 14)
            trust.append(
                1
                - 2
                * np.maximum(high_ranks[j, embedded] - 14, 0).sum()
                / (14 * (2 * n - 3 * 14 - 1))
            )
            continuity.append(
                1
                - 2
                * np.maximum(low_ranks[j, exact] - 14, 0).sum()
                / (14 * (2 * n - 3 * 14 - 1))
            )
            inflation.append(high_d[j, approximate].max() / high_d[j, exact[-1]] - 1)
            missed_ranks.extend(
                high_ranks[j, approximate][high_ranks[j, approximate] > 14].tolist()
            )
    recall, overlap, trust = np.array(recall), np.array(overlap), np.array(trust)
    counts = {group: int((labels == group).sum()) for group in groups}
    evidence = {
        "baseline": "a6f1a2d2",
        "versions": {
            p: importlib.metadata.version(p)
            for p in [
                "scanpy",
                "umap-learn",
                "pynndescent",
                "numpy",
                "scipy",
                "scikit-learn",
            ]
        },
        "inputs": str(args.artifacts / "checkpoints"),
        "shape": list(data.shape),
        "pca_shape": list(x.shape),
        "neighbors_params": data.uns["neighbors"]["params"],
        "umap_params": data.uns["umap"]["params"],
        "row_identity": {
            "obs_names_equal": before.obs_names.equals(data.obs_names),
            "obs_equal": before.obs.equals(data.obs),
            "var_names_equal": before.var_names.equals(data.var_names),
            "pca_equal": bool(np.array_equal(before.obsm["X_pca"], data.obsm["X_pca"])),
            "unique_names": data.obs_names.is_unique,
            "obs_names_sha256": hashlib.sha256(
                "\n".join(data.obs_names).encode()
            ).hexdigest(),
        },
        "distances": {
            "nnz": d.nnz,
            "stored_per_row": summary(np.diff(d.indptr)),
            "positive_per_row": summary(
                np.bincount(rows, weights=d.data > 0, minlength=n)
            ),
            "stored_self_entries": int((~nonself).sum()),
            "nonself_zero_entries": int(((d.data == 0) & nonself).sum()),
            "finite": bool(np.isfinite(d.data).all()),
            "range": summary(d.data),
            "absolute_error_all_entries": summary(abs(d.data - actual_distances)),
            "relative_error_nonself": summary(
                abs(d.data[nonself] - actual_distances[nonself])
                / actual_distances[nonself]
            ),
        },
        "connectivities": {
            "nnz": c.nnz,
            "finite": bool(np.isfinite(c.data).all()),
            "weights": summary(c.data),
            "degree": summary(np.diff(c.indptr)),
            "diagonal_nonzero": int(np.count_nonzero(c.diagonal())),
            "asymmetric_entries": (c - c.T).nnz,
            "components": components,
            "component_sizes": sorted(
                np.bincount(component_labels).tolist(), reverse=True
            ),
            "independent_formula_support_difference": int(
                ((rebuilt != 0) != (c != 0)).nnz
            ),
            "embedding_pruned_directed_entries": pruned_count,
            "embedding_graph_components": embedding_components,
            "embedding_component_sizes": sorted(
                np.bincount(embedding_labels).tolist(), reverse=True
            ),
            "independent_formula_max_absolute_error": float(abs(formula_error).max()),
            "independent_formula_mean_absolute_error_over_union_support": float(
                abs(formula_error).sum() / rebuilt.nnz
            ),
        },
        "mixing": {
            "sample_counts": counts,
            "cross_sample_edge_fraction": float(cross.mean()),
            "cross_sample_weight_fraction": float(np.sum(c.data[cross]) / c.data.sum()),
            "per_cell_cross_weight": summary(mixing),
            "by_sample": {
                group: {
                    "cross_weight": summary(mixing[labels == group]),
                    "random_label_expected_cross_fraction": (n - counts[group])
                    / (n - 1),
                }
                for group in groups
            },
        },
        "evaluation": {
            "seed": 703,
            "sample_scheme": "512 cells without replacement per sample; sorted row indices; all 17041 cells as candidates",
            "query_indices": query.tolist(),
            "query_obs_names": data.obs_names[query].tolist(),
            "k_nonself": 14,
            "approximate_neighbor_recall": summary(recall),
            "perfect_recall_queries": int((recall == 1).sum()),
            "approximate_radius_relative_excess": summary(inflation),
            "nonexact_neighbor_true_ranks": summary(missed_ranks),
            "embedding_exact_neighbor_overlap": summary(overlap),
            "query_trustworthiness": summary(trust),
            "query_continuity": summary(continuity),
            "by_sample": {
                group: {
                    "recall_mean": float(recall[labels[query] == group].mean()),
                    "overlap_mean": float(overlap[labels[query] == group].mean()),
                    "trustworthiness_mean": float(trust[labels[query] == group].mean()),
                }
                for group in groups
            },
        },
        "embedding": {
            "shape": list(y.shape),
            "finite": bool(np.isfinite(y).all()),
            "unique_rows": len(np.unique(y, axis=0)),
            "axis_min": y.min(axis=0).tolist(),
            "axis_max": y.max(axis=0).tolist(),
        },
    }
    args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in evidence.items() if k != "evaluation"}, indent=2)
    )
    print(
        json.dumps(
            {
                k: v
                for k, v in evidence["evaluation"].items()
                if k not in ["query_indices", "query_obs_names"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

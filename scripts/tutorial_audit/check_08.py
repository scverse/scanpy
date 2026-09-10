"""Audit tutorial Leiden partitions and marker summaries on saved checkpoints."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.metadata
import json
import random
from pathlib import Path

import anndata as ad
import igraph as ig
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

import scanpy as sc


class LegacyIgraphRng:
    """Reproduce Scanpy's default legacy seed without its RNG helpers."""

    def __init__(self, seed):
        self.state = np.random.RandomState(seed)

    def random(self):
        """Return a uniform draw."""
        return self.state.random_sample()

    def randint(self, a, b):
        """Return an inclusive integer draw."""
        return self.state.randint(a, b + 1)

    def gauss(self, mu, sigma):
        """Return a normal draw."""
        return self.state.normal(mu, sigma)

    def getrandbits(self, k):
        """Return the legacy masked integer draw."""
        return self.state.tomaxint() & ((1 << k) - 1)


def partition(graph, resolution, iterations, seed=0):
    """Run igraph with the tutorial objective and RNG."""
    ig.set_random_number_generator(LegacyIgraphRng(seed))
    result = graph.community_leiden(
        objective_function="modularity",
        weights="weight",
        resolution=resolution,
        n_iterations=iterations,
        beta=0.01,
    )
    ig.set_random_number_generator(random)
    return result


def objective(matrix, labels, resolution):
    """Calculate generalized modularity directly from the symmetric matrix."""
    _, codes = np.unique(labels, return_inverse=True)
    coo = matrix.tocoo()
    total = coo.data.sum(dtype=np.float64)
    degree = np.asarray(matrix.astype(np.float64).sum(axis=1)).ravel()
    mass = np.bincount(codes, weights=degree) / total
    within = coo.data[codes[coo.row] == codes[coo.col]].sum(dtype=np.float64)
    return float(within / total - resolution * np.square(mass).sum())


def cluster_summary(obs, key):
    """Summarize sample composition and retained QC annotations."""
    result = {}
    for label in obs[key].cat.categories:
        group = obs.loc[obs[key] == label]
        result[label] = {
            "n": len(group),
            "samples": group["sample"].value_counts().to_dict(),
            "predicted_doublet_n": int(group["predicted_doublet"].sum()),
            "predicted_doublet_fraction": float(group["predicted_doublet"].mean()),
            "pct_counts_mt_mean": float(group["pct_counts_mt"].mean()),
            "pct_counts_mt_median": float(group["pct_counts_mt"].median()),
            "pct_counts_mt_p95": float(group["pct_counts_mt"].quantile(0.95)),
            "pct_counts_mt_over_20_fraction": float(
                (group["pct_counts_mt"] > 20).mean()
            ),
        }
    return result


def main():  # noqa: PLR0915
    """Compare saved results with independent calculations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifacts", type=Path, default=Path("/home/fdr/scanpy-tutorial-audit")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("audit/tutorial/evidence/08.json")
    )
    args = parser.parse_args()
    checkpoints = args.artifacts / "checkpoints"
    before = ad.read_h5ad(checkpoints / "07.h5ad", backed="r")
    clustered = ad.read_h5ad(checkpoints / "08.h5ad", backed="r")
    annotated = ad.read_h5ad(checkpoints / "09.h5ad", backed="r")
    matrix = before.obsp["connectivities"].tocsr()
    coo = matrix.tocoo()
    # Preserve CSR traversal order and both symmetric entries, independently of Scanpy.
    graph = ig.Graph(
        n=before.n_obs,
        edges=list(zip(coo.row.tolist(), coo.col.tolist(), strict=True)),
        directed=False,
    )
    graph.es["weight"] = coo.data.astype(np.float64).tolist()
    upper = sparse.triu(matrix, k=1).tocoo()
    simple = ig.Graph(
        n=before.n_obs,
        edges=list(zip(upper.row.tolist(), upper.col.tolist(), strict=True)),
        directed=False,
    )
    simple.es["weight"] = upper.data.astype(np.float64).tolist()
    working = ad.AnnData(obs=before.obs.copy())
    working.obsp["connectivities"] = matrix
    evidence = {
        "baseline": "a6f1a2d2",
        "versions": {
            p: importlib.metadata.version(p)
            for p in ["scanpy", "igraph", "numpy", "pandas", "anndata"]
        },
        "inputs": str(checkpoints),
        "shape": list(before.shape),
        "row_alignment": {
            "07_08_equal": before.obs_names.equals(clustered.obs_names),
            "08_09_equal": clustered.obs_names.equals(annotated.obs_names),
            "07_columns_preserved": before.obs.equals(
                clustered.obs[before.obs.columns]
            ),
            "08_columns_preserved": clustered.obs.equals(
                annotated.obs[clustered.obs.columns]
            ),
            "graph_07_08_equal": (matrix != clustered.obsp["connectivities"]).nnz == 0,
            "names_unique": before.obs_names.is_unique,
            "obs_names_sha256": hashlib.sha256(
                "\n".join(before.obs_names).encode()
            ).hexdigest(),
        },
        "graph": {
            "directed": graph.is_directed(),
            "vertices": graph.vcount(),
            "edges": graph.ecount(),
            "simple_edges": simple.ecount(),
            "parallel_edges": sum(graph.is_multiple()),
            "asymmetric_entries": (matrix - matrix.T).nnz,
            "self_edges": int(np.count_nonzero(matrix.diagonal())),
            "weight_min": float(coo.data.min()),
            "weight_max": float(coo.data.max()),
            "weights_finite": bool(np.isfinite(coo.data).all()),
        },
        "partitions": {},
    }
    specs = [
        ("leiden", 1.0, 2),
        ("leiden_res_0.02", 0.02, -1),
        ("leiden_res_0.50", 0.5, -1),
        ("leiden_res_2.00", 2.0, -1),
    ]
    for key, resolution, iterations in specs:
        saved = clustered.obs[key].astype(str).to_numpy()
        direct = partition(graph, resolution, iterations)
        repeated = partition(graph, resolution, iterations)
        single = partition(simple, resolution, iterations)
        sc.tl.leiden(
            working,
            adjacency=matrix,
            flavor="igraph",
            key_added=key,
            resolution=resolution,
            n_iterations=iterations,
        )
        formula = objective(matrix, saved, resolution)
        codes = saved.astype(int).tolist()
        evidence["partitions"][key] = {
            "params": clustered.uns[key]["params"],
            "categories": clustered.obs[key].cat.categories.tolist(),
            "cluster_count": clustered.obs[key].nunique(),
            "scanpy_exact_saved": bool(
                np.array_equal(saved, working.obs[key].astype(str))
            ),
            "direct_exact_saved": bool(
                np.array_equal(saved, np.array(direct.membership).astype(str))
            ),
            "repeat_exact": direct.membership == repeated.membership,
            "direct_ari_saved": adjusted_rand_score(saved, direct.membership),
            "objective_independent": formula,
            "objective_stored": float(clustered.uns[key]["modularity"]),
            "objective_direct": direct.quality,
            "objective_simple_same_partition": simple.modularity(
                codes, weights="weight", resolution=resolution
            ),
            "objective_duplicate_same_partition": graph.modularity(
                codes, weights="weight", resolution=resolution
            ),
            "simple_refit_ari": adjusted_rand_score(saved, single.membership),
            "simple_refit_objective": objective(matrix, single.membership, resolution),
            "clusters": cluster_summary(clustered.obs, key),
            "seed_sensitivity": {
                str(seed): {
                    "ari_saved": adjusted_rand_score(
                        saved,
                        (
                            part := partition(graph, resolution, iterations, seed)
                        ).membership,
                    ),
                    "n_clusters": len(part),
                    "objective": objective(matrix, part.membership, resolution),
                }
                for seed in [1, 2, 3]
            },
        }
    converged = partition(graph, 1, -1)
    evidence["iteration_sensitivity"] = {
        "resolution": 1,
        "seed": 0,
        "converged_clusters": len(converged),
        "ari_two_vs_converged": adjusted_rand_score(
            clustered.obs["leiden"], converged.membership
        ),
        "two_objective": objective(matrix, clustered.obs["leiden"], 1),
        "converged_objective": objective(matrix, converged.membership, 1),
    }
    evidence["resolution_cross_tabs"] = {
        f"{a}__{b}": pd.crosstab(clustered.obs[a], clustered.obs[b]).to_dict()
        for a, b in [
            ("leiden_res_0.02", "leiden_res_0.50"),
            ("leiden_res_0.50", "leiden_res_2.00"),
        ]
    }
    notebook_path = Path("docs/tutorials/basics/clustering.ipynb")
    notebook = json.loads(notebook_path.read_text())
    markers = ast.literal_eval(
        ast.parse("".join(notebook["cells"][59]["source"])).body[0].value
    )
    genes = [gene for group in markers.values() for gene in group]
    unique = list(dict.fromkeys(genes))
    marker_data = ad.AnnData(
        X=annotated[:, unique].X,
        obs=annotated.obs.copy(),
        var=annotated.var.loc[unique].copy(),
    )
    dense = marker_data.X.toarray().astype(np.float64)
    evidence["markers"] = {
        "notebook_sha256": hashlib.sha256(notebook_path.read_bytes()).hexdigest(),
        "genes_with_plot_duplicates": genes,
        "matrix": "checkpoint09.X (log1p normalized), no raw, no selected layer",
        "by_partition": {},
    }
    for key in ["leiden_res_0.02", "leiden_res_0.50"]:
        categories = annotated.obs[key].cat.categories.tolist()
        means = np.array([
            dense[(annotated.obs[key] == c).to_numpy()].mean(axis=0) for c in categories
        ])
        fractions = np.array([
            (dense[(annotated.obs[key] == c).to_numpy()] > 0).mean(axis=0)
            for c in categories
        ])
        span = np.ptp(means, axis=0)
        scaled = np.divide(
            means - means.min(axis=0), span, out=np.zeros_like(means), where=span != 0
        )
        dot = sc.pl.dotplot(
            marker_data, markers, groupby=key, standard_scale="var", return_fig=True
        )
        unscaled = sc.pl.dotplot(marker_data, markers, groupby=key, return_fig=True)
        ix = [unique.index(gene) for gene in genes]
        evidence["markers"]["by_partition"][key] = {
            "categories": categories,
            "unique_genes": unique,
            "means": means.tolist(),
            "nonzero_fractions": fractions.tolist(),
            "scaled_means": scaled.tolist(),
            "plot_row_order_equal": dot.dot_color_df.index.tolist() == categories,
            "plot_gene_order_equal": dot.dot_color_df.columns.tolist() == genes,
            "max_mean_error": float(
                np.max(abs(unscaled.dot_color_df.to_numpy() - means[:, ix]))
            ),
            "max_scaled_error": float(
                np.max(abs(dot.dot_color_df.to_numpy() - scaled[:, ix]))
            ),
            "max_fraction_error": float(
                np.max(abs(dot.dot_size_df.to_numpy() - fractions[:, ix]))
            ),
        }
    evidence["annotation"] = {
        "mapping": {
            c: (None if pd.isna(v) else v)
            for c, v in annotated.obs
            .groupby("leiden_res_0.02", observed=True)["cell_type_lvl1"]
            .first()
            .items()
        },
        "missing_cells": int(annotated.obs["cell_type_lvl1"].isna().sum()),
        "fraction_missing": float(annotated.obs["cell_type_lvl1"].isna().mean()),
    }
    de = ad.read_h5ad(checkpoints / "10.h5ad", backed="r")
    evidence["cluster7_dge"] = (
        sc.get.rank_genes_groups_df(de, group="7").head(10).to_dict(orient="records")
    )
    args.output.write_text(
        json.dumps(evidence, indent=2, default=lambda x: x.item()) + "\n"
    )
    print(
        json.dumps(
            {
                "annotation": evidence["annotation"],
                "iteration_sensitivity": evidence["iteration_sensitivity"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

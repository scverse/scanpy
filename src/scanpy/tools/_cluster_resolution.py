from __future__ import annotations

import sys
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData


def find_cluster_specific_genes(
    adata: AnnData,
    resolutions: Sequence[float],
    *,
    prefix: str = "leiden_res_",
    method: Literal["wilcoxon"] = "wilcoxon",
    n_top_genes: int = 3,
    min_cells: int = 2,
    deg_mode: Literal["within_parent", "per_resolution"] = "within_parent",
) -> dict[tuple[str, str], list[str]]:
    """Find differentially expressed genes for clusters in two modes."""
    from . import rank_genes_groups

    if deg_mode not in ["within_parent", "per_resolution"]:
        msg = "deg_mode must be 'within_parent' or 'per_resolution'"
        raise ValueError(msg)

    # Validate resolutions and clustering columns
    for res in resolutions:
        col = f"{prefix}{res}"
        if col not in adata.obs:
            msg = f"Column {col} not found in adata.obs"
            raise ValueError(msg)

    top_genes_dict: dict[tuple[str, str], list[str]] = {}

    if deg_mode == "within_parent":
        top_genes_dict.update(
            find_within_parent_degs(
                adata,
                resolutions,
                prefix=prefix,
                n_top_genes=n_top_genes,
                min_cells=min_cells,
                rank_genes_groups=rank_genes_groups,
            )
        )
    elif deg_mode == "per_resolution":
        top_genes_dict.update(
            find_per_resolution_degs(
                adata,
                resolutions,
                prefix=prefix,
                n_top_genes=n_top_genes,
                min_cells=min_cells,
                rank_genes_groups=rank_genes_groups,
            )
        )

    return top_genes_dict


def find_within_parent_degs(
    adata: AnnData,
    resolutions: Sequence[float],
    *,
    prefix: str,
    n_top_genes: int,
    min_cells: int,
    rank_genes_groups,
) -> dict[tuple[str, str], list[str]]:
    top_genes_dict = {}

    for i, res in enumerate(resolutions[:-1]):
        res_key = f"{prefix}{res}"
        next_res_key = f"{prefix}{resolutions[i + 1]}"
        clusters = adata.obs[res_key].cat.categories

        for cluster in clusters:
            cluster_mask = adata.obs[res_key] == cluster
            cluster_adata = adata[cluster_mask, :]

            subclusters = cluster_adata.obs[next_res_key].value_counts()
            valid_subclusters = subclusters[subclusters >= min_cells].index

            if len(valid_subclusters) < 2:
                print(
                    f"Skipping res_{res}_C{cluster}: < 2 subclusters with >= {min_cells} cells."
                )
                continue

            subcluster_mask = cluster_adata.obs[next_res_key].isin(valid_subclusters)
            deg_adata = cluster_adata[subcluster_mask, :]

            try:
                rank_genes_groups(deg_adata, groupby=next_res_key, method="wilcoxon")
                for subcluster in valid_subclusters:
                    names = deg_adata.uns["rank_genes_groups"]["names"][subcluster]
                    scores = deg_adata.uns["rank_genes_groups"]["scores"][subcluster]
                    top_genes = [
                        name
                        for name, score in zip(names, scores, strict=False)
                        if score > 0
                    ][:n_top_genes]
                    parent_node = f"res_{res}_C{cluster}"
                    child_node = f"res_{resolutions[i + 1]}_C{subcluster}"
                    top_genes_dict[(parent_node, child_node)] = top_genes
                    print(f"{parent_node} -> {child_node}: {top_genes}")
            except KeyError as e:
                print(f"Key error when processing {parent_node} -> {child_node}: {e}")
                continue
            except TypeError as e:
                print(
                    f"Type error with the data when processing {parent_node} -> {child_node}: {e}"
                )
                continue

    return top_genes_dict


def find_per_resolution_degs(
    adata: AnnData,
    resolutions: Sequence[float],
    *,
    prefix: str,
    n_top_genes: int,
    min_cells: int,
    rank_genes_groups,
) -> dict[tuple[str, str], list[str]]:
    top_genes_dict = {}

    for i, res in enumerate(resolutions[1:], 1):
        res_key = f"{prefix}{res}"
        prev_res_key = f"{prefix}{resolutions[i - 1]}"
        clusters = adata.obs[res_key].cat.categories
        valid_clusters = [
            c for c in clusters if (adata.obs[res_key] == c).sum() >= min_cells
        ]

        if not valid_clusters:
            print(f"Skipping resolution {res}: no clusters with >= {min_cells} cells.")
            continue

        deg_adata = adata[adata.obs[res_key].isin(valid_clusters), :]
        try:
            rank_genes_groups(
                deg_adata, groupby=res_key, method="wilcoxon", reference="rest"
            )
            for cluster in valid_clusters:
                names = deg_adata.uns["rank_genes_groups"]["names"][cluster]
                scores = deg_adata.uns["rank_genes_groups"]["scores"][cluster]
                top_genes = [
                    name
                    for name, score in zip(names, scores, strict=False)
                    if score > 0
                ][:n_top_genes]
                parent_cluster = adata.obs[deg_adata.obs[res_key] == cluster][
                    prev_res_key
                ].mode()[0]
                parent_node = f"res_{resolutions[i - 1]}_C{parent_cluster}"
                child_node = f"res_{res}_C{cluster}"
                top_genes_dict[(parent_node, child_node)] = top_genes
                print(f"{parent_node} -> {child_node}: {top_genes}")
        except KeyError as e:
            print(f"Key error when processing {parent_node} -> {child_node}: {e}")
            continue
        except TypeError as e:
            print(
                f"Type error with the data when processing {parent_node} -> {child_node}: {e}"
            )
            continue

    return top_genes_dict


def find_cluster_resolution(
    adata: AnnData,
    resolutions: list[float],
    *,
    prefix: str = "leiden_res_",
    method: Literal["wilcoxon"] = "wilcoxon",
    n_top_genes: int = 3,
    min_cells: int = 2,
    deg_mode: Literal["within_parent", "per_resolution"] = "within_parent",
    flavor: Literal["igraph"] = "igraph",
    n_iterations: int = 2,
) -> None:
    """
    Find clusters across multiple resolutions and identify cluster-specific genes.

    This function performs Leiden clustering at specified resolutions, identifies
    differentially expressed genes (DEGs) for clusters, and stores the results in `adata`.

    Params
    ------
    adata
        The annotated data matrix.
    resolutions
        List of resolution values for Leiden clustering (e.g., [0.0, 0.2, 0.5]).
    prefix
        Prefix for clustering keys in `adata.obs` (e.g., "leiden_res_").
    method
        Method for differential expression analysis: only "wilcoxon" is supported.
    n_top_genes
        Number of top genes to identify per child cluster.
    min_cells
        Minimum number of cells required in a subcluster to include it.
    deg_mode
        Mode for DEG analysis: "within_parent" (compare child to parent cluster) or
        "per_resolution" (compare within each resolution).
    flavor
        Flavor of Leiden clustering: only "igraph" is supported.
    n_iterations
        Number of iterations for Leiden clustering.

    Returns
    -------
    None

    The following annotations are added to `adata`:

    leiden_res_{resolution}
        Cluster assignments for each resolution in `adata.obs`.
    cluster_resolution_top_genes
        Dictionary mapping (parent_node, child_node) pairs to lists of top marker genes,
        stored in `adata.uns`.

    Notes
    -----
    This function requires the `igraph` library for Leiden clustering, which is included in the
    `leiden` extra. Install it with: ``pip install scanpy[leiden]``.

    Requires `sc.pp.neighbors` to be run on `adata` beforehand.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.pp.neighbors(adata)
    >>> sc.tl.find_cluster_resolution(adata, resolutions=[0.0, 0.5])
    >>> sc.pl.cluster_decision_tree(adata, resolutions=[0.0, 0.5])
    """
    import io

    from . import leiden

    # Suppress prints if pytest is running
    if "pytest" in sys.modules:
        sys.stdout = io.StringIO()

    _validate_cluster_resolution_inputs(adata, resolutions, method, flavor)

    # Run Leiden clustering
    for resolution in resolutions:
        res_key = f"{prefix}{resolution}"
        try:
            leiden(
                adata,
                resolution=resolution,
                flavor="igraph",
                n_iterations=n_iterations,
                key_added=res_key,
            )
            if "pytest" not in sys.modules and not hasattr(
                sys, "_called_from_test"
            ):  # Suppress print in tests
                print(f"Completed Leiden clustering for resolution {resolution}")
        except ValueError as e:
            msg = f"Leiden clustering failed at resolution {resolution} due to invalid value: {e}"
            raise RuntimeError(msg) from None
        except TypeError as e:
            msg = f"Leiden clustering failed at resolution {resolution} due to incorrect type: {e}"
            raise RuntimeError(msg) from None
        except RuntimeError as e:
            msg = f"Leiden clustering failed at resolution {resolution}: {e}"
            raise RuntimeError(msg) from None

    # Find cluster-specific genes
    top_genes_dict = find_cluster_specific_genes(
        adata=adata,
        resolutions=resolutions,
        prefix=prefix,
        method=method,
        n_top_genes=n_top_genes,
        min_cells=min_cells,
        deg_mode=deg_mode,
    )

    # Create DataFrame for clusterDecisionTree
    try:
        cluster_data = pd.DataFrame(
            {f"{prefix}{r}": adata.obs[f"{prefix}{r}"] for r in resolutions}
        )
    except KeyError as e:
        msg = f"Failed to create cluster_data DataFrame: missing column {e}"
        raise RuntimeError(msg) from None
    except ValueError as e:
        msg = f"Failed to create cluster_data DataFrame due to invalid value: {e}"
        raise RuntimeError(msg) from None
    except TypeError as e:
        msg = f"Failed to create cluster_data DataFrame due to incorrect type: {e}"
        raise RuntimeError(msg) from None

    # Store the results in adata.uns
    adata.uns["cluster_resolution_top_genes"] = top_genes_dict
    adata.uns["cluster_resolution_cluster_data"] = cluster_data


def _validate_cluster_resolution_inputs(
    adata: AnnData,
    resolutions: Sequence[float],
    method: str,
    flavor: str,
) -> None:
    """Validate inputs for the find_cluster_resolution function."""
    if not resolutions:
        msg = "resolutions list cannot be empty"
        raise ValueError(msg)
    if not all(isinstance(r, int | float) and r >= 0 for r in resolutions):
        msg = "All resolutions must be non-negative numbers"
        raise ValueError(msg)
    if method != "wilcoxon":
        msg = "Only method='wilcoxon' is supported"
        raise ValueError(msg)
    if flavor != "igraph":
        msg = "Only flavor='igraph' is supported"
        raise ValueError(msg)
    if "neighbors" not in adata.uns:
        msg = "adata must have precomputed neighbors (run sc.pp.neighbors first)."
        raise ValueError(msg)

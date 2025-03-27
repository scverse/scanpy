# tests/test_cluster_resolution.py
from __future__ import annotations

import re

import pandas as pd
import pytest

import scanpy as sc
from scanpy.tools._cluster_resolution import cluster_resolution_finder


# Test 1: Basic functionality
def test_cluster_resolution_finder_basic(adata):
    """Test that cluster_resolution_finder runs without errors and returns expected output."""
    resolutions = [0.1, 0.5]
    top_genes_dict, cluster_data = cluster_resolution_finder(
        adata,
        resolutions,
        prefix="leiden_res_",
        method="wilcoxon",
        n_top_genes=2,
        min_cells=2,
        deg_mode="within_parent",
        flavor="igraph",
        n_iterations=2,
        copy=False,
    )

    # Check output types
    assert isinstance(top_genes_dict, dict)
    assert isinstance(cluster_data, pd.DataFrame)

    # Check that clustering columns were added to adata.obs
    for res in resolutions:
        assert f"leiden_res_{res}" in adata.obs

    # Check that top_genes_dict has entries for some parent-child pairs
    assert len(top_genes_dict) > 0
    for (parent, child), genes in top_genes_dict.items():
        assert isinstance(parent, str)
        assert isinstance(child, str)
        assert isinstance(genes, list)
        assert len(genes) <= 2  # n_top_genes=2

    # Check that cluster_data has the expected columns
    for res in resolutions:
        assert f"leiden_res_{res}" in cluster_data.columns


# Test 2: Conflicting arguments (invalid deg_mode)
def test_cluster_resolution_finder_invalid_deg_mode(adata):
    """Test that an invalid deg_mode raises a ValueError."""
    with pytest.raises(
        ValueError, match=r"deg_mode must be 'within_parent' or 'per_resolution'"
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
            deg_mode="invalid_mode",
        )


# Test 3: Input values that should cause an error (empty resolutions)
def test_cluster_resolution_finder_empty_resolutions(adata):
    """Test that an empty resolutions list raises a ValueError."""
    with pytest.raises(ValueError, match=r"resolutions list cannot be empty"):
        cluster_resolution_finder(
            adata,
            resolutions=[],
        )


# Test 4: Input values that should cause an error (negative resolutions)
def test_cluster_resolution_finder_negative_resolutions(adata):
    """Test that negative resolutions raise a ValueError."""
    with pytest.raises(
        ValueError, match="All resolutions must be non-negative numbers"
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1, -0.5],
        )


# Test 5: Input values that should cause an error (missing neighbors)
def test_cluster_resolution_finder_missing_neighbors():
    """Test that an adata object without neighbors raises a ValueError."""
    adata = sc.datasets.pbmc68k_reduced()  # Create a fresh adata
    # Remove neighbors if they exist
    if "neighbors" in adata.uns:
        del adata.uns["neighbors"]
    # Also remove connectivities and distances to ensure leiden doesn't recompute
    if "connectivities" in adata.obsp:
        del adata.obsp["connectivities"]
    if "distances" in adata.obsp:
        del adata.obsp["distances"]
    with pytest.raises(
        ValueError,
        match=re.escape(
            "adata must have precomputed neighbors (run sc.pp.neighbors first)."
        ),
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
        )


# Test 6: Helpful error message (unsupported method)
def test_cluster_resolution_finder_unsupported_method(adata):
    """Test that an unsupported method raises a ValueError with a helpful message."""
    with pytest.raises(ValueError, match="Only method='wilcoxon' is supported"):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
            method="t-test",
        )


# Test 7: Bounds on returned values (n_top_genes)
@pytest.mark.parametrize("n_top_genes", [1, 3])
def test_cluster_resolution_finder_n_top_genes(adata, n_top_genes):
    """Test that n_top_genes bounds the number of genes returned."""
    top_genes_dict, _ = cluster_resolution_finder(
        adata,
        resolutions=[0.1, 0.5],
        n_top_genes=n_top_genes,
    )
    for genes in top_genes_dict.values():
        assert len(genes) <= n_top_genes


# Test 8: Orthogonal effects (copy argument)
def test_cluster_resolution_finder_copy_argument(adata):
    """Test that the copy argument doesn't affect the output but protects the input."""
    adata_original = adata.copy()

    # Run with copy=True
    top_genes_dict_copy, cluster_data_copy = cluster_resolution_finder(
        adata,
        resolutions=[0.1],
        copy=True,
    )

    # Check that adata wasn't modified
    assert adata.obs.equals(adata_original.obs)

    # Run with copy=False
    top_genes_dict_nocopy, cluster_data_nocopy = cluster_resolution_finder(
        adata,
        resolutions=[0.1],
        copy=False,
    )

    # Check that adata was modified
    assert "leiden_res_0.1" in adata.obs

    # Check that outputs are the same regardless of copy
    assert top_genes_dict_copy == top_genes_dict_nocopy
    assert cluster_data_copy.equals(cluster_data_nocopy)

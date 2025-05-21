from __future__ import annotations

import networkx as nx
import pytest

from scanpy.plotting._cluster_tree import cluster_decision_tree
from scanpy.tools._cluster_resolution import find_cluster_resolution
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.leidenalg]


@pytest.fixture
def adata_for_test():
    """Fixture to provide a preprocessed AnnData object for testing."""
    import scanpy as sc

    adata = pbmc68k_reduced()
    sc.pp.neighbors(adata)
    return adata


@pytest.fixture
def adata_with_clusters(adata_for_test):
    """Fixture providing clustering data and top_genes_dict for cluster_decision_tree."""
    adata = adata_for_test.copy()
    resolutions = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]
    find_cluster_resolution(
        adata,
        resolutions,
        prefix="leiden_res_",
        n_top_genes=2,
        min_cells=2,
        deg_mode="within_parent",
        flavor="igraph",
        n_iterations=2,
    )
    return adata, resolutions


# Test 1: Basic functionality without gene labels
def test_cluster_decision_tree_basic(adata_with_clusters):
    """Test that cluster_decision_tree runs without errors and returns a graph."""
    adata, resolutions = adata_with_clusters

    G = cluster_decision_tree(
        adata=adata,
        resolutions=resolutions,
    )

    assert isinstance(G, nx.DiGraph)
    assert len(G.nodes) > 0
    assert len(G.edges) > 0

    for node in G.nodes:
        assert "resolution" in G.nodes[node]
        assert "cluster" in G.nodes[node]


# Test 2: Basic functionality with gene labels
def test_cluster_decision_tree_with_gene_labels(adata_with_clusters):
    """Test that cluster_decision_tree handles gene labels when show_gene_labels is True."""
    adata, resolutions = adata_with_clusters

    G = cluster_decision_tree(
        adata=adata,
        resolutions=resolutions,
        gene_label_settings={
            "show_gene_labels": True,
            "n_top_genes": 2,
        },
    )

    assert isinstance(G, nx.DiGraph)
    assert len(G.nodes) > 0
    assert len(G.edges) > 0


# Test 3: Error condition (show_gene_labels=True but top_genes_dict missing in adata.uns)
def test_cluster_decision_tree_missing_top_genes_dict(adata_with_clusters):
    """Test that show_gene_labels=True raises an error if top_genes_dict is missing in adata.uns."""
    adata, resolutions = adata_with_clusters

    del adata.uns["cluster_resolution_top_genes"]

    with pytest.raises(
        ValueError,
        match=r"Gene labels requested but `adata\.uns\['cluster_resolution_top_genes'\]` not found\. Run `sc\.tl\.cluster_resolution_finder` first\.",
    ):
        cluster_decision_tree(
            adata=adata,
            resolutions=resolutions,
            gene_label_settings={"show_gene_labels": True},
        )


# Test 4: Conflicting arguments (negative node_size)
def test_cluster_decision_tree_negative_node_size(adata_with_clusters):
    """Test that a negative node_size raises a ValueError."""
    adata, resolutions = adata_with_clusters

    with pytest.raises(ValueError, match=r"node_size must be a positive number."):
        cluster_decision_tree(
            adata=adata, resolutions=resolutions, node_style={"node_size": -100}
        )


# Test 5: Error conditions (invalid figsize)
def test_cluster_decision_tree_invalid_figsize(adata_with_clusters):
    """Test that an invalid figsize raises a ValueError."""
    adata, resolutions = adata_with_clusters

    with pytest.raises(
        ValueError,
        match=r"figsize must be a tuple of two positive numbers \(width, height\)\.",
    ):
        cluster_decision_tree(
            adata=adata,
            resolutions=resolutions,
            output_settings={"figsize": (0, 5)},  # Invalid width
        )


# Test 6: Helpful error message (missing cluster_data in adata.uns)
def test_cluster_decision_tree_missing_cluster_data(adata_with_clusters):
    """Test that a missing cluster_data in adata.uns raises a ValueError."""
    adata, resolutions = adata_with_clusters

    del adata.uns["cluster_resolution_cluster_data"]

    with pytest.raises(
        ValueError,
        match=r"adata\.uns\['cluster_resolution_cluster_data'\] not found\. Run `sc\.tl\.cluster_resolution_finder` first\.",
    ):
        cluster_decision_tree(
            adata=adata,
            resolutions=resolutions,
        )


# Test 7: Orthogonal effects (draw argument)
def test_cluster_decision_tree_draw_argument(adata_with_clusters):
    """Test that the draw argument doesn't affect the graph output."""
    adata, resolutions = adata_with_clusters

    G_no_draw = cluster_decision_tree(
        adata=adata,
        resolutions=resolutions,
    )

    from unittest import mock

    with mock.patch("matplotlib.pyplot.show"):
        G_draw = cluster_decision_tree(adata=adata, resolutions=resolutions)

    assert nx.is_isomorphic(G_no_draw, G_draw)
    assert G_no_draw.nodes(data=True) == G_draw.nodes(data=True)

    def make_edge_hashable(edges):
        return {
            (
                u,
                v,
                tuple(
                    (k, tuple(v) if isinstance(v, list) else v)
                    for k, v in sorted(d.items())
                ),
            )
            for u, v, d in edges
        }

    assert make_edge_hashable(G_no_draw.edges(data=True)) == make_edge_hashable(
        G_draw.edges(data=True)
    )


# Test 8: Equivalent inputs (node_colormap)
@pytest.mark.parametrize(
    "node_colormap",
    [
        None,
        ["Set3", "Set3"],
    ],
)
def test_cluster_decision_tree_node_colormap(adata_with_clusters, node_colormap):
    """Test that node_colormap=None and a uniform colormap produce similar results."""
    adata, resolutions = adata_with_clusters

    G = cluster_decision_tree(
        adata=adata,
        resolutions=resolutions,
        node_style={"node_colormap": node_colormap},
    )
    assert isinstance(G, nx.DiGraph)
    assert len(G.nodes) > 0


# Test 9: Bounds on gene labels (n_top_genes)
@pytest.mark.parametrize("n_top_genes", [1, 3])
def test_cluster_decision_tree_n_top_genes(adata_with_clusters, n_top_genes):
    """Test that n_top_genes parameter works correctly."""
    adata, resolutions = adata_with_clusters

    G = cluster_decision_tree(
        adata=adata,
        resolutions=resolutions,
        gene_label_settings={"show_gene_labels": True, "n_top_genes": n_top_genes},
    )

    assert isinstance(G, nx.DiGraph)
    assert len(G.nodes) > 0
    assert len(G.edges) > 0

    for node in G.nodes:
        if "top_genes" in G.nodes[node]:
            assert len(G.nodes[node]["top_genes"]) == n_top_genes

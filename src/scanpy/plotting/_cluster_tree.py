from __future__ import annotations

import math
from typing import TYPE_CHECKING

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path

if TYPE_CHECKING:
    from typing import Literal

    import networkx as nx
    import pandas as pd
    from anndata import AnnData
    from pandas import DataFrame


def count_crossings(
    G: nx.DiGraph,
    pos: dict[str, tuple[float, float]],
    edges: list[tuple[str, str]],
    level_nodes: dict[int, list[str]],
) -> int:
    """Count the number of edge crossings in the graph based on node positions.

    Args:
        G: Directed graph with nodes and edges.
        pos: Dictionary mapping nodes to their (x, y) positions.
        edges: List of edge tuples (u, v).
        level_nodes: Dictionary mapping resolution levels to lists of nodes.

    Returns
    -------
        Number of edge crossings.
    """
    crossings = 0
    for i, (u1, v1) in enumerate(edges):
        for j, (u2, v2) in enumerate(edges[i + 1 :], start=i + 1):
            # Skip edges at the same level to avoid counting self-crossings
            level_u1 = G.nodes[u1]["resolution"]
            level_v1 = G.nodes[v1]["resolution"]
            level_u2 = G.nodes[u2]["resolution"]
            level_v2 = G.nodes[v2]["resolution"]
            if level_u1 == level_u2 and level_v1 == level_v2:
                continue

            # Get coordinates of the edge endpoints
            x1_start, y1_start = pos[u1]
            x1_end, y1_end = pos[v1]
            x2_start, y2_start = pos[u2]
            x2_end, y2_end = pos[v2]

            # Compute the direction vectors of the edges
            dx1 = x1_end - x1_start
            dy1 = y1_end - y1_start
            dx2 = x2_end - x2_start
            dy2 = y2_end - y2_start

            # Compute the denominator for the line intersection formula
            denom = dx1 * dy2 - dy1 * dx2
            if abs(denom) < 1e-8:  # Adjusted threshold for numerical stability
                continue

            # Compute intersection parameters s and t
            s = ((x2_start - x1_start) * dy2 - (y2_start - y1_start) * dx2) / denom
            t = ((x2_start - x1_start) * dy1 - (y2_start - y1_start) * dx1) / denom

            # Check if the intersection occurs within both edge segments
            if 0 < s < 1 and 0 < t < 1:
                crossings += 1

    return crossings


def optimize_node_ordering(
    G: nx.DiGraph,
    pos: dict[str, tuple[float, float]],
    edges: list[tuple[str, str]],
    resolutions: list[str],
    max_iterations: int = 10,
) -> None:
    """Optimize node ordering at each level to minimize edge crossings by swapping adjacent nodes.

    Args:
        G: Directed graph with nodes and edges.
        pos: Dictionary mapping nodes to their (x, y) positions.
        edges: List of edge tuples (u, v).
        resolutions: List of resolution identifiers.
        max_iterations: Maximum number of iterations per level to prevent excessive computation.
    """
    # Group nodes by resolution level
    level_nodes = {
        res_idx: [node for node in G.nodes if G.nodes[node]["resolution"] == res_idx]
        for res_idx in range(len(resolutions))
    }

    for res_idx in range(len(resolutions)):
        nodes = level_nodes[res_idx]
        if len(nodes) < 2:
            continue

        # Sort nodes by their x-coordinate to establish an initial order
        nodes.sort(key=lambda node: pos[node][0])

        iteration = 0
        improved = True
        while improved and iteration < max_iterations:
            improved = False
            for i in range(len(nodes) - 1):
                node1, node2 = nodes[i], nodes[i + 1]
                x1, y1 = pos[node1]
                x2, y2 = pos[node2]

                # Compute current number of crossings
                current_crossings = count_crossings(G, pos, edges, level_nodes)

                # Swap positions and compute new crossings
                pos[node1] = (x2, y1)
                pos[node2] = (x1, y2)
                new_crossings = count_crossings(G, pos, edges, level_nodes)

                # If swapping reduces crossings, keep the swap
                if new_crossings < current_crossings:
                    nodes[i], nodes[i + 1] = nodes[i + 1], nodes[i]
                    improved = True
                else:
                    # Revert the swap if it doesn't improve crossings
                    pos[node1] = (x1, y1)
                    pos[node2] = (x2, y2)

            iteration += 1


def evaluate_bezier(
    t: float, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
) -> np.ndarray:
    """Evaluate a cubic Bezier curve at parameter t.

    Args:
        t: Parameter value in [0, 1] where the curve is evaluated.
        p0: Starting point of the Bezier curve.
        p1: First control point.
        p2: Second control point.
        p3: Ending point of the Bezier curve.

    Returns
    -------
        The (x, y) coordinates on the Bezier curve at parameter t.

    Raises
    ------
        ValueError: If t is not in [0, 1].
    """
    if not 0 <= t <= 1:
        msg = "Parameter t must be in the range [0, 1]"
        raise ValueError(msg)

    t2 = t * t
    t3 = t2 * t
    mt = 1 - t
    mt2 = mt * mt
    mt3 = mt2 * mt
    return mt3 * p0 + 3 * mt2 * t * p1 + 3 * mt * t2 * p2 + t3 * p3


def evaluate_bezier_tangent(
    t: float, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
) -> np.ndarray:
    """Compute the tangent vector of a cubic Bezier curve at parameter t.

    Args:
        t: Parameter value in [0, 1] where the tangent is computed.
        p0: Starting point of the Bezier curve.
        p1: First control point.
        p2: Second control point.
        p3: Ending point of the Bezier curve.

    Returns
    -------
        The tangent vector (dx/dt, dy/dt) at parameter t.

    Raises
    ------
        ValueError: If t is not in [0, 1].
    """
    if not 0 <= t <= 1:
        msg = "Parameter t must be in the range [0, 1]"
        raise ValueError(msg)

    t2 = t * t
    mt = 1 - t
    mt2 = mt * mt
    return 3 * mt2 * (p1 - p0) + 6 * mt * t * (p2 - p1) + 3 * t2 * (p3 - p2)


def build_cluster_graph(
    data: DataFrame, prefix: str = "leiden_res_", edge_threshold: float = 0.02
) -> nx.DiGraph:
    """Build a directed graph representing hierarchical clustering from data.

    Args:
        data: DataFrame containing clustering results with columns named as '{prefix}{resolution}'.
        prefix: Prefix for column names (default: "leiden_res_").
        edge_threshold: Minimum fraction of samples to create an edge between clusters.

    Returns
    -------
        graph G: Directed graph representing hierarchical clustering.

    Raises
    ------
        ValueError: If no columns in the DataFrame match the given prefix.
    """
    import networkx as nx

    # Validate input data
    matching_columns = [col for col in data.columns if col.startswith(prefix)]
    if not matching_columns:
        msg = f"No columns found with prefix '{prefix}' in the DataFrame."
        raise ValueError(msg)

    G = nx.DiGraph()

    # Extract resolutions from column names
    resolutions = [col[len(prefix) :] for col in matching_columns]
    resolutions.sort()

    # Add nodes with resolution attribute for layout
    for i, res in enumerate(resolutions):
        clusters = data[f"{prefix}{res}"].unique()
        for cluster in sorted(clusters):
            node = f"{res}_C{cluster}"
            G.add_node(node, resolution=i, cluster=cluster)

    # Build edges between consecutive resolutions
    for i in range(len(resolutions) - 1):
        res1 = f"{prefix}{resolutions[i]}"
        res2 = f"{prefix}{resolutions[i + 1]}"

        grouped = (
            data.loc[:, [res1, res2]]
            .astype(str)
            .groupby(res1, observed=False)[res2]
            .value_counts(normalize=True)
        )

        for key, frac in grouped.items():
            parent, child = key if isinstance(key, tuple) else (key, None)
            parent = str(parent) if parent is not None else ""
            child = str(child)
            parent_node = f"{resolutions[i]}_C{parent}"
            child_node = f"{resolutions[i + 1]}_C{child}"
            G.add_edge(parent_node, child_node, weight=frac)

    return G


def compute_cluster_layout(
    G: nx.DiGraph,
    node_spacing: float = 10.0,
    level_spacing: float = 1.5,
    orientation: str = "vertical",
    barycenter_sweeps: int = 2,
    *,
    use_reingold_tilford: bool = False,
) -> dict[str, tuple[float, float]]:
    """Compute node positions for the cluster decision tree with crossing minimization.

    Args:
        G: Directed graph with nodes and edges.
        node_spacing: Horizontal spacing between nodes at the same level.
        level_spacing: Vertical spacing between resolution levels.
        orientation: Orientation of the tree ("vertical" or "horizontal").
        barycenter_sweeps: Number of barycenter-based reordering sweeps.
        use_reingold_tilford: Whether to use the Reingold-Tilford layout (requires igraph).

    Returns
    -------
        Dictionary mapping nodes to their (x, y) positions.
    """
    import networkx as nx

    # Step 1: Calculate initial node positions
    if use_reingold_tilford:
        try:
            import igraph as ig

            nodes = list(G.nodes)
            edges = [(u, v) for u, v in G.edges()]
            g = ig.Graph()
            g.add_vertices(nodes)
            g.add_edges([(nodes.index(u), nodes.index(v)) for u, v in edges])
            layout = g.layout_reingold_tilford(root=[0])
            pos = {node: coord for node, coord in zip(nodes, layout.coords)}
        except ImportError as e:
            print(
                f"igraph not installed or failed: {e}. Falling back to multipartite_layout."
            )
            pos = nx.multipartite_layout(
                G, subset_key="resolution", scale=int(node_spacing)
            )
        except Exception as e:
            print(
                f"Error in Reingold-Tilford layout: {e}. Falling back to multipartite_layout."
            )
            pos = nx.multipartite_layout(
                G, subset_key="resolution", scale=int(node_spacing)
            )
    else:
        pos = nx.multipartite_layout(
            G, subset_key="resolution", scale=int(node_spacing)
        )

    # Step 2: Adjust orientation (vertical: lower resolutions at top, higher at bottom)
    if orientation == "vertical":
        pos = {node: (y, -x) for node, (x, y) in pos.items()}

    # Step 3: Increase vertical spacing between levels
    new_pos = {}
    for node, (x, y) in pos.items():
        new_y = y * level_spacing
        new_pos[node] = (x, new_y)
    pos = new_pos

    # Step 4: Barycenter-based reordering to minimize edge crossings
    resolutions = sorted(set(node.split("_")[0] for node in G.nodes))
    for sweep in range(barycenter_sweeps):
        # Downward sweep: Adjust nodes based on parent positions
        for i in range(1, len(resolutions)):
            res = resolutions[i]
            nodes_at_level = [node for node in G.nodes if node.startswith(f"{res}_C")]
            node_to_barycenter = {}
            for node in nodes_at_level:
                parents = list(G.predecessors(node))
                barycenter = (
                    np.mean([pos[parent][0] for parent in parents]) if parents else 0
                )
                cluster_id = int(node.split("_C")[1])
                node_to_barycenter[node] = (barycenter, cluster_id)
            sorted_nodes = sorted(
                node_to_barycenter.keys(), key=lambda x: node_to_barycenter[x]
            )
            y_level = pos[sorted_nodes[0]][1]
            n_nodes = len(sorted_nodes)
            x_positions = (
                np.linspace(
                    -node_spacing * (n_nodes - 1) / 2,
                    node_spacing * (n_nodes - 1) / 2,
                    n_nodes,
                )
                if n_nodes > 1
                else [0]
            )
            for node, x in zip(sorted_nodes, x_positions):
                pos[node] = (x, y_level)

        # Upward sweep: Adjust nodes based on child positions
        for i in range(len(resolutions) - 2, -1, -1):
            res = resolutions[i]
            nodes_at_level = [node for node in G.nodes if node.startswith(f"{res}_C")]
            node_to_barycenter = {}
            for node in nodes_at_level:
                children = list(G.successors(node))
                barycenter = (
                    np.mean([pos[child][0] for child in children]) if children else 0
                )
                cluster_id = int(node.split("_C")[1])
                node_to_barycenter[node] = (barycenter, cluster_id)
            sorted_nodes = sorted(
                node_to_barycenter.keys(), key=lambda x: node_to_barycenter[x]
            )
            y_level = pos[sorted_nodes[0]][1]
            n_nodes = len(sorted_nodes)
            x_positions = (
                np.linspace(
                    -node_spacing * (n_nodes - 1) / 2,
                    node_spacing * (n_nodes - 1) / 2,
                    n_nodes,
                )
                if n_nodes > 1
                else [0]
            )
            for node, x in zip(sorted_nodes, x_positions):
                pos[node] = (x, y_level)

    # Step 5: Optimize node ordering to further reduce crossings
    filtered_edges = [
        (u, v, d["weight"]) for u, v, d in G.edges(data=True) if d["weight"] >= 0.02
    ]
    edges = [(u, v) for u, v, w in filtered_edges]
    edges_set = set(edges)
    if len(edges_set) < len(edges):
        print(
            f"Warning: Found {len(edges) - len(edges_set)} duplicate edges in the visualization."
        )
    edges = list(edges_set)
    optimize_node_ordering(G, pos, edges, resolutions)

    return pos


def draw_curved_edge(
    ax,
    start_x: float,
    start_y: float,
    end_x: float,
    end_y: float,
    *,
    linewidth: float,
    color: str,
    edge_curvature: float = 0.1,
    arrow_size: float = 12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Draw a gentle S-shaped curved edge between two points with an arrowhead.

    Args:
        ax: Matplotlib axis to draw on.
        start_x, start_y: Starting coordinates of the edge.
        end_x, end_y: Ending coordinates of the edge.
        linewidth: Width of the edge.
        color: Color of the edge.
        edge_curvature: Controls the intensity of the S-shape (smaller values for subtler curves).
        arrow_size: Size of the arrowhead.

    Returns
    -------
        Tuple of Bézier control points (p0, p1, p2, p3) for label positioning.
    """
    # Define the start and end points
    p0 = np.array([start_x, start_y])
    p3 = np.array([end_x, end_y])

    # Calculate the vector from start to end
    vec = p3 - p0
    length = np.sqrt(vec[0] ** 2 + vec[1] ** 2)

    if length == 0:
        empty_array = np.array([[], []])
        return empty_array, empty_array, empty_array, empty_array

    # Unit vector along the edge
    unit_vec = vec / length
    # Perpendicular vector for creating the S-shape
    perp_vec = np.array([-unit_vec[1], unit_vec[0]])

    # Define control points for a single cubic Bézier curve with an S-shape
    # Place control points at 1/3 and 2/3 along the edge, with small perpendicular offsets
    offset = length * edge_curvature
    p1 = p0 + (p3 - p0) / 3 + perp_vec * offset  # First control point (bend outward)
    p2 = (
        p0 + 2 * (p3 - p0) / 3 - perp_vec * offset
    )  # Second control point (bend inward)

    # Define the path vertices and codes for a single cubic Bézier curve
    vertices = [
        (start_x, start_y),  # Start point
        (p1[0], p1[1]),  # First control point
        (p2[0], p2[1]),  # Second control point
        (end_x, end_y),  # End point
    ]
    codes = [
        Path.MOVETO,  # Move to start
        Path.CURVE4,  # Cubic Bézier curve (needs 3 points: p0, p1, p2)
        Path.CURVE4,  # Continuation of the Bézier curve
        Path.CURVE4,  # End of the Bézier curve
    ]

    # Create the path
    path = Path(vertices, codes)

    # Draw the curve
    patch = PathPatch(
        path, facecolor="none", edgecolor=color, linewidth=linewidth, alpha=0.8
    )
    ax.add_patch(patch)

    # Add an arrowhead at the end
    # t = 0.9  # Near the end of the curve
    # tangent = evaluate_bezier_tangent(t, p0, p1, p2, p3)
    # tangent_angle = np.arctan2(tangent[1], tangent[0])
    arrow = FancyArrowPatch(
        (end_x, end_y),
        # (end_x - 0.01 * np.cos(tangent_angle), end_y - 0.01 * np.sin(tangent_angle)),
        (end_x, end_y),
        arrowstyle="->",
        mutation_scale=arrow_size,
        color=color,
        linewidth=linewidth,
        alpha=0.8,
    )
    ax.add_patch(arrow)

    return p0, p1, p2, p3


def draw_gene_labels(
    ax,
    pos: dict[str, tuple[float, float]],
    gene_labels: dict[str, str],
    *,
    node_sizes: dict[str, float],
    node_colors: dict[str, str],
    offset: float = 0.2,
    fontsize: float = 8,
) -> dict[str, float]:
    """Draw gene labels in boxes below nodes with matching boundary colors.

    Args:
        ax: Matplotlib axis to draw on.
        pos: Dictionary mapping nodes to their (x, y) positions.
        gene_labels: Dictionary mapping nodes to their gene labels.
        node_sizes: Dictionary mapping nodes to their sizes.
        node_colors: Dictionary mapping nodes to their colors.
        offset: Distance below the node to place the label (in data coordinates).

    Returns
    -------
        Dictionary mapping nodes to the bottom y-coordinate of their label boxes.
    """
    gene_label_bottoms = {}
    for node, label in gene_labels.items():
        if label:
            x, y = pos[node]
            # Compute the node radius in data coordinates
            radius = math.sqrt(node_sizes[node] / math.pi)
            fig_width, fig_height = ax.figure.get_size_inches()
            radius_fig = radius / (72 * fig_height)
            # xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            data_height = ylim[0] - ylim[1]
            radius_data = radius_fig * data_height

            # Position the top of the label box just below the node
            box_top_y = y - radius_data - offset

            # Compute the height of the label box based on the number of lines
            num_lines = label.count("\n") + 1
            line_height = 0.03  # Reduced line height for better scaling
            label_height = num_lines * line_height + 0.04  # Reduced padding
            box_center_y = box_top_y - label_height / 2

            # Draw the label
            ax.text(
                x,
                box_center_y,
                label,
                fontsize=fontsize,
                ha="center",
                va="center",
                color="black",
                bbox=dict(
                    facecolor="white",
                    edgecolor=node_colors[node],
                    boxstyle="round,pad=0.2",  # Reduced padding for the box
                ),
            )
            gene_label_bottoms[node] = box_top_y - label_height
    return gene_label_bottoms


def draw_cluster_tree(
    # Core Inputs
    G: nx.DiGraph,
    pos: dict[str, tuple[float, float]],
    data: pd.DataFrame,
    prefix: str,
    resolutions: list[float],
    *,
    # Output and Display Options
    output_path: str | None = None,
    draw: bool = True,
    figsize: tuple[float, float] = (10, 8),
    dpi: float = 300,
    # Node Appearance
    node_size: float = 500,
    node_color: str = "prefix",
    node_colormap: list[str] | None = None,
    node_label_fontsize: float = 12,
    # Edge Appearance
    edge_color: str = "parent",
    edge_curvature: float = 0.01,
    edge_threshold: float = 0.05,
    show_weight: bool = True,
    edge_label_threshold: float = 0.1,
    edge_label_position: float = 0.5,
    edge_label_fontsize: float = 8,
    # Gene Label Options
    top_genes_dict: dict[tuple[str, str], list[str]] | None = None,
    show_gene_labels: bool = False,
    n_top_genes: int = 2,
    gene_label_offset: float = 0.3,
    gene_label_fontsize: float = 10,
    gene_label_threshold: float = 0.05,
    # Level Label Options
    level_label_offset: float = 5,
    level_label_fontsize: float = 12,
    # Title Options
    title: str = "Hierarchical Leiden Clustering",
    title_fontsize: float = 16,
) -> None:
    """
    Draw a hierarchical clustering decision tree with nodes, edges, and optional gene labels.

    This function visualizes a hierarchical clustering tree where nodes represent clusters at different
    resolutions, edges represent transitions between clusters, and edge weights indicate the proportion
    of cells transitioning from a parent cluster to a child cluster. The tree can include gene labels
    showing differentially expressed genes (DEGs) between parent and child clusters.

    Args:
        G (nx.DiGraph):
            Directed graph representing the clustering hierarchy. Nodes should have a 'resolution'
            attribute, and edges should have a 'weight' attribute indicating the proportion of cells
            transitioning from the parent to the child cluster.
        pos (Dict[str, Tuple[float, float]]):
            Dictionary mapping node names (e.g., "res_0.0_C0") to their (x, y) positions in the plot.
        data (pd.DataFrame):
            DataFrame containing clustering results, with columns named as '{prefix}{resolution}'
            (e.g., 'leiden_res_0.0', 'leiden_res_0.5') indicating cluster assignments for each cell.
        prefix (str):
            Prefix for column names in the DataFrame (e.g., "leiden_res_"). Used to identify clustering
            columns and label resolution levels in the plot.
        resolutions (List[float]):
            List of resolution values to include in the visualization (e.g., [0.0, 0.5, 1.0]). Determines
            the levels of the tree, with each resolution corresponding to a level from top to bottom.

        output_path (Optional[str], optional):
            Path to save the figure (e.g., 'cluster_tree.png'). Supports formats like PNG, PDF, SVG.
            If None, the figure is not saved. Defaults to None.
        draw (bool, optional):
            Whether to display the plot using plt.show(). If False, the plot is created but not displayed.
            Defaults to True.
        figsize (Tuple[float, float], optional):
            Figure size as (width, height) in inches. Controls the overall size of the plot.
            Defaults to (10, 8).
        dpi (float, optional):
            Resolution for saving the figure (dots per inch). Higher values result in higher-quality output.
            Defaults to 300.

        node_size (float, optional):
            Base size for nodes in points^2 (area of the node). Node sizes are scaled within each level
            based on cluster sizes, using this value as the maximum size. Defaults to 500.
        node_color (str, optional):
            Color specification for nodes. If "prefix", nodes are colored by resolution level using a
            distinct color palette for each level. Alternatively, a single color can be specified
            (e.g., "red", "#FF0000"). Defaults to "prefix".
        node_colormap (Optional[List[str]], optional):
            Custom colormap for nodes, as a list of colors or colormaps (one per resolution level).
            Each entry can be a color (e.g., "red", "#FF0000") or a colormap name (e.g., "viridis").
            If None, the default "Set3" palette is used for "prefix" coloring. Defaults to None.
        node_label_fontsize (float, optional):
            Font size for node labels (e.g., cluster numbers like "0", "1"). Defaults to 12.

        edge_color (str, optional):
            Color specification for edges. Options are:
            - "parent": Edges inherit the color of the parent node.
            - "samples": Edges are colored by weight using the "viridis" colormap.
            - A single color (e.g., "blue", "#0000FF").
            Defaults to "parent".
        edge_curvature (float, optional):
            Curvature of edges, controlling the intensity of the S-shape. Smaller values result in subtler
            curves, while larger values create more pronounced S-shapes. Defaults to 0.1.
        edge_threshold (float, optional):
            Minimum weight (proportion of cells) required to draw an edge. Edges with weights below this
            threshold are not drawn, reducing clutter. Defaults to 0.5.
        show_weight (bool, optional):
            Whether to show edge weights as labels on the edges. If True, weights above `edge_label_threshold`
            are displayed. Defaults to True.
        edge_label_threshold (float, optional):
            Minimum weight required to label an edge with its weight. Only edges with weights above this
            threshold will have labels (if `show_weight` is True). Defaults to 0.7.
        edge_label_position (float, optional):
            Position of the edge weight label along the edge, as a ratio from 0.0 (near the parent node) to
            1.0 (near the child node). A value of 0.5 places the label at the midpoint. A small buffer is
            applied to avoid overlap with nodes. Defaults to 0.5.
        edge_label_fontsize (float, optional):
            Font size for edge weight labels (e.g., "0.86"). Defaults to 8.

        top_genes_dict (Optional[Dict[Tuple[str, str], List[str]]], optional):
            Dictionary mapping (parent, child) node pairs to lists of differentially expressed genes (DEGs).
            Keys are tuples of node names (e.g., ("res_0.0_C0", "res_0.5_C1")), and values are lists of gene
            names (e.g., ["GeneA", "GeneB"]). If provided and `show_gene_labels` is True, DEGs are displayed
            below child nodes. Defaults to None.
        show_gene_labels (bool, optional):
            Whether to show gene labels (DEGs) below child nodes. Requires `top_genes_dict` to be provided.
            Defaults to False.
        n_top_genes (int, optional):
            Number of top genes to display for each (parent, child) pair. Genes are taken from `top_genes_dict`
            in the order provided. Defaults to 2.
        gene_label_offset (float, optional):
            Vertical offset (in data coordinates) for gene labels below nodes. Controls the distance between
            the node and its gene label. Defaults to 0.2.
        gene_label_fontsize (float, optional):
            Font size for gene labels (e.g., gene names like "GeneA"). Defaults to 10.
        gene_label_threshold (float, optional):
            Minimum weight (proportion of cells) required to display a gene label for a (parent, child) pair.
            Gene labels are only shown for edges with weights above this threshold. Defaults to 0.05.

        level_label_offset (float, optional):
            Horizontal buffer space (in data coordinates) between the level labels (e.g., "leiden_res_0.0")
            and the leftmost node at the bottom level. Controls the spacing of level labels on the left side
            of the plot. Defaults to 0.5.
        level_label_fontsize (float, optional):
            Font size for level labels (e.g., "leiden_res_0.0"). Defaults to 12.

        title (str, optional):
            Title of the plot, displayed at the top. Defaults to "Hierarchical Leiden Clustering".
        title_fontsize (float, optional):
            Font size for the plot title. Defaults to 16.
    """
    import networkx as nx
    import seaborn as sns

    # Step 1: Compute cluster sizes
    cluster_sizes = {}
    for res in resolutions:
        res_key = f"{prefix}{res}"
        counts = data[res_key].value_counts()
        for cluster, count in counts.items():
            node = f"{res}_C{cluster}"
            cluster_sizes[node] = count

    # Step 2: Scale node sizes within each level
    node_sizes = {}
    for i, res in enumerate(resolutions):
        nodes_at_level = [
            f"{res}_C{cluster}" for cluster in data[f"{prefix}{res}"].unique()
        ]
        sizes = np.array([cluster_sizes[node] for node in nodes_at_level])
        if len(sizes) > 1:
            min_size, max_size = sizes.min(), sizes.max()
            if min_size != max_size:
                normalized_sizes = 0.5 + (sizes - min_size) / (max_size - min_size)
            else:
                normalized_sizes = np.ones_like(sizes)
            scaled_sizes = normalized_sizes * node_size
        else:
            scaled_sizes = np.array([node_size])
        for node, scaled_size in zip(nodes_at_level, scaled_sizes):
            node_sizes[node] = scaled_size

    # Step 3: Generate color schemes for nodes
    if node_color == "prefix":
        if node_colormap is None:
            color_schemes = {
                r: sns.color_palette("Set3", n_colors=data[f"{prefix}{r}"].nunique())
                for r in resolutions
            }
        else:
            if len(node_colormap) < len(resolutions):
                print(
                    f"Warning: node_colormap has {len(node_colormap)} entries, but there are {len(resolutions)} resolutions. Cycling colors."
                )
                node_colormap = list(node_colormap) + [
                    node_colormap[i % len(node_colormap)]
                    for i in range(len(resolutions) - len(node_colormap))
                ]
            color_schemes = {}
            for i, r in enumerate(resolutions):
                color_spec = node_colormap[i]
                if (
                    isinstance(color_spec, str)
                    and mcolors.is_color_like(color_spec)
                    or isinstance(color_spec, tuple)
                    and len(color_spec) in (3, 4)
                    and all(isinstance(x, int | float) for x in color_spec)
                ):
                    color_schemes[r] = [color_spec]
                else:
                    try:
                        color_schemes[r] = sns.color_palette(
                            color_spec, n_colors=data[f"{prefix}{r}"].nunique()
                        )
                    except ValueError:
                        print(
                            f"Warning: '{color_spec}' is not a valid color or colormap for {r}. Using 'Set3'."
                        )
                        color_schemes[r] = sns.color_palette(
                            "Set3", n_colors=data[f"{prefix}{r}"].nunique()
                        )
    else:
        color_schemes = None

    # Step 4: Assign colors to nodes
    node_colors = {}

    for res in resolutions:
        clusters = data[f"{prefix}{res}"].unique()
        for cluster in clusters:
            node = f"{res}_C{cluster}"
            if node_color == "prefix":
                # Defensive check to satisfy linters/type checkers
                if color_schemes is None:
                    msg = "color_schemes is None when node_color is 'prefix', which should not happen."
                    raise RuntimeError(msg)
                if len(color_schemes[res]) == 1:
                    node_colors[node] = color_schemes[res][0]
                else:
                    node_colors[node] = color_schemes[res][
                        int(cluster) % len(color_schemes[res])
                    ]
            else:
                node_colors[node] = node_color

    # Step 5: Initialize the plot
    plt.figure(figsize=figsize, dpi=dpi)
    ax = plt.gca()

    # Step 6: Compute edge weights and colors
    edges = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] >= edge_threshold]
    weights = [
        max(d["weight"] * 5, 1.0)
        for u, v, d in G.edges(data=True)
        if d["weight"] >= edge_threshold
    ]
    edge_colors = []
    # for u, v in [(u, v) for u, v in G.edges()]:
    for u, v in edges:
        d = G[u][v]
        if edge_color == "parent":
            edge_colors.append(node_colors[u])
        elif edge_color == "samples":
            edge_colors.append(plt.cm.get_cmap("viridis")(d["weight"] / 5))
        else:
            edge_colors.append(edge_color)

    # Step 7: Draw nodes and node labels
    node_labels = {}
    gene_labels = {}
    for res in resolutions:
        clusters = data[f"{prefix}{res}"].unique()
        for cluster in clusters:
            node = f"{res}_C{cluster}"
            color = node_colors[node]
            size = node_sizes[node]
            nx.draw_networkx_nodes(
                G,
                pos,
                nodelist=[node],
                node_size=size,
                node_color=color,
                edgecolors="none",
            )
            node_labels[node] = str(cluster)
            if show_gene_labels and top_genes_dict:
                # Find the resolution of the parent level
                res_idx = resolutions.index(float(res))
                if res_idx == 0:
                    continue  # No parent level for the top resolution
                parent_res = resolutions[res_idx - 1]
                parent_clusters = data[f"{prefix}{parent_res}"].unique()
                for parent_cluster in parent_clusters:
                    parent_node = f"{parent_res}_C{parent_cluster}"
                    try:
                        edge_weight = G[parent_node][node]["weight"]
                    except KeyError:
                        continue
                    if edge_weight >= gene_label_threshold:
                        key = (f"res_{parent_node}", f"res_{node}")
                        if key in top_genes_dict:
                            genes = top_genes_dict[key][:n_top_genes]
                            gene_labels[node] = "\n".join(genes) if genes else ""

    nx.draw_networkx_labels(
        G,
        pos,
        labels=node_labels,
        font_size=int(node_label_fontsize),
        font_color="black",
    )

    # Step 8: Draw gene labels below nodes
    gene_label_bottoms = {}
    if show_gene_labels and gene_labels:
        gene_label_bottoms = draw_gene_labels(
            ax,
            pos,
            gene_labels,
            node_sizes=node_sizes,
            node_colors=node_colors,
            offset=gene_label_offset,
            fontsize=gene_label_fontsize,
        )

    # Step 9: Draw edges with labels using the new S-shaped edge function
    edge_labels = {
        (u, v): f"{w:.2f}"
        for u, v, w in [
            (u, v, d["weight"])
            for u, v, d in G.edges(data=True)
            if d["weight"] >= edge_threshold
        ]
        if w >= edge_label_threshold
    }

    # for (u, v), w, e_color in zip([(u, v) for u, v in G.edges()], weights, edge_colors):
    for (u, v), w, e_color in zip(edges, weights, edge_colors):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        radius_parent = math.sqrt(node_sizes[u] / math.pi)
        radius_child = math.sqrt(node_sizes[v] / math.pi)
        fig_width, fig_height = figsize
        radius_parent_fig = radius_parent / (72 * fig_height)
        radius_child_fig = radius_child / (72 * fig_height)
        # xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        data_height = ylim[0] - ylim[1]
        radius_parent_data = radius_parent_fig * data_height
        radius_child_data = radius_child_fig * data_height
        start_y = (
            gene_label_bottoms[u]
            if (show_gene_labels and u in gene_label_bottoms and gene_labels.get(u))
            else y1 - radius_parent_data
        )
        start_x = x1
        end_x, end_y = x2, y2 - radius_child_data

        # Draw the S-shaped edge
        p0, p1, p2, p3 = draw_curved_edge(
            ax,
            start_x,
            start_y,
            end_x,
            end_y,
            linewidth=w,
            color=e_color,
            edge_curvature=edge_curvature,
        )

        # Add edge label if required
        if show_weight and (u, v) in edge_labels and p0 is not None:
            t = edge_label_position
            point = evaluate_bezier(t, p0, p1, p2, p3)
            label_x, label_y = point[0], point[1]
            tangent = evaluate_bezier_tangent(t, p0, p1, p2, p3)
            tangent_angle = np.arctan2(tangent[1], tangent[0])
            rotation = np.degrees(tangent_angle)
            if rotation > 90:
                rotation -= 180
            elif rotation < -90:
                rotation += 180
            ax.text(
                label_x,
                label_y,
                edge_labels[(u, v)],
                fontsize=edge_label_fontsize,
                rotation=rotation,
                ha="center",
                va="center",
                bbox=None,
            )

    # Step 10: Draw level labels
    level_positions = {}
    for node, (x, y) in pos.items():
        res = node.split("_")[0]
        level_positions[res] = y

    # Count the number of clusters at each resolution
    cluster_counts = {}
    for res in resolutions:
        res_str = f"{res:.1f}"
        col_name = f"{prefix}{res_str}"
        if col_name not in data.columns:
            msg = f"Column {col_name} not found in data. Ensure clustering results are present."
            raise ValueError(msg)
        # Count unique clusters at this resolution
        num_clusters = len(data[col_name].dropna().unique())
        cluster_counts[res_str] = num_clusters

    # Draw the level labels
    min_x = min(p[0] for p in pos.values())
    label_offset = min_x - level_label_offset
    for i, res in enumerate(resolutions):
        res_str = f"{res:.1f}"
        label_pos = level_positions[res_str]
        num_clusters = cluster_counts[res_str]
        label_text = f"Resolution {res_str}:\n {num_clusters} clusters"
        plt.text(
            label_offset,
            label_pos,
            label_text,
            fontsize=level_label_fontsize,
            verticalalignment="center",
            bbox=dict(facecolor="white", edgecolor="black", alpha=0.7),
        )

    # Step 11: Finalize the plot
    plt.axis("off")
    plt.title(title, fontsize=title_fontsize)
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if draw:
        plt.show()
    plt.close()


def cluster_decision_tree(
    # Core Inputs
    adata: AnnData,
    prefix: str = "leiden_res_",
    resolutions: list[float] = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0],
    *,
    # Layout Options
    orientation: Literal["vertical", "horizontal"] = "vertical",
    node_spacing: float = 5.0,
    level_spacing: float = 1.5,
    barycenter_sweeps: int = 2,
    use_reingold_tilford: bool = False,
    # Output and Display Options
    output_path: str | None = None,
    draw: bool = True,
    figsize: tuple[float, float] = (10, 8),
    dpi: float = 300,
    # Node Appearance
    node_size: float = 500,
    node_color: str = "prefix",
    node_colormap: list[str] | None = None,
    node_label_fontsize: float = 12,
    # Edge Appearance
    edge_color: Literal["parent", "samples"] | str = "parent",
    edge_curvature: float = 0.01,
    edge_threshold: float = 0.05,
    show_weight: bool = True,
    edge_label_threshold: float = 0.1,
    edge_label_position: float = 0.5,
    edge_label_fontsize: float = 8,
    # Gene Label Options
    show_gene_labels: bool = False,
    n_top_genes: int = 2,
    gene_label_offset: float = 0.3,
    gene_label_fontsize: float = 10,
    gene_label_threshold: float = 0.05,
    # Level Label Options
    level_label_offset: float = 0.5,
    level_label_fontsize: float = 12,
    # Title Options
    title: str = "Hierarchical Leiden Clustering",
    title_fontsize: float = 16,
) -> nx.DiGraph:
    """
    Plot a hierarchical clustering decision tree based on multiple resolutions.

    This function performs Leiden clustering at different resolutions (if not already computed),
    constructs a decision tree representing the hierarchical relationships between clusters,
    and visualizes it as a directed graph. Nodes represent clusters at different resolutions,
    edges represent transitions between clusters, and edge weights indicate the proportion of
    cells transitioning from a parent cluster to a child cluster.

    Params
    ------
    adata
        The annotated data matrix containing clustering results in `adata.uns["cluster_resolution_cluster_data"]`
        and top genes in `adata.uns["cluster_resolution_top_genes"]`. Typically populated by
        `sc.tl.cluster_resolution_finder`.
    prefix
        Prefix for clustering keys in `adata.obs` (e.g., "leiden_res_").
    resolutions
        List of resolution values for Leiden clustering.
    orientation
        Orientation of the tree: "vertical" or "horizontal".
    node_spacing
        Horizontal spacing between nodes at the same level (in data coordinates).
    level_spacing
        Vertical spacing between resolution levels (in data coordinates).
    barycenter_sweeps
        Number of barycenter-based reordering sweeps to minimize edge crossings.
    use_reingold_tilford
        Whether to use the Reingold-Tilford layout algorithm (requires `igraph`).
    output_path
        Path to save the figure (e.g., "cluster_tree.png"). Supports PNG, PDF, SVG.
    draw
        Whether to display the plot using `plt.show()`.
    figsize
        Figure size as (width, height) in inches.
    dpi
        Resolution for saving the figure (dots per inch).
    node_size
        Base size for nodes in points^2 (area of the node).
    node_color
        Color specification for nodes: "prefix" (color by resolution level) or a single color.
    node_colormap
        Custom colormap for nodes, as a list of colors (one per resolution level).
    node_label_fontsize
        Font size for node labels (e.g., cluster numbers).
    edge_color
        Color specification for edges: "parent" (inherit parent node color), "samples" (by weight), or a single color.
    edge_curvature
        Curvature of edges (intensity of the S-shape).
    edge_threshold
        Minimum weight (proportion of cells) required to draw an edge.
    show_weight
        Whether to show edge weights as labels on the edges.
    edge_label_threshold
        Minimum weight required to label an edge with its weight.
    edge_label_position
        Position of the edge weight label along the edge (0.0 to 1.0).
    edge_label_fontsize
        Font size for edge weight labels.
    show_gene_labels
        Whether to show gene labels below child nodes.
    n_top_genes
        Number of top genes to display for each (parent, child) pair.
    gene_label_offset
        Vertical offset for gene labels below nodes (in data coordinates).
    gene_label_fontsize
        Font size for gene labels.
    gene_label_threshold
        Minimum weight required to display a gene label for a (parent, child) pair.
    level_label_offset
        Horizontal buffer space between level labels and the leftmost node.
    level_label_fontsize
        Font size for level labels (e.g., "leiden_res_0.0").
    title
        Title of the plot.
    title_fontsize
        Font size for the plot title.

    Returns
    -------
    G
        The directed graph representing the hierarchical clustering, with nodes and edges
        annotated with resolution levels and weights.

    Notes
    -----
    This function requires the `igraph` library for Leiden clustering, which is included in the
    `leiden` extra. Install it with: ``pip install scanpy[leiden]``.

    If clustering results are not already present in `adata.obs`, the function will run
    `sc.tl.leiden` for the specified resolutions, which requires `sc.pp.neighbors` to be
    run first.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=0.0, key_added="leiden_res_0.0")
        sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res_0.5")
        sc.pl.cluster_decision_tree(adata, resolutions=[0.0, 0.5])
    """
    # Validate input parameters
    if (
        not isinstance(figsize, tuple | list)
        or len(figsize) != 2
        or any(dim <= 0 for dim in figsize)
    ):
        msg = "figsize must be a tuple of two positive numbers (width, height)."
        raise ValueError(msg)
    if dpi <= 0:
        msg = "dpi must be a positive number."
        raise ValueError(msg)
    if node_size <= 0:
        msg = "node_size must be a positive number."
        raise ValueError(msg)
    if edge_threshold < 0 or edge_label_threshold < 0:
        msg = "edge_threshold and edge_label_threshold must be non-negative."
        raise ValueError(msg)

    # Retrieve clustering data from adata.uns
    if "cluster_resolution_cluster_data" not in adata.uns:
        msg = "adata.uns['cluster_resolution_cluster_data'] not found. Run sc.tl.cluster_resolution_finder first."
        raise ValueError(msg)
    data = adata.uns["cluster_resolution_cluster_data"]

    # Validate that data has the required columns
    cluster_columns = [f"{prefix}{res}" for res in resolutions]
    missing_columns = [col for col in cluster_columns if col not in data.columns]
    if missing_columns:
        msg = f"Clustering results for resolutions {missing_columns} not found in adata.uns['cluster_resolution_cluster_data']."
        raise ValueError(msg)

    # Retrieve top genes from adata.uns
    if show_gene_labels:
        if "cluster_resolution_top_genes" not in adata.uns:
            msg = "adata.uns['cluster_resolution_top_genes'] not found. Run sc.tl.cluster_resolution_finder first or disable show_gene_labels."
            raise ValueError(msg)
        top_genes_dict = adata.uns["cluster_resolution_top_genes"]
    else:
        top_genes_dict = None

    # Build the graph
    G = build_cluster_graph(data, prefix, edge_threshold)

    # Compute node positions
    pos = compute_cluster_layout(
        G,
        node_spacing,
        level_spacing,
        orientation,
        barycenter_sweeps,
        use_reingold_tilford=use_reingold_tilford,
    )

    # Draw the visualization if requested
    if draw or output_path:
        draw_cluster_tree(
            G,
            pos,
            data,
            prefix,
            resolutions,
            output_path=output_path,
            draw=draw,
            figsize=figsize,
            dpi=dpi,
            node_size=node_size,
            node_color=node_color,
            node_colormap=node_colormap,
            node_label_fontsize=node_label_fontsize,
            edge_color=edge_color,
            edge_curvature=edge_curvature,
            edge_threshold=edge_threshold,
            show_weight=show_weight,
            edge_label_threshold=edge_label_threshold,
            edge_label_position=edge_label_position,
            edge_label_fontsize=edge_label_fontsize,
            top_genes_dict=top_genes_dict,
            show_gene_labels=show_gene_labels,
            n_top_genes=n_top_genes,
            gene_label_offset=gene_label_offset,
            gene_label_fontsize=gene_label_fontsize,
            gene_label_threshold=gene_label_threshold,
            level_label_offset=level_label_offset,
            level_label_fontsize=level_label_fontsize,
            title=title,
            title_fontsize=title_fontsize,
        )

    return G

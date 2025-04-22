from __future__ import annotations

import math
from typing import TYPE_CHECKING, TypedDict, cast

import igraph as ig
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path

if TYPE_CHECKING:
    from typing import NotRequired

    import networkx as nx
    import pandas as pd
    from anndata import AnnData


class OutputSettings(TypedDict):
    output_path: NotRequired[str | None]
    draw: NotRequired[bool]
    figsize: NotRequired[tuple[float, float] | None]
    dpi: NotRequired[int | None]


class NodeStyle(TypedDict):
    node_size: NotRequired[float]
    node_color: NotRequired[str]
    node_colormap: NotRequired[list[str] | None]
    node_label_fontsize: NotRequired[float]


class EdgeStyle(TypedDict):
    edge_color: NotRequired[str]
    edge_curvature: NotRequired[float]
    edge_threshold: NotRequired[float]
    show_weight: NotRequired[bool]
    edge_label_threshold: NotRequired[float]
    edge_label_position: NotRequired[float]
    edge_label_fontsize: NotRequired[float]


class GeneLabelSettings(TypedDict):
    show_gene_labels: NotRequired[bool]
    n_top_genes: NotRequired[int]
    gene_label_threshold: NotRequired[float]
    gene_label_style: NotRequired[dict[str, float]]
    top_genes_dict: NotRequired[dict[tuple[str, str], list[str]] | None]


class LevelLabelStyle(TypedDict):
    level_label_offset: NotRequired[float]
    level_label_fontsize: NotRequired[float]


class TitleStyle(TypedDict):
    title: NotRequired[str]
    title_fontsize: NotRequired[float]


class LayoutSettings(TypedDict):
    node_spacing: NotRequired[float]
    level_spacing: NotRequired[float]
    orientation: NotRequired[str]
    barycenter_sweeps: NotRequired[int]
    use_reingold_tilford: NotRequired[bool]


class ClusteringSettings(TypedDict):
    prefix: NotRequired[str]
    edge_threshold: NotRequired[float]


class ClusterTreePlotter:
    def __init__(
        self,
        adata: AnnData,
        resolutions: list[float],
        *,
        output_settings: OutputSettings | None = None,
        node_style: NodeStyle | None = None,
        edge_style: EdgeStyle | None = None,
        gene_label_settings: GeneLabelSettings | None = None,
        level_label_style: LevelLabelStyle | None = None,
        title_style: TitleStyle | None = None,
        layout_settings: LayoutSettings | None = None,
        clustering_settings: ClusteringSettings | None = None,
    ):
        """
        Initialize the cluster tree plotter.

        Parameters
        ----------
        adata
            AnnData object with clustering results.
        resolutions
            List of resolution values.
        output_settings
            Output settings (output_path, draw, figsize, dpi).
        node_style
            Node styling (node_size, node_color, node_colormap, node_label_fontsize).
        edge_style
            Edge styling (edge_color, edge_curvature, edge_threshold, ...).
        gene_label_settings
            Gene label settings (show_gene_labels, n_top_genes, ...).
        level_label_style
            Level label settings (level_label_offset, level_label_fontsize).
        title_style
            Title settings (title, title_fontsize).
        layout_settings
            Layout settings (node_spacing, level_spacing).
        clustering_settings
            Clustering settings (prefix).
        """
        self.adata = adata
        self.resolutions = resolutions
        self.output_settings = self._merge_with_default(
            output_settings, self.default_output_settings()
        )
        self.node_style = self._merge_with_default(
            node_style, self.default_node_style()
        )
        self.edge_style = self._merge_with_default(
            edge_style, self.default_edge_style()
        )
        self.gene_label_settings = self._merge_with_default(
            gene_label_settings, self.default_gene_label_settings()
        )
        self.level_label_style = self._merge_with_default(
            level_label_style, self.default_level_label_style()
        )
        self.title_style = self._merge_with_default(
            title_style, self.default_title_style()
        )
        self.layout_settings = self._merge_with_default(
            layout_settings, self.default_layout_settings()
        )
        self.clustering_settings = self._merge_with_default(
            clustering_settings, self.default_clustering_settings()
        )

        self.settings = {}
        self.settings["output"] = self.output_settings
        self.settings["node"] = self.node_style
        self.settings["edge"] = self.edge_style
        self.settings["gene_label"] = self.gene_label_settings
        self.settings["level_label"] = self.level_label_style
        self.settings["title"] = self.title_style
        self.settings["layout"] = self.layout_settings
        self.settings["clustering"] = self.clustering_settings

        # Initialize attributes
        self.G = None
        self.pos = None
        self.ax = plt.gca()  # Initialize self.ax with the current axis
        self.fig = None

    def _merge_with_default(self, user_dict, default_dict):
        return {**default_dict, **(user_dict or {})}

    @staticmethod
    def default_output_settings() -> OutputSettings:
        return {"output_path": None, "draw": False, "figsize": (12, 6), "dpi": 300}

    @staticmethod
    def default_node_style() -> NodeStyle:
        return {
            "node_size": 500,
            "node_color": "prefix",
            "node_colormap": None,
            "node_label_fontsize": 12,
        }

    @staticmethod
    def default_edge_style() -> EdgeStyle:
        return {
            "edge_color": "parent",
            "edge_curvature": 0.01,
            "edge_threshold": 0.01,
            "show_weight": True,
            "edge_label_threshold": 0.05,
            "edge_label_position": 0.8,
            "edge_label_fontsize": 8,
        }

    @staticmethod
    def default_gene_label_settings() -> GeneLabelSettings:
        return {
            "show_gene_labels": False,
            "n_top_genes": 2,
            "gene_label_threshold": 0.001,
            "gene_label_style": {"offset": 0.5, "fontsize": 8},
            "top_genes_dict": None,
        }

    @staticmethod
    def default_level_label_style() -> LevelLabelStyle:
        return {"level_label_offset": 15, "level_label_fontsize": 12}

    @staticmethod
    def default_title_style() -> TitleStyle:
        return {"title": "Hierarchical Leiden Clustering", "title_fontsize": 20}

    @staticmethod
    def default_layout_settings() -> LayoutSettings:
        return {
            "node_spacing": 5.0,
            "level_spacing": 1.5,
            "orientation": "vertical",
            "barycenter_sweeps": 2,
            "use_reingold_tilford": False,
        }

    @staticmethod
    def default_clustering_settings() -> ClusteringSettings:
        return {"prefix": "leiden_res_", "edge_threshold": 0.05}

    def build_cluster_graph(self) -> None:
        """
        Build a directed graph representing hierarchical clustering.

        Uses self.adata.obs, self.settings["clustering"]["prefix"], and self.settings["clustering"]["edge_threshold"].
        Stores the graph in self.G and updates top_genes_dict.
        """
        import networkx as nx

        prefix = self.settings["clustering"]["prefix"]
        edge_threshold = self.settings["clustering"]["edge_threshold"]
        data = self.adata.obs

        # Validate input data
        matching_columns = [col for col in data.columns if col.startswith(prefix)]
        if not matching_columns:
            msg = f"No columns found with prefix '{prefix}' in the DataFrame."
            raise ValueError(msg)

        self.G = nx.DiGraph()

        # Extract resolutions from column names
        resolutions_col = [col[len(prefix) :] for col in matching_columns]
        resolutions_col = sorted(
            [float(r) for r in resolutions_col if r.replace(".", "", 1).isdigit()]
        )

        # Add nodes with resolution attribute for layout
        for i, res in enumerate(resolutions_col):
            clusters = data[f"{prefix}{res}"].unique()
            for cluster in sorted(clusters):
                node = f"{res}_C{cluster}"
                self.G.add_node(node, resolution=i, cluster=cluster)

        # Build edges between consecutive resolutions
        for i in range(len(resolutions_col) - 1):
            res1 = f"{prefix}{resolutions_col[i]}"
            res2 = f"{prefix}{resolutions_col[i + 1]}"

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
                parent_node = f"{resolutions_col[i]}_C{parent}"
                child_node = f"{resolutions_col[i + 1]}_C{child}"
                if frac >= edge_threshold:
                    self.G.add_edge(parent_node, child_node, weight=frac)

        self.settings["gene_label"]["top_genes_dict"] = self.adata.uns.get(
            "top_genes_dict", {}
        )

    def compute_cluster_layout(self) -> dict[str, tuple[float, float]]:
        """Compute node positions for the cluster decision tree with crossing minimization."""
        import networkx as nx

        if self.G is None:
            msg = "Graph is not initialized. Call build_graph() first."
            raise ValueError(msg)

        use_reingold_tilford = self.settings["layout"]["use_reingold_tilford"]
        node_spacing = self.settings["layout"]["node_spacing"]
        level_spacing = self.settings["layout"]["level_spacing"]
        orientation = self.settings["layout"]["orientation"]
        barycenter_sweeps = self.settings["layout"]["barycenter_sweeps"]
        # Step 1: Apply Reingold-Tilford layout or fallback to multipartite layout
        if use_reingold_tilford:
            pos = self._apply_reingold_tilford_layout(self.G, node_spacing)
        else:
            pos = nx.multipartite_layout(
                self.G, subset_key="resolution", scale=int(node_spacing)
            )

        # Step 2: Adjust orientation
        pos = self._adjust_orientation(
            pos=cast("dict[str, tuple[float, float]]", pos), orientation=orientation
        )

        # Step 3: Increase vertical spacing
        pos = self._adjust_vertical_spacing(pos, level_spacing)

        # Step 4: Barycenter-based reordering to minimize edge crossings
        pos = self._barycenter_sweep(
            self.G, pos, self.resolutions, node_spacing, barycenter_sweeps
        )

        # Step 5: Optimize node ordering
        filtered_edges = [
            (u, v, d["weight"])
            for u, v, d in self.G.edges(data=True)
            if d["weight"] >= 0.02
        ]
        edges = [(u, v) for u, v, w in filtered_edges]
        edges_set = set(edges)
        if len(edges_set) < len(edges):
            print(
                f"Warning: Found {len(edges) - len(edges_set)} duplicate edges in the visualization."
            )
        edges = list(edges_set)
        self._optimize_node_ordering(self.G, pos, edges, self.resolutions)
        self.pos = pos
        return self.pos

    def _apply_reingold_tilford_layout(
        self, G: nx.DiGraph, node_spacing: float
    ) -> dict[str, tuple[float, float]]:
        """Apply Reingold-Tilford layout to the graph."""
        import networkx as nx

        try:
            nodes = list(G.nodes)
            edges = [(u, v) for u, v in G.edges()]
            g = ig.Graph()
            g.add_vertices(nodes)
            g.add_edges([(nodes.index(u), nodes.index(v)) for u, v in edges])
            layout = g.layout_reingold_tilford(root=[0])
            return dict(zip(nodes, layout.coords, strict=False))
        except ImportError as e:
            print(
                f"igraph not installed or failed: {e}. Falling back to multipartite_layout."
            )
            return dict(
                nx.multipartite_layout(
                    G, subset_key="resolution", scale=int(node_spacing)
                )
            )

    def _adjust_orientation(
        self, pos: dict[str, tuple[float, float]], orientation: str
    ) -> dict[str, tuple[float, float]]:
        """Adjust the node positions for the specified orientation."""
        if orientation == "vertical":
            return {node: (y, -x) for node, (x, y) in pos.items()}
        return pos

    def _adjust_vertical_spacing(
        self, pos: dict[str, tuple[float, float]], level_spacing: float
    ) -> dict[str, tuple[float, float]]:
        """Increase vertical spacing between nodes at different levels."""
        new_pos = {}
        for node, (x, y) in pos.items():
            new_y = y * level_spacing
            new_pos[node] = (x, new_y)
        return new_pos

    def _barycenter_sweep(
        self,
        G: nx.DiGraph,
        pos: dict[str, tuple[float, float]],
        resolutions: list,
        node_spacing: float,
        barycenter_sweeps: int,
    ) -> dict[str, tuple[float, float]]:
        """Perform barycenter-based reordering to minimize edge crossings."""
        for _sweep in range(barycenter_sweeps):
            # Downward sweep: Adjust nodes based on parent positions
            pos = self._downward_sweep(G, pos, resolutions, node_spacing)
            # Upward sweep: Adjust nodes based on child positions
            pos = self._upward_sweep(G, pos, resolutions, node_spacing)
        self.pos = pos
        return pos

    def _downward_sweep(
        self, G: nx.DiGraph, pos: dict, resolutions: list, node_spacing: float
    ) -> dict[str, tuple[float, float]]:
        """Perform downward sweep in barycenter reordering."""
        for res in resolutions[1:]:
            nodes_at_level = [node for node in G.nodes if node.startswith(f"{res}_C")]
            node_to_barycenter = {}
            for node in nodes_at_level:
                parents = list(G.predecessors(node))
                barycenter = (
                    np.mean([pos[parent][0] for parent in parents]) if parents else 0
                )
                node_to_barycenter[node] = barycenter
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
            for node, x in zip(sorted_nodes, x_positions, strict=True):
                pos[node] = (x, y_level)
        return pos

    def _upward_sweep(
        self,
        G: nx.DiGraph,
        pos: dict[str, tuple[float, float]],
        resolutions: list,
        node_spacing: float,
    ) -> dict[str, tuple[float, float]]:
        """Perform upward sweep in barycenter reordering."""
        for res in reversed(resolutions[:-1]):
            nodes_at_level = [node for node in G.nodes if node.startswith(f"{res}_C")]
            node_to_barycenter = {}
            for node in nodes_at_level:
                children = list(G.successors(node))
                barycenter = (
                    np.mean([pos[child][0] for child in children]) if children else 0
                )
                node_to_barycenter[node] = barycenter
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
            for node, x in zip(sorted_nodes, x_positions, strict=True):
                pos[node] = (x, y_level)
        return pos

    def _optimize_node_ordering(
        self,
        G: nx.DiGraph,
        pos: dict[str, tuple[float, float]],
        edges: list[tuple[str, str]],
        resolutions: list,
        max_iterations=10,
    ) -> None:
        """Optimize node ordering at each level to minimize edge crossings by swapping adjacent nodes."""
        # Group nodes by resolution level
        level_nodes = {
            res_idx: [
                node for node in G.nodes if G.nodes[node]["resolution"] == res_idx
            ]
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
                    current_crossings = self._count_crossings(G, pos, edges)

                    # Swap positions and compute new crossings
                    pos[node1] = (x2, y1)
                    pos[node2] = (x1, y2)
                    new_crossings = self._count_crossings(G, pos, edges)

                    # If swapping reduces crossings, keep the swap
                    if new_crossings < current_crossings:
                        nodes[i], nodes[i + 1] = nodes[i + 1], nodes[i]
                        improved = True
                    else:
                        # Revert the swap if it doesn't improve crossings
                        pos[node1] = (x1, y1)
                        pos[node2] = (x2, y2)

                iteration += 1

    def _count_crossings(
        self,
        G: nx.DiGraph,
        pos: dict[str, tuple[float, float]],
        edges: list[tuple[str, str]],
    ) -> int:
        """Count the number of edge crossings in the graph based on node positions."""
        crossings = 0
        for i, (u1, v1) in enumerate(edges):
            for _j, (u2, v2) in enumerate(edges[i + 1 :], start=i + 1):
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

    def draw_cluster_tree(self) -> None:
        """Draw a hierarchical cluster tree with nodes, edges, and labels."""
        if self.G is None or self.pos is None:
            msg = "Graph or positions not initialized. Call build_graph() and compute_cluster_layout() first."
            raise ValueError(msg)
        if "cluster_resolution_cluster_data" not in self.adata.uns:
            msg = "adata.uns['cluster_resolution_cluster_data'] not found."
            raise ValueError(msg)

        import networkx as nx

        # Retrieve settings
        settings = self._get_draw_settings()
        data = settings["data"]
        prefix = settings["prefix"]

        # Step 1: Compute Cluster Sizes, Node Sizes, and Node Colors
        cluster_sizes = self._compute_cluster_sizes(data, prefix, self.resolutions)
        node_sizes = self._scale_node_sizes(
            data, prefix, self.resolutions, cluster_sizes, settings["node_size"]
        )
        color_schemes = self._generate_node_color_schemes(
            data,
            prefix,
            self.resolutions,
            settings["node_color"],
            settings["node_colormap"],
        )
        node_colors = self._assign_node_colors(
            data, prefix, self.resolutions, settings["node_color"], color_schemes
        )
        # Step 2: Set up the plot figure and axis
        self.fig = plt.figure(figsize=settings["figsize"], dpi=settings["dpi"])
        self.ax = self.fig.add_subplot(111)
        # Step 3: Compute Edge Weights, Edge Colors
        edges, weights, edge_colors = self._compute_edge_weights_colors(
            self.G, settings["edge_threshold"], settings["edge_color"], node_colors
        )
        # Step 4: Draw Nodes and Node Labels
        node_styles = {"colors": node_colors, "sizes": node_sizes}
        node_labels, gene_labels = self._draw_nodes_and_labels(
            self.G,
            self.pos,
            self.resolutions,
            node_styles=node_styles,
            data=data,
            prefix=prefix,
            top_genes_dict=self.adata.uns.get("cluster_resolution_top_genes", {}),
            show_gene_labels=settings["show_gene_labels"],
            n_top_genes=settings["n_top_genes"],
            gene_label_threshold=settings["gene_label_threshold"],
        )
        nx.draw_networkx_labels(
            self.G,
            self.pos,
            labels=node_labels,
            font_size=int(settings["node_label_fontsize"]),
            font_color="black",
        )
        # Step 5: Draw Gene Labels
        gene_label_bottoms = {}
        if settings["show_gene_labels"] and gene_labels:
            gene_label_bottoms = self._draw_gene_labels(
                self.ax,
                self.pos,
                gene_labels,
                node_sizes=node_sizes,
                node_colors=node_colors,
                offset=settings["gene_label_offset"],
                fontsize=settings["gene_label_fontsize"],
            )
        # Step 6: Build and Draw Edge Labels
        edge_labels = self._build_edge_labels(
            self.G, settings["edge_threshold"], settings["edge_label_threshold"]
        )
        edge_label_style = {
            "position": settings["edge_label_position"],
            "fontsize": settings["edge_label_fontsize"],
        }
        self._draw_edges_with_labels(
            self.ax,
            self.pos,
            edges,
            weights,
            edge_colors=edge_colors,
            node_sizes=node_sizes,
            gene_label_bottoms=gene_label_bottoms,
            show_gene_labels=settings["show_gene_labels"],
            edge_labels=edge_labels,
            edge_label_style=edge_label_style,
        )
        # Step 7: Draw Level Labels
        self._draw_level_labels(
            resolutions=self.resolutions,
            pos=self.pos,
            data=self.adata.uns["cluster_resolution_cluster_data"],
            prefix=prefix,
            level_label_offset=settings["level_label_offset"],
            level_label_fontsize=settings["level_label_fontsize"],
        )
        # Step 8: Final Plot Settings
        self.ax.set_title(settings["title"], fontsize=settings["title_fontsize"])
        self.ax.axis("off")
        # Save or show the plot
        if settings["output_path"]:
            plt.savefig(settings["output_path"], bbox_inches="tight")
        if settings["draw"]:
            plt.show()

    def _get_draw_settings(self) -> dict:
        """Retrieve settings for drawing the cluster tree."""
        data = self.adata.uns["cluster_resolution_cluster_data"]
        return {
            "data": data,
            "prefix": self.settings["clustering"]["prefix"],
            "node_size": self.settings["node"]["node_size"],
            "node_color": self.settings["node"]["node_color"],
            "node_colormap": self.settings["node"]["node_colormap"],
            "figsize": self.settings["output"]["figsize"],
            "dpi": self.settings["output"]["dpi"],
            "edge_threshold": self.settings["edge"]["edge_threshold"],
            "edge_color": self.settings["edge"]["edge_color"],
            "show_gene_labels": self.settings["gene_label"]["show_gene_labels"],
            "n_top_genes": self.settings["gene_label"]["n_top_genes"],
            "gene_label_threshold": self.settings["gene_label"]["gene_label_threshold"],
            "node_label_fontsize": self.settings["node"]["node_label_fontsize"],
            "gene_label_offset": self.settings["gene_label"]["gene_label_style"][
                "offset"
            ],
            "gene_label_fontsize": self.settings["gene_label"]["gene_label_style"][
                "fontsize"
            ],
            "edge_label_threshold": self.settings["edge"]["edge_label_threshold"],
            "edge_label_position": self.settings["edge"]["edge_label_position"],
            "edge_label_fontsize": self.settings["edge"]["edge_label_fontsize"],
            "level_label_offset": self.settings["level_label"]["level_label_offset"],
            "level_label_fontsize": self.settings["level_label"][
                "level_label_fontsize"
            ],
            "title": self.settings["title"]["title"],
            "title_fontsize": self.settings["title"]["title_fontsize"],
            "output_path": self.settings["output"]["output_path"],
            "draw": self.settings["output"]["draw"],
        }

    def _compute_cluster_sizes(
        self, data: pd.DataFrame, prefix: str, resolutions: list
    ) -> dict[str, int]:
        """Compute cluster sizes for each node."""
        cluster_sizes = {}
        for res in resolutions:
            res_key = f"{prefix}{res}"
            counts = data[res_key].value_counts()
            for cluster, count in counts.items():
                node = f"{res}_C{cluster}"
                cluster_sizes[node] = count
        return cluster_sizes

    def _scale_node_sizes(
        self,
        data: pd.DataFrame,
        prefix: str,
        resolutions: list,
        cluster_sizes: dict[str, int],
        node_size: float,
    ) -> dict[str, float]:
        """Scale node sizes based on cluster sizes and node_size setting."""
        node_sizes = {}
        for res in resolutions:
            nodes_at_level = [
                f"{res}_C{cluster}" for cluster in data[f"{prefix}{res}"].unique()
            ]
            sizes = np.array([cluster_sizes[node] for node in nodes_at_level])
            if len(sizes) > 1:
                min_size, max_size = sizes.min(), sizes.max()
                if min_size != max_size:
                    normalized_sizes = 0.5 + (sizes - min_size) / (max_size - min_size)
                else:
                    normalized_sizes = np.ones_like(sizes) * 0.5
                scaled_sizes = normalized_sizes * node_size
            else:
                scaled_sizes = np.array([node_size])
            if len(nodes_at_level) != len(scaled_sizes):
                msg = (
                    f"Length mismatch at resolution {res}: "
                    f"{len(nodes_at_level)} nodes vs {len(scaled_sizes)} sizes"
                )
                raise ValueError(msg)
            node_sizes.update(dict(zip(nodes_at_level, scaled_sizes, strict=False)))
        return node_sizes

    def _generate_node_color_schemes(
        self,
        data: pd.DataFrame,
        prefix: str,
        resolutions: list,
        node_color: str | None,
        node_colormap: list[str] | None,
    ) -> list[str] | dict[str, list] | None:
        """Generate color schemes for nodes."""
        if node_color != "prefix":
            return None

        if node_colormap is None:
            return {
                r: sns.color_palette("Set3", n_colors=data[f"{prefix}{r}"].nunique())
                for r in resolutions
            }

        if len(node_colormap) < len(resolutions):
            node_colormap = list(node_colormap) + [
                node_colormap[i % len(node_colormap)]
                for i in range(len(resolutions) - len(node_colormap))
            ]

        color_schemes = {}
        for i, r in enumerate(resolutions):
            color_spec = node_colormap[i]
            if (isinstance(color_spec, str) and mcolors.is_color_like(color_spec)) or (
                isinstance(color_spec, tuple)
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
                        f"Warning: '{color_spec}' is not valid for {r}. Using 'Set3'."
                    )
                    color_schemes[r] = sns.color_palette(
                        "Set3", n_colors=data[f"{prefix}{r}"].nunique()
                    )
        return color_schemes

    def _assign_node_colors(
        self,
        data: pd.DataFrame,
        prefix: str,
        resolutions: list,
        node_color: str,
        color_schemes: list[str] | dict[str, list] | None,
    ) -> dict[str, str]:
        node_colors = {}
        for res in resolutions:
            clusters = data[f"{prefix}{res}"].unique()
            for cluster in clusters:
                node = f"{res}_C{cluster}"
                if node_color == "prefix":
                    if color_schemes is None:
                        msg = "color_schemes is None but node_color='prefix'"
                        raise RuntimeError(msg)
                    colors = color_schemes[res]
                    node_colors[node] = (
                        colors[0]
                        if len(colors) == 1
                        else colors[int(cluster) % len(colors)]
                    )
                else:
                    node_colors[node] = node_color
        return node_colors

    def _compute_edge_weights_colors(
        self,
        G: nx.DiGraph,
        edge_threshold: float,
        edge_color: str,
        node_colors: dict,
    ) -> tuple[list, list, list]:
        """Compute edge weights and colors based on the graph and edge_threshold."""
        edges = [
            (u, v) for u, v, d in G.edges(data=True) if d["weight"] >= edge_threshold
        ]
        weights = [
            max(d["weight"] * 5, 1.0)
            for u, v, d in G.edges(data=True)
            if d["weight"] >= edge_threshold
        ]
        edge_colors = []
        for u, v in edges:
            d = G[u][v]
            if edge_color == "parent":
                edge_colors.append(node_colors[u])
            elif edge_color == "samples":
                edge_colors.append(plt.cm.get_cmap("viridis")(d["weight"] / 5))
            else:
                edge_colors.append(edge_color)
        return edges, weights, edge_colors

    def _draw_nodes_and_labels(
        self,
        G: nx.DiGraph,
        pos: dict[str, tuple[float, float]],
        resolutions: list,
        *,
        node_styles: dict,
        data: pd.DataFrame,
        prefix: str,
        top_genes_dict: dict[tuple[str, str], list[str]],
        show_gene_labels: bool,
        n_top_genes: int,
        gene_label_threshold: float,
    ) -> tuple[dict, dict]:
        """Draw the nodes and their labels."""
        import networkx as nx

        node_colors = node_styles["colors"]
        node_sizes = node_styles["sizes"]
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
        return node_labels, gene_labels

    def _draw_gene_labels(
        self,
        ax,
        pos: dict[str, tuple[float, float]],
        gene_labels: dict[str, str],
        *,
        node_sizes: dict[str, float],
        node_colors: dict[str, str],
        offset: float = 0.2,
        fontsize: float = 8,
    ) -> dict[str, float]:
        """Draw gene labels in boxes below nodes with matching boundary colors."""
        gene_label_bottoms = {}
        for node, label in gene_labels.items():
            if label:
                x, y = pos[node]
                # Compute the node radius in data coordinates
                radius = math.sqrt(node_sizes[node] / math.pi)
                _fig_width, fig_height = ax.figure.get_size_inches()
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

    def _build_edge_labels(
        self, G: nx.DiGraph, edge_threshold: float, edge_label_threshold: float
    ) -> dict:
        """Build the edge labels to display on the plot."""
        edge_labels = {
            (u, v): f"{w:.2f}"
            for u, v, w in [
                (u, v, d["weight"])
                for u, v, d in G.edges(data=True)
                if d["weight"] >= edge_threshold
            ]
            if w >= edge_label_threshold
        }
        return edge_labels

    def _draw_edges_with_labels(
        self,
        ax,
        pos: dict[str, tuple[float, float]],
        edges: list,
        weights: list,
        *,
        edge_colors: list,
        node_sizes: dict,
        gene_label_bottoms: dict,
        show_gene_labels: bool,
        edge_labels: dict,
        edge_label_style: dict,
    ) -> None:
        """Draw edges with labels using Bezier curves."""
        edge_label_position = edge_label_style["position"]
        edge_label_fontsize = edge_label_style["fontsize"]
        for (u, v), w, e_color in zip(edges, weights, edge_colors, strict=False):
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            radius_parent = math.sqrt(node_sizes[u] / math.pi)
            radius_child = math.sqrt(node_sizes[v] / math.pi)
            _fig_width, fig_height = ax.figure.get_size_inches()
            radius_parent_fig = radius_parent / (72 * fig_height)
            radius_child_fig = radius_child / (72 * fig_height)
            ylim = ax.get_ylim()
            data_height = ylim[0] - ylim[1]
            radius_parent_data = radius_parent_fig * data_height
            radius_child_data = radius_child_fig * data_height
            start_y = (
                gene_label_bottoms[u]
                if (show_gene_labels and u in gene_label_bottoms and edge_labels.get(u))
                else y1 - radius_parent_data
            )
            start_x = x1
            end_x, end_y = x2, y2 - radius_child_data

            p0, p1, p2, p3 = self._draw_curved_edge(
                ax,
                start_x,
                start_y,
                end_x,
                end_y,
                linewidth=w,
                color=e_color,
                edge_curvature=0.01,
            )

            if (u, v) in edge_labels and p0 is not None:
                t = edge_label_position
                point = self._evaluate_bezier(t, p0, p1, p2, p3)
                label_x, label_y = point[0], point[1]
                tangent = self._evaluate_bezier_tangent(t, p0, p1, p2, p3)
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

    def _draw_curved_edge(
        self,
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
        """Draw a gentle S-shaped curved edge between two points with an arrowhead. Retun a tuple of Bézier control points (p0, p1, p2, p3) for label positioning."""
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

        # Define control points for a single cubic Bézier curve with an S-shape, Place control points at 1/3 and 2/3 along the edge, with small perpendicular offsets
        offset = length * edge_curvature
        p1 = (
            p0 + (p3 - p0) / 3 + perp_vec * offset
        )  # First control point (bend outward)
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
        arrow = FancyArrowPatch(
            (end_x, end_y),
            (end_x, end_y),
            arrowstyle="->",
            mutation_scale=arrow_size,
            color=color,
            linewidth=linewidth,
            alpha=0.8,
        )
        ax.add_patch(arrow)

        return p0, p1, p2, p3

    def _evaluate_bezier(
        self, t: float, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
    ) -> np.ndarray:
        """Evaluate a cubic Bezier curve at parameter t."""
        if not 0 <= t <= 1:
            msg = "Parameter t must be in the range [0, 1]"
            raise ValueError(msg)

        t2 = t * t
        t3 = t2 * t
        mt = 1 - t
        mt2 = mt * mt
        mt3 = mt2 * mt
        return mt3 * p0 + 3 * mt2 * t * p1 + 3 * mt * t2 * p2 + t3 * p3

    def _evaluate_bezier_tangent(
        self, t: float, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
    ) -> np.ndarray:
        """Compute the tangent vector of a cubic Bezier curve at parameter t."""
        if not 0 <= t <= 1:
            msg = "Parameter t must be in the range [0, 1]"
            raise ValueError(msg)

        t2 = t * t
        mt = 1 - t
        mt2 = mt * mt
        return 3 * mt2 * (p1 - p0) + 6 * mt * t * (p2 - p1) + 3 * t2 * (p3 - p2)

    def _draw_level_labels(
        self,
        resolutions: list,
        pos: dict[str, tuple[float, float]],
        data: pd.DataFrame,
        *,
        prefix: str,
        level_label_offset: float,
        level_label_fontsize: float,
    ) -> None:
        """Draw level labels for each resolution in the plot."""
        level_positions = {}
        for node, (_x, y) in pos.items():
            res = node.split("_")[0]
            level_positions[res] = y

        cluster_counts = {}
        for res in resolutions:
            res_str = f"{res:.1f}"
            col_name = f"{prefix}{res_str}"
            if col_name not in data.columns:
                msg = f"Column {col_name} not found in data. Ensure clustering results are present."
                raise ValueError(msg)
            num_clusters = len(data[col_name].dropna().unique())
            cluster_counts[res_str] = num_clusters

        min_x = min(p[0] for p in pos.values())
        label_offset = min_x - level_label_offset
        for res in resolutions:
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

    @staticmethod
    def cluster_decision_tree(
        adata: AnnData,
        resolutions: list[float],
        *,
        output_settings: dict | OutputSettings | None = None,
        node_style: dict | NodeStyle | None = None,
        edge_style: dict | EdgeStyle | None = None,
        gene_label_settings: dict | GeneLabelSettings | None = None,
        level_label_style: dict | LevelLabelStyle | None = None,
        title_style: dict | TitleStyle | None = None,
        layout_settings: dict | LayoutSettings | None = None,
        clustering_settings: dict | ClusteringSettings | None = None,
    ) -> nx.DiGraph:
        """Plot a hierarchical clustering decision tree based on multiple resolutions.

        This static method performs Leiden clustering at different resolutions (if not already computed),
        constructs a decision tree representing hierarchical relationships between clusters,
        and visualizes it as a directed graph. Nodes represent clusters at different resolutions,
        edges represent transitions between clusters, and edge weights indicate the proportion of
        cells transitioning from a parent to a child cluster.

        Parameters
        ----------
        adata
            Annotated data matrix with clustering results in adata.uns["cluster_resolution_cluster_data"].
        resolutions
            List of resolution values for Leiden clustering.
        output_settings
            Dictionary with output options (output_path, draw, figsize, dpi).
        node_style
            Dictionary with node appearance (node_size, node_color, node_colormap, node_label_fontsize).
        edge_style
            Dictionary with edge appearance (edge_color, edge_curvature, edge_threshold, etc.).
        gene_label_settings
            Dictionary with gene label options (show_gene_labels, n_top_genes, etc.).
        level_label_style
            Dictionary with level label options (level_label_offset, level_label_fontsize).
        title_style
            Dictionary with title options (title, title_fontsize).
        layout_settings
            Dictionary with layout options (orientation, node_spacing, level_spacing, etc.).
        clustering_settings
            Dictionary with clustering options (prefix, edge_threshold).

        Returns
        -------
        Directed graph representing the hierarchical clustering.

        """
        # Run all validations
        ClusterTreePlotter._validate_parameters(output_settings, node_style, edge_style)
        ClusterTreePlotter._validate_clustering_data(
            adata, resolutions, clustering_settings
        )
        ClusterTreePlotter._validate_gene_labels(adata, gene_label_settings)

        # Initialize ClusterTreePlotter
        plotter = ClusterTreePlotter(
            adata,
            resolutions,
            output_settings=cast("OutputSettings", output_settings),
            node_style=cast("NodeStyle", node_style),
            edge_style=cast("EdgeStyle", edge_style),
            gene_label_settings=cast("GeneLabelSettings", gene_label_settings),
            level_label_style=cast("LevelLabelStyle", level_label_style),
            title_style=cast("TitleStyle", title_style),
            layout_settings=cast("LayoutSettings", layout_settings),
            clustering_settings=cast("ClusteringSettings", clustering_settings),
        )
        # Build graph and compute layout
        plotter.build_cluster_graph()
        plotter.compute_cluster_layout()

        # Draw if requested
        if (output_settings or {}).get("draw", True) or (output_settings or {}).get(
            "output_path"
        ):
            plotter.draw_cluster_tree()

        if plotter.G is None:
            msg = "Graph is not initialized. Ensure build_cluster_graph() has been called."
            raise ValueError(msg)
        return plotter.G

    @staticmethod
    def _validate_parameters(output_settings, node_style, edge_style):
        if output_settings:
            figsize = output_settings.get("figsize")
            if (
                not isinstance(figsize, tuple | list)
                or len(figsize) != 2
                or any(dim <= 0 for dim in figsize)
            ):
                msg = "figsize must be a tuple of two positive numbers (width, height)."
                raise ValueError(msg)

            dpi = output_settings.get("dpi", 0)
            if not isinstance(dpi, int | float) or dpi <= 0:
                msg = "dpi must be a positive number."
                raise ValueError(msg)

            if output_settings.get("draw") not in [True, False, None]:
                msg = "draw must be True, False, or None."
                raise ValueError(msg)

        if node_style:
            node_size_val = node_style.get("node_size")
            if node_size_val is not None and node_size_val <= 0:
                msg = "node_size must be a positive number."
                raise ValueError(msg)

        if edge_style and (
            (edge_style.get("edge_threshold", 0)) < 0
            or edge_style.get("edge_label_threshold", 0) < 0
        ):
            msg = "edge_threshold and edge_label_threshold must be non-negative."
            raise ValueError(msg)

    @staticmethod
    def _validate_clustering_data(adata, resolutions, clustering_settings):
        if "cluster_resolution_cluster_data" not in adata.uns:
            msg = "adata.uns['cluster_resolution_cluster_data'] not found. Run `sc.tl.cluster_resolution_finder` first."
            raise ValueError(msg)
        if not resolutions:
            msg = "You must provide a list of resolutions."
            raise ValueError(msg)

        prefix = (clustering_settings or {}).get("prefix", "leiden_res_")
        cluster_columns = [f"{prefix}{res}" for res in resolutions]
        data = adata.uns["cluster_resolution_cluster_data"]
        missing = [col for col in cluster_columns if col not in data.columns]
        if missing:
            msg = f"Missing clustering columns: {missing}"
            raise ValueError(msg)

    @staticmethod
    def _validate_gene_labels(adata, gene_label_settings):
        if (
            gene_label_settings
            and gene_label_settings.get("show_gene_labels", False)
            and "cluster_resolution_top_genes" not in adata.uns
        ):
            msg = "Gene labels requested but `adata.uns['cluster_resolution_top_genes']` not found. Run `sc.tl.cluster_resolution_finder` first."
            raise ValueError(msg)


cluster_decision_tree = ClusterTreePlotter.cluster_decision_tree

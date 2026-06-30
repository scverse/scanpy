from __future__ import annotations

from typing import TYPE_CHECKING, cast

import holoviews as hv
import numpy as np
import pandas as pd
from anndata.acc import GraphAcc, MultiAcc
from hv_anndata import A

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from hv_anndata import AdDim
    from scipy.sparse import csr_matrix  # noqa: TID251


def draw_graph(
    adata: AnnData,
    kdims: list[AdDim] | MultiAcc,
    edge_vdim: Literal["distances", "connectivities"] | GraphAcc = "connectivities",
    node_vdims: AdDim | list[AdDim] | None = None,
    *,
    neighbors_key: str = "neighbors",
) -> hv.Graph:
    """Draw a graph.

    Parameters
    ----------
    adata
        Annotated data matrix.
    kdims
        Key dimensions of the graph.
        Can just be ``A.obsm[<k>]`` or ``A.varm[<k>]`` to auto-select components.
    edge_vdim
        Edge value dimension.
        If ``"distances"``/``"connectivities"``,
        the graph data is retrieved like :func:`scanpy.pp.neighbors` stores it:
        ``A.uns[neighbors_key][f"{edge_vdim}_key"]``.
        Therefore ``.opts(edge_color=calculated_edge_vdim)`` is set by default.
    node_vdims
        Node value dimensions.
    neighbors_key
        Key in ``adata.uns`` where neighbors are stored.
        Used only if ``edge_vdim`` is ``"distances"``/``"connectivities"``.

    Returns
    -------
    Graph with colored edges.

    Examples
    --------

    ..  holoviews::
        :backends: bokeh

        import scanpy as sc
        import hv_anndata.plotting.scanpy as hv_sc
        from hv_anndata import data, register, A

        register()

        adata = data.pbmc68k_processed()
        sc.pp.neighbors(adata)

        hv_sc.draw_graph(
            adata, A.obsm["X_umap"], "distances", [A.obs["bulk_labels"]]
        ).opts(
            node_color=A.obs["bulk_labels"],
            node_cmap="tab10",
            aspect="square",
            show_legend=True,
            legend_position="right",
        )

    """
    adata = adata.copy()
    adata.obs["cell index"] = range(adata.n_obs)
    if isinstance(kdims, MultiAcc):
        kdims = [kdims[0], kdims[1]]
    if isinstance(edge_vdim, str):
        edge_vdim = A.obsp[adata.uns[neighbors_key][f"{edge_vdim}_key"]]
    elif not isinstance(edge_vdim, GraphAcc):
        msg = f"edge_vdim must be a string or `A.obsp[key]`, got {edge_vdim!r}."
        raise TypeError(msg)

    edges = cast("csr_matrix", getattr(adata, f"{edge_vdim.dim}p")[edge_vdim.k]).tocoo()
    nodes = hv.Nodes(adata, [*kdims, A.obs["cell index"]], node_vdims)
    return hv.Graph(((*edges.coords, edges.data), nodes), vdims=edge_vdim[:, :]).opts(
        edge_color=edge_vdim[:, :]
    )


def ranking(
    adata: AnnData,
    ref: AdDim,
    /,
    n_points: int = 10,
    *,
    include_lowest: bool = True,
    label_dim: AdDim | None = None,
) -> hv.Labels:
    """Plot (e.g. PCA) score ranking.

    Parameters
    ----------
    adata
        Annotated data matrix.
    ref
        Dimension containing scores to rank.
    n_points
        Number of points to plot.
    include_lowest
        Whether to include the lowest-scored names in addition to the highest-scored ones.
    label_dim
        Dimension to use for labels.
        The default is ``dim``’s axis index (e.g. ``A.obs.index`` for ``A.obs["scores"]``).

    Returns
    -------
    Holoviews plot with labels and points.

    Examples
    --------

    ..  holoviews::

        import hv_anndata.plotting.scanpy as hv_sc
        from hv_anndata import data, register, A

        register()

        adata = data.pbmc68k_processed()
        hv.Layout([
            hv_sc.ranking(adata, A.varm["PCs"][0]).opts(aspect=1.2),
            hv_sc.ranking(adata, A.varm["PCs"][0], include_lowest=False).opts(aspect=0.6),
        ]).opts(shared_axes=False)

    """
    [dim] = ref.dims
    if label_dim is None:
        label_dim = getattr(A, dim).index
    # full arrays
    scores = adata[ref]
    labels = adata[label_dim]

    # subset
    idx = np.argsort(scores)
    idx_top, idx_bot = idx[-n_points:][::-1], idx[:n_points][::-1]
    scores = np.r_[scores[idx_top], np.nan, scores[idx_bot]]
    labels = np.r_[labels[idx_top], ["⋯"], labels[idx_bot]]

    # prepare
    data = pd.DataFrame(
        dict(
            rank=np.arange(n_points * 2 + 1),
            score=np.where(np.isnan(scores), np.nanmean(scores), scores),
            dot=np.r_[[True] * n_points, False, [True] * n_points].astype(float),
            text=labels,
            align=np.r_[["right"] * n_points, ["center"], ["left"] * n_points],
        )
    )
    if not include_lowest:
        i = data[:n_points]["score"].diff().abs().argmax()
        data = data[:n_points].assign(
            align=np.r_[["right"] * i, ["left"] * (n_points - i)]
        )

    return hv.Labels(data, ["rank", "score"], ["text", "align"]).opts(
        angle=90, text_align="align", xticks=0
    ) * hv.Points(data[data["dot"] > 0], ["rank", "score"], ["text"])


def embedding_density(
    adata: AnnData,
    basis: MultiAcc,
    *,
    groupby: str | None = None,
    key: str | None = None,
) -> hv.Scatter | hv.NdLayout:
    """Plot embedding density.

    Parameters
    ----------
    adata
        Annotated data matrix.
    basis
        Embedding to plot (e.g. ``A.obsm["X_umap"]``).
    groupby
        ``groupby`` as specified in :func:`scanpy.tl.embedding_density`.
    key
        ``key_added`` as specified in :func:`scanpy.tl.embedding_density`.

    Returns
    -------
    Scatter plot if ``groupby is None``, else a layout of scatter plots.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        import hv_anndata.plotting.scanpy as hv_sc
        from hv_anndata import data, register, A

        register()

        adata = data.pbmc68k_processed()
        sc.tl.embedding_density(adata, basis="umap", groupby="phase")
        hv_sc.embedding_density(adata, A.obsm["X_umap"], groupby="phase")

    """
    basis_name = basis.k.removeprefix("X_")
    if key is None:
        key = f"{basis_name}_density{'' if groupby is None else f'_{groupby}'}"

    groupby_vdims = [] if groupby is None else [A.obs[groupby]]
    scatter = hv.Scatter(adata, basis[0], [basis[1], *groupby_vdims, A.obs[key]]).opts(
        color=A.obs[key],
        aspect="square",
        xlabel=f"{basis_name} 1",
        ylabel=f"{basis_name} 2",
        legend_position="right",
    )
    if groupby is None:
        return scatter
    return scatter.groupby(A.obs[groupby], hv.NdLayout)

"""Scanpy plots."""

from __future__ import annotations

import inspect
from collections.abc import Collection
from functools import partial, reduce, update_wrapper
from types import MappingProxyType
from typing import TYPE_CHECKING, assert_never, overload

import holoviews as hv
import numpy as np
import pandas as pd
from anndata.acc import GraphAcc
from fast_array_utils import stats
from hv_anndata import A, AdDim

import scanpy as sc
from scanpy.plotting._common import dot_area

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Mapping
    from typing import Literal

    from anndata import AnnData
    from anndata.acc import LayerAcc, MultiAcc
    from pandas.api.extensions import ExtensionArray

    # TODO: export in scanpy: https://github.com/scverse/scanpy/issues/3826
    AggType = Literal["count_nonzero", "mean", "sum", "var", "median"]


__all__ = [
    "diffmap",
    "dotplot",
    "heatmap",
    "matrixplot",
    "scatter",
    "stacked_violin",
    "tracksplot",
    "tsne",
    "umap",
    "violin",
]


def scatter(
    adata: AnnData,
    /,
    kdims: Collection[AdDim],
    vdims: Collection[AdDim] = (),
    *,
    color: AdDim | Collection[AdDim] | None = None,
) -> hv.Scatter | hv.Layout:
    """Shortcut for a scatter plot.

    Basically just

    >>> hv.Scatter(adata, kdims[0], [kdims[1], *vdims]).opts(aspect="square", ...)

    If ``color`` is set, it’s both added to ``vdims`` and in ``.opts(color=...)``.

    Parameters
    ----------
    adata
        The AnnData object.
    kdims
        A sequence containing the x and y dimension.
    vdims
        The value dimensions (``color`` will be added automatically).
    color
        The color dimension.
        If a collection of refs is given (e.g. ``A.obs[["v1", "v2"]]``),
        returns a :class:`~holoviews.Layout` with one scatter plot per ref.

    Returns
    -------
    A scatter plot object,
    or a :class:`~holoviews.Layout` of them if ``color`` is a collection.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()

        sc.pl.scatter(adata, A.X[:, ["PSAP", "C1QA"]], color=A.obs["bulk_labels"]).opts(
            cmap="tab10", show_legend=False
        )

    If you use faceting, make sure to apply options to the plots and not the :class:`~holoviews.Layout`:

    ..  holoviews::

        import holoviews as hv

        sc.pl.scatter(adata, A.obs[["n_counts", "n_genes"]], color=A.obs[["S_score", "G2M_score"]]).opts(
            hv.opts.Scatter(cmap="viridis", show_legend=False)
        )

    """
    try:
        i, j = kdims
    except ValueError:
        msg = "kdims must have length 2"
        raise ValueError(msg) from None

    if color is not None and not isinstance(color, AdDim):
        return _facet(color, lambda c: scatter(adata, kdims, vdims, color=c))

    if color is not None:
        vdims = [*vdims, color]
    sc = hv.Scatter(adata, i, [j, *vdims])
    if color is not None:
        sc = sc.opts(color=color)

    return sc.opts(aspect="square", legend_position="right")


def _scatter(
    kdims: Collection[AdDim],
    adata: AnnData,
    /,
    vdims: Collection[AdDim] = (),
    *,
    color: AdDim | Collection[AdDim] | None = None,
) -> hv.Scatter | hv.Layout:
    __tracebackhide__ = True
    return scatter(adata, kdims, vdims, color=color)


def _embedding(key: str, name: str, /) -> partial:
    p = partial(_scatter, A.obsm[key][:, [0, 1]])
    update_wrapper(p, scatter, updated=())
    p.__name__ = p.__qualname__ = key

    scatter_doc = inspect.getdoc(scatter) or ""
    examples = scatter_doc[scatter_doc.index("Examples\n-------") :]
    for needle in ['A.X[:, ["PSAP", "C1QA"]]', 'A.obs[["n_counts", "n_genes"]]']:
        needle = f"scatter(adata, {needle}"  # noqa: PLW2901
        assert needle in examples, (needle, examples)
        examples = examples.replace(needle, f"{key}(adata")

    setup = "    adata = sc.datasets.pbmc68k_reduced()\n"
    assert setup in examples, (setup, examples)
    examples = examples.replace(setup, f"{setup}    sc.tl.{key}(adata)\n")

    p.__doc__ = f"""\
Shortcut for a {name} scatter plot.

See :func:`scanpy.pl.scatter`.

{examples}
"""
    return p


pca = _embedding("pca", "PCA")
umap = _embedding("umap", "UMAP")
tsne = _embedding("tsne", "t-SNE")
diffmap = _embedding("diffmap", "diffusion map")


def heatmap(
    adata: AnnData,
    base: LayerAcc | GraphAcc = A.X,
    /,
    vdims: Collection[AdDim] = (),
    *,
    transpose: bool = False,
    add_dendrogram: bool | Literal["obs", "var"] = False,
) -> hv.HeatMap:
    """Shortcut for a heatmap.

    Basically just

    >>> hv.HeatMap(adata, [A.obs.index, A.var.index], [base[:, :], *vdims]).opts(...)

    Set ``base`` to e.g. ``A`` or ``A.layers[key]``,
    and ``transpose=True`` to switch the order of the dims.

    If ``add_dendrogram`` is True, the dendrogram is added.
    Call it directly to customize the dendrogram:

    >>> hv.operation.dendrogram(heatmap, adjoint_dims=..., main_dim=base[:, :])

    Parameters
    ----------
    adata
        The AnnData object.
    base
        The base layer/graph of the heatmap.
    vdims
        The value dimensions.
    transpose
        Whether to transpose the dims.
    add_dendrogram
        Where to add dendrograms to the heatmap: ``True`` for both,
        ``"obs"``/``"var"`` for one, and ``False`` for none.

    Returns
    -------
    A heatmap object

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        import holoviews as hv

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
        sc.pl.heatmap(
            adata[:, markers], A.X, [A.obs["n_counts"]]
        ).opts(hv.opts.HeatMap(xticks=0, aspect=2))

    With dendrogram:

    ..  holoviews::

        sc.pl.heatmap(
            adata[:, markers], A.X, [A.obs["n_counts"]], add_dendrogram="obs"
        ).opts(hv.opts.HeatMap(xticks=0, aspect=2))

    """
    kdims = (
        [getattr(A, base.dim).index] * 2
        if isinstance(base, GraphAcc)
        else [A.obs.index, A.var.index]
    )
    if transpose:
        kdims.reverse()
    hm = hv.HeatMap(adata, kdims, [base[:, :], *vdims])
    shape = hm.gridded.interface.shape(hm.gridded, gridded=True)
    hm = hm.opts(_supported_opts(hv.HeatMap, show_values=sum(shape) < 80))
    if isinstance(base, GraphAcc):
        hm = hm.opts(aspect="square")
    if add_dendrogram:
        dims = kdims
        if isinstance(add_dendrogram, str):
            dims = [dims[0] if add_dendrogram == "obs" else dims[1]]
        hm = hv.operation.dendrogram(hm, adjoint_dims=dims, main_dim=base[:, :])
    return hm


def tracksplot(
    adata: AnnData,
    /,
    vdims: Collection[AdDim],
    *,
    kdim: AdDim | None = None,
    color: AdDim | None = None,
) -> hv.NdLayout:
    """Tracksplot.

    Parameters
    ----------
    adata
        The AnnData object.
    vdims
        The value dimensions (one per curve).
    kdim
        The key dimension (inferred if ``None``).
    color
        The color dimension.

    Returns
    -------
    A :class:`~holoviews.NdLayout` containing :class:`~holoviews.Curve` objects.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        import holoviews as hv

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
        sc.pl.tracksplot(
            adata, A.X[:, markers], color=A.obs["bulk_labels"]
        ).opts(hv.opts.Curve(aspect=20))

    """
    if kdim is None:
        [dim] = {dim for vdim in vdims for dim in vdim.dims}
        kdim = getattr(A, dim).index
    more_vdims = [] if color is None else [color]
    curves = {
        vdim: hv.Curve(adata, [kdim], [vdim, *more_vdims]).opts(
            xticks=0,
            xlabel="",
            ylabel=vdim.label,
            title="",
            show_legend=False,  # TODO: switch to below impl after fixing https://github.com/holoviz/holoviews/issues/5438
            aspect=2 * len(vdims),
        )
        for vdim in vdims
    }
    if color is not None:
        curves = {vdim: c.groupby(color, hv.NdOverlay) for vdim, c in curves.items()}
    layout = hv.NdLayout(curves, kdims=["marker"]).cols(1)
    return layout.opts(_supported_opts(hv.Curve, frame_width=300))


def _tracksplot2(
    adata: AnnData,
    /,
    vdims: Collection[AdDim],
    *,
    kdim: AdDim | None = None,
    color: AdDim | None = None,
) -> hv.GridSpace:
    """Tracksplot variant. Faster but Gridspace is generally buggy.

    We can switch after <https://github.com/holoviz/holoviews/issues/5438> is fixed.
    """
    if kdim is None:
        [dim] = {dim for vdim in vdims for dim in vdim.dims}
        kdim = getattr(A, dim).index
    more_vdims = [] if color is None else [color]

    _shared = {}

    def link_x(plot, element):
        """Work around <https://github.com/holoviz/holoviews/issues/6948>."""
        fig = plot.state
        if "x_range" not in _shared:
            _shared["x_range"] = fig.x_range
        else:
            fig.x_range = _shared["x_range"]

    return (
        hv
        .GridSpace(
            {
                (0, vdim): hv.Curve(adata, [A.obs.index], [vdim, *more_vdims]).groupby(
                    color, hv.NdOverlay
                )
                for vdim in vdims
            },
            kdims=["_", "marker"],
        )
        .opts(show_legend=True, xaxis=None)
        # .opts(hv.opts.Curve(aspect=2 * len(vdims)))
        .opts(hv.opts.GridSpace(plot_size=(300, 20)))
        .opts(hv.opts.NdOverlay(hooks=[link_x]))
    )


@overload
def violin(
    adata: AnnData,
    /,
    vdims: AdDim,
    *,
    kdims: Collection[AdDim] = (),
    color: AdDim | None = None,
) -> hv.Violin: ...
@overload
def violin(
    adata: AnnData,
    /,
    vdims: AdDim,
    *,
    kdims: Collection[AdDim] = (),
    color: Collection[AdDim],
) -> hv.Layout: ...
@overload
def violin(
    adata: AnnData,
    /,
    vdims: Collection[AdDim],
    *,
    kdims: Collection[AdDim] = (),
    color: AdDim | None = None,
) -> hv.Layout: ...
@overload
def violin(
    adata: AnnData,
    /,
    vdims: Collection[AdDim],
    *,
    kdims: Collection[AdDim] = (),
    color: Collection[AdDim],
) -> hv.NdLayout: ...
def violin(
    adata: AnnData,
    /,
    vdims: Collection[AdDim] | AdDim,
    *,
    kdims: Collection[AdDim] = (),
    color: AdDim | Collection[AdDim] | None = None,
) -> hv.Violin | hv.Layout | hv.NdLayout:
    """Shortcut for a violin plot.

    If ``vdims`` is an ``AdDim``, a single violin is returned:

    >>> hv.Violin(adata, kdims, [vdims, color]).opts(violin_fill_color=color, ...)

    If either ``vdims`` or ``color`` is a collection (not both), a layout is returned:

    >>> hv.Layout([violin(adata, kdims, vdim, ...) for vdim in vdims]).opts(...)

    If both are collections, a 2D :class:`~holoviews.NdLayout` grid is returned,
    one violin per ``(vdim, color)`` combination.

    Parameters
    ----------
    adata
        The AnnData object.
    vdims
        The value dimension(s).
        If a collection is passed, multiple plots are created (see above).
    kdims
        The key dimensions (``color`` will be added automatically).
    color
        The (categorical) color dimension: for each category, a violin is drawn.
        If a collection is passed, multiple plots are created (see above).

    Returns
    -------
    A :class:`~holoviews.Violin` plot or a :class:`~holoviews.Layout` / :class:`~holoviews.NdLayout`
    containing multiple violin plots.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        import holoviews as hv

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.violin(adata, A.obs[["percent_mito", "n_counts", "n_genes"]]).opts(
            hv.opts.Violin(ylim=(0, None))
        )

    ..  holoviews::

        sc.pl.violin(adata, A.obs["S_score"], color=A.obs["bulk_labels"]).opts(
            width=500, xrotation=30
        )

    ..  holoviews::

        sc.pl.violin(
            adata,
            vdims=A.obs[["percent_mito", "n_counts", "n_genes"]],
            color=A.obs[["phase", "louvain"]],
        )

    """
    match vdims, color:
        case AdDim(), AdDim() | None:
            kdims = list(kdims)
            if color and color not in kdims:
                kdims.append(color)
            opts = dict(violin_fill_color=color) if color else {}
            return hv.Violin(adata, kdims, vdims).opts(**opts, ylabel=vdims.label)
        case AdDim(), Collection():
            color = _check_categorical_color(adata, color)
            return _facet(color, lambda c: violin(adata, vdims, kdims=kdims, color=c))
        case Collection(), AdDim() | None:
            vdims = _as_vdims_list(vdims)
            return hv.Layout([
                violin(adata, vdim, color=color).opts(title=vdim.label, ylabel="")
                for vdim in vdims
            ]).opts(axiswise=True)
        case Collection(), Collection():
            vdims = _as_vdims_list(vdims)
            color = _check_categorical_color(adata, color)
            return hv.NdLayout(
                {
                    (vdim.label, c.label): violin(
                        adata, vdim, kdims=kdims, color=c
                    ).opts(title="", ylabel=vdim.label)
                    for vdim in vdims
                    for c in color
                },
                kdims=["vdim", "color"],
            ).cols(len(color))
        case _:
            msg = f"vdims and color must be AdDim or a collection of AdDims, got {vdims!r} and {color!r}."
            raise TypeError(msg)


def stacked_violin(adata: AnnData, /, xdim: AdDim, ydim: AdDim) -> hv.GridSpace:
    """Stacked violin plot.

    Groups data by `xdim` and `ydim` and then plots a single violin for each group.

    Parameters
    ----------
    adata
        The AnnData object.
    xdim
        The x dimension.
    ydim
        The y dimension.

    Returns
    -------
    A :class:`~holoviews.GridSpace` containing :class:`~holoviews.Violin` objects.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc
        import holoviews as hv

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
        sc.pl.stacked_violin(adata[:, markers], A.var.index, A.obs["bulk_labels"])

    """
    if len(xdim.dims) != 1 or len(ydim.dims) != 1:
        msg = "xdim and ydim must map to the same axis."
        raise ValueError(msg)
    xvals = adata[xdim]
    yvals = adata[ydim]

    match next(iter(xdim.dims)), next(iter(ydim.dims)):
        case "obs", "obs":
            idx = lambda x, y: adata[(xvals == x) & (yvals == y), :]  # noqa: E731
        case "var", "var":
            idx = lambda x, y: adata[:, (xvals == x) & (yvals == y)]  # noqa: E731
        case "obs", "var":
            idx = lambda x, y: adata[xvals == x, yvals == y]  # noqa: E731
        case "var", "obs":
            idx = lambda x, y: adata[yvals == y, xvals == x]  # noqa: E731
        case _:
            assert_never()

    return hv.GridSpace(
        {
            # TODO: should Violin vdim be able to be 2D?
            (x, y): hv.Violin(idx(x, y), vdims=[A.X[:, :]]).opts(
                _supported_opts(hv.Violin, inner=None)
            )
            for x in _get_categories(xvals)
            for y in _get_categories(yvals)
        },
        [xdim, ydim],
    )


def dotplot(
    adata: AnnData,
    /,
    group_by: AdDim,
    *,
    funcs: Mapping[str, AggType] = MappingProxyType(
        dict(color="mean", size="count_nonzero")
    ),
) -> hv.Points:
    """Dot plot of marker expression per group.

    For each group/marker combination,
    calculate all `funcs` and map them to a :class:`~holoviews.Points` object’s `opts`.

    By default encodes the mean expression as color,
    and the fraction of expressing cells as dot size.

    Uses `dot_area` with default parameters to calculate size.
    To override, use e.g. `.opts(s=dot_area(hv.dim('median'), ...))`
    (or `size=dot_area(hv.dim('median'), ...) ** 0.5` on Bokeh/Plotly).

    Parameters
    ----------
    adata
        The AnnData object.
    group_by
        The groupby expression.
    funcs
        The aggregation functions.

    Returns
    -------
    A :class:`~holoviews.Points` object.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
        sc.pl.dotplot(adata[:, markers], A.obs["bulk_labels"])

    """
    stats_wide = sc.get.aggregate(adata, group_by, funcs.values())
    stats_long = reduce(
        pd.merge,
        (
            stats_wide
            .to_df(x)
            .reset_index(names="group")
            .melt("group", var_name="marker", value_name=x)
            for x in funcs.values()
        ),
    )

    opts: dict[str, object] = dict(funcs)
    if (d := opts.pop("size", None)) is not None:
        area = dot_area(hv.dim(d))
        opts.update(_supported_opts(hv.Points, s=area, size=area**0.5).options)

    return hv.Points(stats_long, ["group", "marker"], list(funcs.values())).opts(
        xrotation=30, **opts
    )


def matrixplot(
    adata: AnnData,
    /,
    group_by: AdDim,
    *,
    func: AggType = "mean",
    data: LayerAcc | MultiAcc = A.X,
    add_totals: bool = False,
) -> hv.HeatMap | hv.AdjointLayout:
    """Heatmap with totals per column.

    Parameters
    ----------
    adata
        The AnnData object.
    group_by
        The groupby expression.
    func
        The aggregation function.
    data
        The data to plot.
    add_totals
        Whether to add totals per group.

    Returns
    -------
    A heatmap.
    If ``add_totals`` is True, a :class:`~holoviews.AdjointLayout` is returned
    containing the heatmap and a :class:`~holoviews.Bars` object.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
        sc.pl.matrixplot(
            adata[:, markers], A.obs["bulk_labels"], data=A.layers["counts"],
            add_totals=True
        )

    """
    agg = sc.get.aggregate(adata, group_by, func, acc=data)
    agg.var["totals"] = stats.sum(agg.layers[func], axis=0)
    heatmap = hv.HeatMap(
        agg,
        [A.obs.index, A.var.index],
        [A.layers[func][:, :]],
    ).opts(xrotation=30)
    if not add_totals:
        return _add_hover(heatmap)
    bars = hv.Bars(agg, A.var.index, A.var["totals"]).opts(
        yticks=0,
        xlabel="",  # TODO: holoviews issue
    )
    return hv.AdjointLayout([_add_hover(heatmap), _add_hover(bars)])


def _as_vdims_list(vdims: Collection[AdDim], /) -> list[AdDim]:
    vdims = list(vdims)
    if not all(isinstance(vdim, AdDim) for vdim in vdims):
        msg = f"vdims must be an AdDim or a collection of AdDims, got {vdims!r}."
        raise TypeError(msg)
    return vdims


def _check_categorical_color(
    adata: AnnData, color: Collection[AdDim], /
) -> Collection[AdDim]:
    """Validate that all ``color`` dims are categorical (violin groups by them)."""
    non_cat = [c for c in color if not isinstance(adata[c], pd.Categorical)]
    if non_cat:
        msg = (
            "violin's `color` groups violins by category when faceted over a "
            f"collection; got non-categorical dimension(s): {non_cat!r}."
        )
        raise TypeError(msg)
    return color


def _facet[D: hv.core.dimension.Dimensioned](
    dims: Collection[AdDim], plot: Callable[[AdDim], D], /
) -> hv.Layout:
    """Facet a plot over multiple ``color`` dimensions into a `hv.Layout`."""
    return hv.Layout([plot(c).opts(title=c.label) for c in dims]).opts(axiswise=True)


def _get_categories(
    vals: ExtensionArray | np.ndarray,
) -> Iterable[str | int | float]:
    if isinstance(vals, np.ndarray):
        return np.unique(vals)
    if isinstance(vals, pd.Categorical):
        return vals.categories[vals.categories.isin(vals)]
    return vals.unique()


def _add_hover[D: hv.core.dimension.Dimensioned](obj: D) -> D:
    return obj.opts(**_supported_opts(type(obj), tools=["hover"]))


def _supported_opts(
    cls: type[hv.core.dimension.Dimensioned], **opts: object
) -> hv.Options:
    """Filter `opts` down to holoviews plot/style options valid for `cls` on the current backend."""
    plot_cls = hv.Store.registry[hv.Store.current_backend].get(cls)
    if plot_cls is None:
        return hv.Options()
    return getattr(hv.opts, cls.__name__)({
        k: v
        for k, v in opts.items()
        if k in plot_cls.param or k in getattr(plot_cls, "style_opts", ())
    })

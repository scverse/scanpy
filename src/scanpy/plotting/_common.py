"""Shared data preparation utilities for scanpy plotting functions."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from fast_array_utils import stats

from .._utils import own_deprecations
from ..preprocessing._normalization import normalize_total

if TYPE_CHECKING:
    import holoviews as hv
    from anndata import AnnData
    from numpy.typing import NDArray

    from ..get._types import RepAcc, VecRef


__all__ = ["dot_area", "highest_expr_genes"]


def highest_expr_genes(
    adata: AnnData,
    n_top: int,
    *,
    use: RepAcc | str | None = None,
    layer: str | None = None,
    gene_symbols: VecRef | None = None,
) -> pd.DataFrame:
    """Return a wide-format DataFrame of the top ``n_top`` genes by mean expression.

    Normalizes each cell to 100% total counts, then selects the ``n_top`` genes
    with the highest mean fraction.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_top
        Number of top genes to return.
    use
        Matrix to pull data from.
    layer
        Layer to pull data from (for the legacy API).
    gene_symbols
        Gene symbols: an accessor like `A.var['symbol']`, or (for the legacy API)
        a column in ``.var``; uses ``.var_names`` if ``None``.

    Returns
    -------
    DataFrame with cells as rows and top genes as columns, values are
    percent of total counts.
    """
    with own_deprecations():
        norm_expr = normalize_total(
            adata, target_sum=100, use=use, out=True, layer=layer
        ).normalized
    mean_percent = stats.mean(norm_expr, axis=0)
    top_idx = np.argsort(mean_percent)[::-1][:n_top]
    if gene_symbols is None:
        columns = adata.var_names[top_idx]
    else:
        symbols = (
            adata.var[gene_symbols]
            if isinstance(gene_symbols, str)
            else pd.Series(adata[gene_symbols], index=adata.var_names)
        )
        columns = symbols.iloc[top_idx].astype("string")
    counts_top = norm_expr[:, top_idx]
    if hasattr(counts_top, "toarray"):
        counts_top = counts_top.toarray()
    return pd.DataFrame(counts_top, index=adata.obs_names, columns=columns)


def dot_area[V: NDArray[np.floating] | hv.dim](
    vec: V,
    *,
    dot_min: float = 0,
    dot_max: float = 1,
    smallest_dot: float = 0,
    largest_dot: float = 200,
    size_exponent: float = 1.5,
) -> V:
    """Dot area transform, e.g. for :func:`~scanpy.pl.dotplot`.

    Normalizes ``vec`` by its maximum, clips/rescales it between
    ``dot_min``/``dot_max``, raises it to ``size_exponent``,
    and finally scales it into ``[smallest_dot, largest_dot]``.

    ``smallest_dot``/``largest_dot`` are marker *areas*, matching Matplotlib’s ``s``
    (as used by the legacy, Matplotlib-only :func:`~scanpy.pl.dotplot`).
    HoloViews' ``size`` style option means a diameter on Bokeh/Plotly,
    unlike Matplotlib's area-based ``s``,
    so callers targeting those backends need to take the square root of the result.

    Parameters
    ----------
    vec
        Vector or dimension to transform.
    dot_min
        Normalized values below this are clipped to it.
    dot_max
        Normalized values above this are clipped to it.
    smallest_dot
        The dot area mapped to ``dot_min``.
    largest_dot
        The dot area mapped to ``dot_max``.
    size_exponent
        Exponent applied before rescaling to ``[smallest_dot, largest_dot]``.

    Returns
    -------
    Transformed vector or :class:`~holoviews.dim` expression to pass as ``.opts(s=...)``.
    """
    frac = vec / vec.max()
    if dot_min != 0 or dot_max != 1:
        frac = (frac.clip(dot_min, dot_max) - dot_min) / (dot_max - dot_min)
    return frac**size_exponent * (largest_dot - smallest_dot) + smallest_dot

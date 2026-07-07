"""Shared data preparation utilities for scanpy plotting functions."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from fast_array_utils import stats

from ..preprocessing._normalization import normalize_total

if TYPE_CHECKING:
    import holoviews as hv
    from anndata import AnnData
    from numpy.typing import NDArray


__all__ = ["dot_size", "highest_expr_genes"]


def highest_expr_genes(
    adata: AnnData,
    n_top: int,
    *,
    layer: str | None = None,
    gene_symbols: str | None = None,
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
    layer
        Layer to pull data from.
    gene_symbols
        Column in ``.var`` containing gene symbols; uses ``.var_names`` if ``None``.

    Returns
    -------
    DataFrame with cells as rows and top genes as columns, values are
    percent of total counts.
    """
    norm_expr = normalize_total(adata, target_sum=100, layer=layer, inplace=False)["X"]
    mean_percent = stats.mean(norm_expr, axis=0)
    top_idx = np.argsort(mean_percent)[::-1][:n_top]
    columns = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols].iloc[top_idx].astype("string")
    )
    counts_top = norm_expr[:, top_idx]
    if hasattr(counts_top, "toarray"):
        counts_top = counts_top.toarray()
    return pd.DataFrame(counts_top, index=adata.obs_names, columns=columns)


def dot_size[V: NDArray[np.floating] | hv.dim](
    vec: V,
    *,
    dot_min: float = 0,
    dot_max: float = 1,
    smallest_dot: float = 0,
    largest_dot: float = 20,
    size_exponent: float = 1.5,
) -> V:
    """Dot size transform, e.g. for :func:`dotplot`.

    Normalizes ``func`` by its maximum, clips/rescales it between
    ``dot_min``/``dot_max``, raises it to ``size_exponent``,
    and finally scales it into ``[smallest_dot, largest_dot]``.

    Parameters
    ----------
    vec
        Vector or dimension to transform.
    dot_min
        Normalized values below this are clipped to it.
    dot_max
        Normalized values above this are clipped to it.
    smallest_dot
        The dot size mapped to ``dot_min``.
    largest_dot
        The dot size mapped to ``dot_max``.
    size_exponent
        Exponent applied before rescaling to ``[smallest_dot, largest_dot]``.

    Returns
    -------
    A :class:`~holoviews.core.dimension.dim` expression to pass as ``.opts(size=...)``.
    """
    frac = vec / vec.max()
    if dot_min != 0 or dot_max != 1:
        frac = frac.clip(dot_min, dot_max).norm((dot_min, dot_max))
    return frac**size_exponent * (largest_dot - smallest_dot) + smallest_dot

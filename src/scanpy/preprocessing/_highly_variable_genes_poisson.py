"""Analytical Poisson gene selection (``flavor="poisson_gene_selection"``).

Genes are ranked by how strongly their zeros are enriched relative to a
Poisson null model :cite:p:`Andrews2019`. The method is three closed-form
calculations; there is no iterative model fitting::

    p_g     = gene_sum_g / total_counts   (relative abundance)
    p_exp_g = mean_c exp(-p_g * L_c)      (expected zero fraction under Poisson)
    S_g     = p_obs_g * (1 - p_exp_g)     (zero-enrichment probability)

where ``L_c`` is the library size of cell ``c`` and ``p_obs_g`` the observed
zero fraction of gene ``g``.

scvi-tools estimates ``S_g`` by Monte-Carlo sampling: it counts how often a
Bernoulli(``p_obs``) draw exceeds a Bernoulli(``p_exp``) draw. ``S_g`` is the
exact probability of that event, i.e. the value the Monte-Carlo estimate
converges to. rapids-singlecell computes the same closed form on GPU. Output
columns and their semantics match both implementations.

User-facing documentation lives in :func:`~scanpy.pp.highly_variable_genes`.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numba
import numpy as np
import pandas as pd
from fast_array_utils import stats
from fast_array_utils.numba import njit

from .. import logging as logg
from .._compat import warn
from .._utils import (
    axis_nnz,
    check_nonnegative_integers,
    raise_if_dask_feature_axis_chunked,
)
from ..get import _get_arr
from ._distributed import materialize_as_ndarray

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import CSBase, DaskArray

    type _Array = NDArray | CSBase | DaskArray


@njit
def _expected_zero_fraction(
    library_sizes: NDArray[np.float64],
    n_cells_per_size: NDArray[np.float64],
    scaled_mean: NDArray[np.float64],
    n_cells: int,
) -> NDArray[np.float64]:
    """Per gene, the expected fraction of zeros under the Poisson null.

    Computes ``(1 / n_cells) * sum_c exp(-p_g * L_c)``, grouping cells that
    share a library size: ``sum_u n_u * exp(-p_g * L_u)``. This is exact, and
    since raw-count library sizes are integers with many repeats, it reduces
    the work from ``n_genes * n_cells`` to ``n_genes * n_unique_sizes``.
    No ``(genes x cells)`` intermediate is materialized; each gene writes only
    its own output, so the parallel loop is deterministic.
    """
    n_genes = scaled_mean.shape[0]
    out = np.empty(n_genes, dtype=np.float64)
    for g in numba.prange(n_genes):
        p = scaled_mean[g]
        acc = 0.0
        for k in range(library_sizes.shape[0]):
            acc += n_cells_per_size[k] * math.exp(-p * library_sizes[k])
        out[g] = acc / n_cells
    return out


def _zero_enrichment_single_batch(x: _Array) -> dict[str, NDArray[np.float64]]:
    """Zero-enrichment statistics for one raw-count matrix (cells x genes)."""
    library_size, gene_sum, n_cells_expressing = materialize_as_ndarray((
        stats.sum(x, axis=1, dtype=np.float64),
        stats.sum(x, axis=0, dtype=np.float64),
        axis_nnz(x, axis=0),  # ignores explicitly stored zeros
    ))
    n_cells = library_size.shape[0]  # shape of a (materialized) result: dask-safe

    total_counts = gene_sum.sum()
    # all-zero batch: every gene then has p_exp = p_obs = 1 and score 0
    scaled_mean = (
        gene_sum / total_counts if total_counts > 0 else np.zeros_like(gene_sum)
    )

    p_obs = 1.0 - n_cells_expressing / n_cells
    sizes, n_per_size = np.unique(library_size, return_counts=True)
    p_exp = _expected_zero_fraction(
        sizes, n_per_size.astype(np.float64), scaled_mean, n_cells
    )
    return {
        "observed_fraction_zeros": p_obs,
        "expected_fraction_zeros": p_exp,
        "prob_zero_enrichment": p_obs * (1.0 - p_exp),
    }


def _rank_ascending(score: NDArray[np.float64]) -> NDArray[np.int64]:
    """0-based ranks where a higher rank means more enriched.

    Matches scvi-tools and rapids-singlecell (``argsort().argsort()``).
    Ties are broken deterministically in favor of the earlier gene.
    """
    n = score.shape[0]
    rank = np.empty(n, dtype=np.int64)
    rank[np.argsort(-score, kind="stable")] = np.arange(n - 1, -1, -1)
    return rank


def _batch_masks(adata: AnnData, batch_key: str) -> list[NDArray[np.bool]]:
    """One boolean cell mask per non-empty batch, in category order."""
    batch = adata.obs[batch_key]
    if batch.isna().any():
        msg = f"`adata.obs[{batch_key!r}]` contains missing batch labels."
        raise ValueError(msg)
    if not isinstance(batch.dtype, pd.CategoricalDtype):
        batch = batch.astype("category")
    codes = batch.cat.codes.to_numpy()
    return [codes == i for i in np.unique(codes)]  # skips empty categories


def _highly_variable_genes_poisson(
    adata: AnnData,
    *,
    layer: str | None = None,
    n_top_genes: int | None = None,
    batch_key: str | None = None,
    check_values: bool = True,
    subset: bool = False,
    inplace: bool = True,
) -> pd.DataFrame | None:
    """See `highly_variable_genes`.

    Deliberately does not route through `_highly_variable_genes_batched`,
    which combines means and dispersions this flavor does not have.
    """
    if n_top_genes is None:
        n_top_genes = 2000
    if n_top_genes < 1:
        msg = f"`n_top_genes` must be a positive integer, got {n_top_genes}."
        raise ValueError(msg)
    n_top_genes = min(n_top_genes, adata.n_vars)

    x = _get_arr(adata, layer=layer)
    raise_if_dask_feature_axis_chunked(x)
    if check_values and not check_nonnegative_integers(x):
        msg = (
            "`flavor='poisson_gene_selection'` expects raw count data, "
            "but non-integers were found."
        )
        warn(msg, UserWarning)

    masks = [None] if batch_key is None else _batch_masks(adata, batch_key)
    per_batch = [
        _zero_enrichment_single_batch(x if mask is None else x[mask, :])
        for mask in masks
    ]

    # Combine batches exactly like scvi-tools / rapids-singlecell: medians of
    # the statistics and of the within-batch ranks, and the number of batches
    # in which a gene is among that batch's `n_top_genes`.
    def median(key: str) -> NDArray[np.float64]:
        return np.median(np.vstack([b[key] for b in per_batch]), axis=0)

    ranks = np.vstack([_rank_ascending(b["prob_zero_enrichment"]) for b in per_batch])
    median_rank = np.median(ranks, axis=0)
    nbatches = (ranks >= adata.n_vars - n_top_genes).sum(axis=0)

    # sort by nbatches, then median rank (both descending); `lexsort` is
    # stable, so exact ties keep gene order like `DataFrame.nlargest`
    order = np.lexsort((-median_rank, -nbatches))
    highly_variable = np.zeros(adata.n_vars, dtype=bool)
    highly_variable[order[:n_top_genes]] = True

    df = pd.DataFrame(
        {
            "highly_variable": highly_variable,
            "observed_fraction_zeros": median("observed_fraction_zeros"),
            "expected_fraction_zeros": median("expected_fraction_zeros"),
            "prob_zero_enrichment": median("prob_zero_enrichment"),
            "prob_zero_enrichment_rank": median_rank,
        },
        index=adata.var_names,
    )
    if batch_key is not None:
        df["prob_zero_enriched_nbatches"] = nbatches

    if not inplace:
        return df.loc[df["highly_variable"]] if subset else df

    adata.uns["hvg"] = {"flavor": "poisson_gene_selection"}
    hint = (
        "added\n"
        "    'highly_variable', boolean vector (adata.var)\n"
        "    'observed_fraction_zeros', float vector (adata.var)\n"
        "    'expected_fraction_zeros', float vector (adata.var)\n"
        "    'prob_zero_enrichment', float vector (adata.var)\n"
        "    'prob_zero_enrichment_rank', float vector (adata.var)"
    )
    if batch_key is not None:
        hint += "\n    'prob_zero_enriched_nbatches', int vector (adata.var)"
    logg.hint(hint)
    for col in df.columns:
        adata.var[col] = df[col].to_numpy()
    if subset:
        adata._inplace_subset_var(highly_variable)
    return None

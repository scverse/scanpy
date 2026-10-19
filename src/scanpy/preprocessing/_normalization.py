from __future__ import annotations

from operator import truediv
from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from fast_array_utils import stats
from fast_array_utils.numba import njit
from scverse_misc import Deprecation, deprecated_arg

from .. import logging as logg
from .._compat import CSBase, CSCBase, CSRBase, DaskArray, warn
from .._docs import DEPR_KEY_ADDED, doc_out_multi, doc_use
from .._utils import _doc_params, axis_mul_or_truediv, dematrix, view_to_actual
from ..get._out import (
    _dest_key,
    _dest_kind,
    _legacy_out,
    _out_dests,
    _use_dest,
    _write_dests,
    deprecated_copy,
)
from ..get._use import _get_arr, _resolve_obs

if TYPE_CHECKING:
    from collections.abc import Mapping

    from anndata import AnnData

    from ..get._out import OutKind, Prefix
    from ..get._types import ArrAcc, RepAcc, VecRef


class NormalizeTotalOut[
    M = np.ndarray | CSRBase | DaskArray,
    V = np.ndarray | DaskArray,
](NamedTuple):
    """Outputs of :func:`~scanpy.pp.normalize_total`."""

    normalized: M
    """The normalized matrix."""
    norm_factor: V
    """The normalization factor (counts per cell before normalization, divided by `target_sum`)."""


NORMALIZE_TOTAL_KINDS: NormalizeTotalOut[set[OutKind], set[OutKind]] = (
    NormalizeTotalOut({"layer", "obsm"}, {"obs"})
)
"""Where `out` can write the outputs of `pp.normalize_total`."""


def _compute_nnz_median(counts: np.ndarray | DaskArray) -> np.floating:
    """Given a 1D array of counts, compute the median of the non-zero counts."""
    if isinstance(counts, DaskArray):
        counts = counts.compute()
    counts_greater_than_zero = counts[counts > 0]
    median = np.median(counts_greater_than_zero)
    return median


@njit
def _normalize_csr(
    mat: CSRBase,
    *,
    rows,
    columns,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
):
    """For sparse CSR matrix, compute the normalization factors."""
    counts_per_cell = np.zeros(rows, dtype=mat.data.dtype)
    counts_per_cols = np.zeros(columns, dtype=np.int32)
    for i in numba.prange(rows):
        count = 0.0
        for j in range(mat.indptr[i], mat.indptr[i + 1]):
            count += mat.data[j]
        counts_per_cell[i] = count
    if exclude_highly_expressed:
        counts_per_cols_t = np.zeros((n_threads, columns), dtype=np.int32)

        for i in numba.prange(n_threads):
            for r in range(i, rows, n_threads):
                for j in range(mat.indptr[r], mat.indptr[r + 1]):
                    if mat.data[j] > max_fraction * counts_per_cell[r]:
                        minor_index = mat.indices[j]
                        counts_per_cols_t[i, minor_index] += 1
        for c in numba.prange(columns):
            counts_per_cols[c] = counts_per_cols_t[:, c].sum()

        for i in numba.prange(rows):
            count = 0.0
            for j in range(mat.indptr[i], mat.indptr[i + 1]):
                if counts_per_cols[mat.indices[j]] == 0:
                    count += mat.data[j]
            counts_per_cell[i] = count

    return counts_per_cell, counts_per_cols


def _normalize_total_helper(
    x: np.ndarray | CSBase | DaskArray,
    *,
    exclude_highly_expressed: bool,
    max_fraction: float,
    target_sum: float | None,
) -> tuple[np.ndarray | CSBase | DaskArray, np.ndarray, np.ndarray | None]:
    """Calculate the normalized data, counts per cell, and gene subset.

    Parameters
    ----------
    See `normalize_total` for details.

    Returns
    -------
    X
        The normalized data matrix.
    counts_per_cell
        The normalization factors used for each cell (counts / target_sum).
    gene_subset
        If `exclude_highly_expressed=True`, a boolean mask indicating which genes
        were not considered highly expressed. Otherwise, `None`.
    """
    gene_subset = None
    counts_per_cell = None
    if isinstance(x, CSRBase):
        n_threads = numba.get_num_threads()
        counts_per_cell, counts_per_cols = _normalize_csr(
            x,
            rows=x.shape[0],
            columns=x.shape[1],
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            n_threads=n_threads,
        )
        if target_sum is None:
            target_sum = _compute_nnz_median(counts_per_cell)
        if exclude_highly_expressed:
            gene_subset = ~np.where(counts_per_cols)[0]
    else:
        counts_per_cell = stats.sum(x, axis=1)
        if exclude_highly_expressed:
            # at least one cell as more than max_fraction of counts per cell
            hi_exp = dematrix(x > counts_per_cell[:, None] * max_fraction)
            gene_subset = stats.sum(hi_exp, axis=0) == 0

            counts_per_cell = stats.sum(x[:, gene_subset], axis=1)
        if target_sum is None:
            target_sum = _compute_nnz_median(counts_per_cell)

    counts_per_cell = counts_per_cell / target_sum
    out = x if isinstance(x, np.ndarray | CSBase) else None
    x = axis_mul_or_truediv(
        x, counts_per_cell, op=truediv, out=out, allow_divide_by_zero=False, axis=0
    )
    return x, counts_per_cell, gene_subset


@_doc_params(
    use=doc_use("Which matrix to normalize."),
    out=doc_out_multi(
        dict(
            normalized=(
                "The normalized matrix. "
                "By default written to the place `use` reads from, "
                "with `sc.Prefix('p')` to `A.layers['p_normalized']` "
                "(`A.obsm['p_normalized']` if `use` points to `.obsm`)."
            ),
            norm_factor=(
                "The normalization factor of each cell. "
                "By default not written, "
                "with `sc.Prefix('p')` to `A.obs['p_norm_factor']`."
            ),
        ),
        example="{'normalized': A.layers['normalized']}",
        returns=":class:`~scanpy.pp.NormalizeTotalOut`",
        key_added="`key_added='k'` is `out={'normalized': <use>, 'norm_factor': A.obs['k']}`.",
    ),
)
@deprecated_arg("layer", Deprecation("1.13.0", "Use `use`/`out` instead."))
@deprecated_arg("obsm", Deprecation("1.13.0", "Use `use`/`out` instead."))
@deprecated_arg(
    "inplace", Deprecation("1.13.0", "Use `out=True` instead of `inplace=False`.")
)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_copy
def normalize_total(  # noqa: PLR0913
    adata: AnnData,
    *,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    use: RepAcc | str | None = None,
    out: Mapping[str, ArrAcc | VecRef | str] | Prefix | bool = False,
    # deprecated
    key_added: str | None = None,
    layer: str | None = None,
    obsm: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | NormalizeTotalOut | dict[str, np.ndarray] | None:
    """Normalize counts per cell.

    Normalize each cell by total counts over all genes,
    so that every cell has the same total count after normalization.
    If choosing `target_sum=1e6`, this is CPM normalization.

    If `exclude_highly_expressed=True`, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes :cite:p:`Weinreb2017`.

    Similar functions are used, for example, by Seurat :cite:p:`Satija2015`, Cell Ranger
    :cite:p:`Zheng2017` or SPRING :cite:p:`Weinreb2017`.

    .. array-support:: pp.normalize_total

    .. note::
        When used with a :class:`~dask.array.Array` in `adata.X`, this function will have to
        call functions that trigger `.compute()` on the :class:`~dask.array.Array` if `exclude_highly_expressed`
        is `True`, or if `norm_factor` is written.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    target_sum
        If `None`, after normalization, each observation (cell) has a total
        count equal to the median of total counts for observations (cells)
        before normalization.
    exclude_highly_expressed
        Exclude (very) highly expressed genes for the computation of the
        normalization factor (size factor) for each cell. A gene is considered
        highly expressed, if it has more than `max_fraction` of the total counts
        in at least one cell. The not-excluded genes will sum up to
        `target_sum`.  Providing this argument when `adata.X` is a :class:`~dask.array.Array`
        will incur blocking `.compute()` calls on the array.
    max_fraction
        If `exclude_highly_expressed=True`, consider cells as highly expressed
        that have more counts than `max_fraction` of the original total counts
        in at least one cell.
    {use}
    {out}
    inplace
        `inplace=False` is `out=True`.
    copy
        Whether to return a modified copy instead of modifying `adata` in place.

    Returns
    -------
    Returns the outputs as a :class:`~scanpy.pp.NormalizeTotalOut` if `out=True`, else `None`.
    The deprecated `inplace=False` returns them as a :class:`dict` instead.

    Example
    -------
    >>> import sys
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> sc.settings.verbosity = "info"
    >>> sc.settings.logfile = sys.stdout  # for doctests
    >>> np.set_printoptions(precision=2) and None
    >>> adata = AnnData(
    ...     np.array(
    ...         [
    ...             [3, 3, 3, 6, 6],
    ...             [1, 1, 1, 2, 2],
    ...             [1, 22, 1, 2, 2],
    ...         ],
    ...         dtype="float32",
    ...     )
    ... )
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, out=True).normalized
    normalizing counts per cell
        finished (0:00:00)
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(
    ...     adata,
    ...     target_sum=1,
    ...     exclude_highly_expressed=True,
    ...     max_fraction=0.2,
    ...     out=True,
    ... ).normalized
    normalizing counts per cell
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
        finished (0:00:00)
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)

    """
    out = _legacy_out(out, "inplace", None if inplace else True)
    use = _resolve_obs(use)
    # the normalized matrix goes back to where it was read from by default
    src = _use_dest(use, layer=layer, obsm=obsm)
    dests = _out_dests(
        out,
        key_added,
        dests=lambda key: NormalizeTotalOut(
            src, None if key is None else dict(obs=key)
        ),
        kinds=NORMALIZE_TOTAL_KINDS,
        bases=NormalizeTotalOut(
            {_dest_kind(src): "normalized"}, dict(obs="norm_factor")
        ),
    )

    if max_fraction < 0 or max_fraction > 1:
        msg = "Choose max_fraction between 0 and 1."
        raise ValueError(msg)

    view_to_actual(adata)

    x = _get_arr(adata, use, layer=layer, obsm=obsm)
    if isinstance(x, CSCBase):
        x = x.tocsr()
    if dests is True or src != dests[1].normalized:
        # normalization is in-place, so leave the source alone
        # unless we write back over it
        x = x.copy()
    if issubclass(x.dtype.type, int | np.integer):
        x = x.astype(np.float32)  # TODO: Check if float64 should be used

    start = logg.info("normalizing counts per cell")

    x, counts_per_cell, gene_subset = _normalize_total_helper(
        x,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        target_sum=target_sum,
    )

    if exclude_highly_expressed:
        logg.info(
            "The following highly-expressed genes are not considered during normalization factor computation:\n"
            f"{adata.var_names[~np.asarray(gene_subset)].tolist()}"
        )

    cell_subset = counts_per_cell > 0
    if not isinstance(cell_subset, DaskArray) and not np.all(cell_subset):
        warn("Some cells have zero counts", UserWarning)

    vals = NormalizeTotalOut(x, counts_per_cell)
    logg.info("    finished ({time_passed})", time=start)
    if dests is True:
        if not inplace:  # deprecated, returns what it did in 1.12
            return dict(X=vals.normalized, norm_factor=vals.norm_factor)
        return vals
    _write_dests(adata, vals, dests[1])
    if (key := _dest_key(dests[1].norm_factor)) is not None:
        logg.debug(
            f"and added {key!r}, counts per cell before normalization (adata.obs)"
        )
    return None

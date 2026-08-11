"""Simple Preprocessing Functions.

Compositions of these functions are found in sc.preprocess.recipes.
"""

from __future__ import annotations

from contextlib import AbstractContextManager
from dataclasses import dataclass, field
from functools import singledispatch
from typing import TYPE_CHECKING, overload

import numba
import numpy as np
from anndata import AnnData
from fast_array_utils import stats
from fast_array_utils.conv import to_dense
from fast_array_utils.numba import njit
from numpy._typing._array_like import NDArray
from pandas.api.types import CategoricalDtype
from sklearn.utils import check_array

from .. import logging as logg
from .._compat import CSBase, CSRBase, DaskArray
from .._docs import doc_rng
from .._settings import settings
from .._utils import (
    _doc_params,
    _numba_thread_limit,
    _resolve_axis,
    check_array_function_arguments,
    is_backed_type,
    raise_not_implemented_error_if_backed_type,
    sanitize_anndata,
    view_to_actual,
)
from .._utils.random import _accepts_legacy_random_state, _if_legacy_apply_global
from ..get import _check_mask, _get_arr, _set_obs_rep
from ..get._aggregated import Aggregate, _group_counts, _kept_observations
from ._distributed import materialize_as_ndarray

if TYPE_CHECKING:
    from collections.abc import Collection, Sequence
    from numbers import Number
    from typing import Literal, Self

    import pandas as pd
    from numpy.typing import NDArray

    from .._utils.random import RNGLike, SeedLike


def filter_cells(
    data: AnnData | CSBase | np.ndarray | DaskArray,
    *,
    min_counts: int | None = None,
    min_genes: int | None = None,
    max_counts: int | None = None,
    max_genes: int | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | tuple[np.ndarray, np.ndarray] | None:
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters `min_counts`, `min_genes`,
    `max_counts`, `max_genes` per call.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    min_counts
        Minimum number of counts required for a cell to pass filtering.
    min_genes
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts
        Maximum number of counts required for a cell to pass filtering.
    max_genes
        Maximum number of genes expressed required for a cell to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on `inplace`, returns the following arrays or directly subsets
    and annotates the data matrix:

    cells_subset
        Boolean index mask that does filtering. `True` means that the
        cell is kept. `False` means the cell is removed.
    number_per_cell
        Depending on what was thresholded (`counts` or `genes`),
        the array stores `n_counts` or `n_genes` per cell.

    .. array-support:: pp.filter_cells

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.krumsiek11()  # doctest: +ELLIPSIS
    UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
        ...
    >>> adata.obs_names_make_unique()
    >>> adata.n_obs
    640
    >>> adata.var_names.tolist()  # doctest: +NORMALIZE_WHITESPACE
    ['Gata2', 'Gata1', 'Fog1', 'EKLF', 'Fli1', 'SCL',
     'Cebpa', 'Pu.1', 'cJun', 'EgrNab', 'Gfi1']
    >>> # add some true zeros
    >>> adata.X[adata.X < 0.3] = 0
    >>> # simply compute the number of genes per cell
    >>> sc.pp.filter_cells(adata, min_genes=0)
    >>> adata.n_obs
    640
    >>> int(adata.obs["n_genes"].min())
    1
    >>> # filter manually
    >>> adata_copy = adata[adata.obs["n_genes"] >= 3]
    >>> adata_copy.n_obs
    554
    >>> int(adata_copy.obs["n_genes"].min())
    3
    >>> # actually do some filtering
    >>> sc.pp.filter_cells(adata, min_genes=3)
    >>> adata.n_obs
    554
    >>> int(adata.obs["n_genes"].min())
    3

    """
    if copy:
        logg.warning("`copy` is deprecated, use `inplace` instead.")
    n_given_options = sum(
        option is not None for option in [min_genes, min_counts, max_genes, max_counts]
    )
    if n_given_options != 1:
        msg = (
            "Provide exactly one of the optional parameters `min_counts`, "
            "`min_genes`, `max_counts`, `max_genes` per call."
        )
        raise ValueError(msg)
    if isinstance(data, AnnData):
        raise_not_implemented_error_if_backed_type(data.X, "filter_cells")
        adata = data.copy() if copy else data
        cell_subset, number = materialize_as_ndarray(
            filter_cells(
                adata.X,
                min_counts=min_counts,
                min_genes=min_genes,
                max_counts=max_counts,
                max_genes=max_genes,
            ),
        )
        if not inplace:
            return cell_subset, number
        if min_genes is None and max_genes is None:
            adata.obs["n_counts"] = number
        else:
            adata.obs["n_genes"] = number
        adata._inplace_subset_obs(cell_subset)
        return adata if copy else None

    min_number = min_counts if min_genes is None else min_genes
    max_number = max_counts if max_genes is None else max_genes
    number_per_cell = stats.sum(
        data if min_genes is None and max_genes is None else data > 0, axis=1
    )
    if min_number is not None:
        cell_subset = number_per_cell >= min_number
    if max_number is not None:
        cell_subset = number_per_cell <= max_number

    s = stats.sum(~cell_subset)
    if s > 0:
        msg = f"filtered out {s} cells that have "
        if min_genes is not None or min_counts is not None:
            msg += "less than "
            msg += (
                f"{min_genes} genes expressed"
                if min_counts is None
                else f"{min_counts} counts"
            )
        if max_genes is not None or max_counts is not None:
            msg += "more than "
            msg += (
                f"{max_genes} genes expressed"
                if max_counts is None
                else f"{max_counts} counts"
            )
        logg.info(msg)
    return cell_subset, number_per_cell


def filter_genes(
    data: AnnData | CSBase | np.ndarray | DaskArray,
    *,
    min_counts: int | None = None,
    min_cells: int | None = None,
    max_counts: int | None = None,
    max_cells: int | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | tuple[np.ndarray, np.ndarray] | None:
    """Filter genes based on number of cells or counts.

    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells or have at most `max_counts` counts or are expressed
    in at most `max_cells` cells.

    Only provide one of the optional parameters `min_counts`, `min_cells`,
    `max_counts`, `max_cells` per call.

    .. array-support:: pp.filter_genes

    Parameters
    ----------
    data
        An annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts
        Minimum number of counts required for a gene to pass filtering.
    min_cells
        Minimum number of cells expressed required for a gene to pass filtering.
    max_counts
        Maximum number of counts required for a gene to pass filtering.
    max_cells
        Maximum number of cells expressed required for a gene to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on `inplace`, returns the following arrays or directly subsets
    and annotates the data matrix

    gene_subset
        Boolean index mask that does filtering. `True` means that the
        gene is kept. `False` means the gene is removed.
    number_per_gene
        Depending on what was thresholded (`counts` or `cells`), the array stores
        `n_counts` or `n_cells` per gene.

    """
    if copy:
        logg.warning("`copy` is deprecated, use `inplace` instead.")
    n_given_options = sum(
        option is not None for option in [min_cells, min_counts, max_cells, max_counts]
    )
    if n_given_options != 1:
        msg = (
            "Provide exactly one of the optional parameters `min_counts`, "
            "`min_cells`, `max_counts`, `max_cells` per call."
        )
        raise ValueError(msg)

    if isinstance(data, AnnData):
        raise_not_implemented_error_if_backed_type(data.X, "filter_genes")
        adata = data.copy() if copy else data
        gene_subset, number = materialize_as_ndarray(
            filter_genes(
                adata.X,
                min_cells=min_cells,
                min_counts=min_counts,
                max_cells=max_cells,
                max_counts=max_counts,
            )
        )
        if not inplace:
            return gene_subset, number
        if min_cells is None and max_cells is None:
            adata.var["n_counts"] = number
        else:
            adata.var["n_cells"] = number
        adata._inplace_subset_var(gene_subset)
        return adata if copy else None

    min_number = min_counts if min_cells is None else min_cells
    max_number = max_counts if max_cells is None else max_cells
    number_per_gene = stats.sum(
        data if min_cells is None and max_cells is None else data > 0, axis=0
    )
    if min_number is not None:
        gene_subset = number_per_gene >= min_number
    if max_number is not None:
        gene_subset = number_per_gene <= max_number

    s = stats.sum(~gene_subset)
    if s > 0:
        msg = f"filtered out {s} genes that are detected "
        if min_cells is not None or min_counts is not None:
            msg += "in less than "
            msg += (
                f"{min_cells} cells" if min_counts is None else f"{min_counts} counts"
            )
        if max_cells is not None or max_counts is not None:
            msg += "in more than "
            msg += (
                f"{max_cells} cells" if max_counts is None else f"{max_counts} counts"
            )
        logg.info(msg)
    return gene_subset, number_per_gene


@singledispatch
def log1p(
    data: AnnData | np.ndarray | CSBase,
    *,
    base: Number | None = None,
    copy: bool = False,
    chunked: bool | None = None,
    chunk_size: int | None = None,
    layer: str | None = None,
    obsm: str | None = None,
) -> AnnData | np.ndarray | CSBase | None:
    r"""Logarithmize the data matrix.

    Computes :math:`X = \log(X + 1)`,
    where :math:`log` denotes the natural logarithm unless a different base is given.

    .. array-support:: pp.log1p

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    base
        Base of the logarithm. Natural logarithm is used by default.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    chunked
        Process the data matrix in chunks, which will save memory.
        Applies only to :class:`~anndata.AnnData`.
    chunk_size
        `n_obs` of the chunks to process the data in.
    layer
        Entry of layers to transform.
    obsm
        Entry of obsm to transform.

    Returns
    -------
    Returns or updates `data`, depending on `copy`.

    """
    check_array_function_arguments(
        chunked=chunked, chunk_size=chunk_size, layer=layer, obsm=obsm
    )
    return log1p_array(data, copy=copy, base=base)


@log1p.register(CSBase)
def log1p_sparse(x: CSBase, *, base: Number | None = None, copy: bool = False):
    x = check_array(
        x, accept_sparse=("csr", "csc"), dtype=(np.float64, np.float32), copy=copy
    )
    x.data = log1p(x.data, copy=False, base=base)
    return x


@log1p.register(np.ndarray)
def log1p_array(x: np.ndarray, *, base: Number | None = None, copy: bool = False):
    # Can force arrays to be np.ndarrays, but would be useful to not
    # X = check_array(X, dtype=(np.float64, np.float32), ensure_2d=False, copy=copy)
    if copy:
        x = x.astype(float) if not np.issubdtype(x.dtype, np.floating) else x.copy()
    elif not (np.issubdtype(x.dtype, np.floating) or np.issubdtype(x.dtype, complex)):
        x = x.astype(float)
    np.log1p(x, out=x)
    if base is not None:
        np.divide(x, np.log(base), out=x)
    return x


@log1p.register(AnnData)
def log1p_anndata(
    adata: AnnData,
    *,
    base: Number | None = None,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: int | None = None,
    layer: str | None = None,
    obsm: str | None = None,
) -> AnnData | None:
    if "log1p" in adata.uns:
        logg.warning("adata.X seems to be already log-transformed.")

    adata = adata.copy() if copy else adata
    view_to_actual(adata)

    if chunked:
        if (layer is not None) or (obsm is not None):
            msg = (
                "Currently cannot perform chunked operations on arrays not stored in X."
            )
            raise NotImplementedError(msg)
        if adata.isbacked and adata.file._filemode != "r+":
            msg = "log1p is not implemented for backed AnnData with backed mode not r+"
            raise NotImplementedError(msg)
        for chunk, start, end in adata.chunked_X(chunk_size):
            adata.X[start:end] = log1p(chunk, base=base, copy=False)
    else:
        x = _get_arr(adata, layer=layer, obsm=obsm)
        if is_backed_type(x):
            msg = f"log1p is not implemented for matrices of type {type(x)}"
            if layer is not None:
                msg = f"{msg} from layers"
                raise NotImplementedError(msg)
            msg = f"{msg} without `chunked=True`"
            raise NotImplementedError(msg)
        x = log1p(x, copy=False, base=base)
        _set_obs_rep(adata, x, layer=layer, obsm=obsm)

    adata.uns["log1p"] = {"base": base}
    if copy:
        return adata


def sqrt(
    data: AnnData | CSBase | np.ndarray,
    *,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: int | None = None,
) -> AnnData | CSBase | np.ndarray | None:
    r"""Take square root of the data matrix.

    Computes :math:`X = \sqrt(X)`.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    copy
        If an :class:`~anndata.AnnData` object is passed,
        determines whether a copy is returned.
    chunked
        Process the data matrix in chunks, which will save memory.
        Applies only to :class:`~anndata.AnnData`.
    chunk_size
        `n_obs` of the chunks to process the data in.

    Returns
    -------
    Returns or updates `data`, depending on `copy`.

    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        if chunked:
            for chunk, start, end in adata.chunked_X(chunk_size):
                adata.X[start:end] = sqrt(chunk)
        else:
            adata.X = sqrt(data.X)
        return adata if copy else None
    x = data  # proceed with data matrix
    return x.sqrt() if isinstance(x, CSBase) else np.sqrt(x)


@njit
def _create_regressor_categorical(
    x: np.ndarray, /, number_categories: int, cat_array: np.ndarray
) -> np.ndarray:
    # create regressor matrix for categorical variables
    regressors = np.zeros(x.shape, dtype=np.float64)
    # iterate over categories
    for category in range(number_categories):
        # iterate over genes and calculate mean expression
        # for each gene per category
        mask = category == cat_array
        for ix in numba.prange(x.T.shape[0]):
            regressors[mask, ix] = x.T[ix, mask].mean()
    return regressors


@njit
def get_resid(
    data: np.ndarray,
    regressor: np.ndarray,
    coeff: np.ndarray,
) -> np.ndarray:
    for i in numba.prange(data.shape[0]):
        data[i] -= regressor[i] @ coeff
    return data


def numpy_regress_out(
    data: np.ndarray,
    regressor: np.ndarray,
) -> np.ndarray:
    """Regress out unwanted sources of variation, in place.

    Solves via the pseudo-inverse rather than ``inv(RᵀR)``: forming the Gram matrix
    squares the condition number of ``regressor``, which for a barely-varying covariate
    (nearly collinear with the intercept column) destroys the residuals. The
    pseudo-inverse also makes a rank-deficient design well defined, so full-rank,
    ill-conditioned and singular designs all take this one path.

    The covariates are centered first. That leaves the residual unchanged — the column
    space of ``[1, r]`` and of ``[1, r - mean(r)]`` is the same — but it removes the
    offset that otherwise makes a barely-varying covariate ill-conditioned. It does not
    remove :func:`~numpy.linalg.pinv`’s ``rcond`` cutoff: a covariate whose *spread*
    still dwarfs the intercept column by ~1/``rcond`` gets dropped either way.
    """
    regressor = regressor.copy()
    regressor[:, 1:] -= regressor[:, 1:].mean(axis=0)
    coeff = np.linalg.pinv(regressor) @ data
    return get_resid(data, regressor, coeff)


@njit
def _regress_out_per_gene_design(
    x: NDArray[np.floating],
    regressors: NDArray[np.floating],
    is_const: NDArray[np.bool],
) -> NDArray[np.floating]:
    """Residuals of a per-gene ``[1, regressors[:, gene]]`` least-squares fit.

    Only needed when an observation is unassigned: it then gets an all-zero regressor
    and still takes part in the fit, so the reduction in :func:`_regress_out_categorical`
    does not hold and every gene needs its own solve.

    Returns the ``(n_genes, n_obs)`` layout that the caller transposes, so each parallel
    iteration writes one contiguous row.
    """
    n_obs, n_genes = x.shape
    out = np.empty((n_genes, n_obs), dtype=x.dtype)
    for gene in numba.prange(n_genes):
        y = x[:, gene]
        # A gene with no variance needs no regression, and the GLM this replaces
        # refused to fit one – pass it through untouched.
        if is_const[gene]:
            for i in range(n_obs):
                out[gene, i] = y[i]
            continue
        r = regressors[:, gene]
        y_mean = 0.0
        r_mean = 0.0
        for i in range(n_obs):
            y_mean += y[i]
            r_mean += r[i]
        y_mean /= n_obs
        r_mean /= n_obs
        # Centered normal equations: numerically equivalent to the pseudo-inverse
        # solution, without needing LAPACK inside a `prange`.
        s_rr = 0.0
        s_ry = 0.0
        for i in range(n_obs):
            dr = r[i] - r_mean
            s_rr += dr * dr
            s_ry += dr * (y[i] - y_mean)
        if s_rr > 0.0:
            slope = s_ry / s_rr
            intercept = y_mean - slope * r_mean
            for i in range(n_obs):
                out[gene, i] = y[i] - (intercept + slope * r[i])
        else:
            # Constant regressor ⇒ rank-deficient design. The pseudo-inverse
            # solution is the intercept-only fit.
            for i in range(n_obs):
                out[gene, i] = y[i] - y_mean
    return out


def _as_float(x: NDArray | CSBase, *, order: Literal["C", "F"]) -> NDArray | CSBase:
    """Cast integer input to float, as the residuals are written back into `x`’s slot."""
    if np.issubdtype(x.dtype, np.integer):
        target_dtype = np.float32 if x.dtype.itemsize <= 4 else np.float64
        kwargs = {"order": order} if isinstance(x, np.ndarray) else {}
        x = x.astype(target_dtype, **kwargs)
    elif isinstance(x, np.ndarray) and not x.flags.writeable:
        # the residuals are written into `x`, which Numba cannot do to a read-only array
        x = x.copy(order=order)
    return x


def _zero_variance_genes(x: NDArray[np.floating] | CSBase) -> NDArray[np.bool]:
    """Genes with nothing to regress out, which are passed through untouched.

    A single observation makes every gene vacuously constant, which would turn
    `regress_out` into a no-op; two observations already give a residual of zero.
    """
    if x.shape[0] < 2:
        return np.zeros(x.shape[1], dtype=np.bool)
    return stats.is_constant(x, axis=0)


def _regress_out_categorical(
    x: NDArray[np.floating] | CSBase, cats: pd.Categorical
) -> NDArray[np.floating]:
    """Regress each gene on its own per-category mean expression.

    Fitting a gene against its own group means leaves exactly the gene minus that mean:
    the centered normal equations give Σ(rᵢ-r̄)(yᵢ-ȳ) == Σ(rᵢ-r̄)², and r̄ == ȳ. (When r
    is constant the design is rank 1 and the fitted coefficients are not (0, 1), but
    then r ≡ ȳ, so the residual is y - ȳ either way.) So no fit is needed, and the
    means are an (n_categories × n_genes) table rather than one value per cell.
    """
    # a gene with no variance needs no regression; the GLM this replaces refused to
    # fit one, so it is passed through untouched
    is_const = _zero_variance_genes(x)

    if not _kept_observations(cats, None).all():
        # An unassigned observation gets an all-zero regressor yet still takes part in
        # the fit, so the reduction above does not hold – solve each gene separately.
        x = _as_float(x, order="F")
        x = to_dense(x, order="F") if isinstance(x, CSBase) else np.asfortranarray(x)
        codes = cats.codes
        regressors = _create_regressor_categorical(
            x, codes.dtype.type(len(cats.categories)), codes
        )
        return _regress_out_per_gene_design(x, regressors, is_const).T

    # `Aggregate` reads CSR and CSC natively and accumulates in float64. Divide by the
    # counts here rather than via `Aggregate.mean()`: a category with no observations
    # would divide 0/0 and warn, and warnings are errors in this test suite.
    counts = _group_counts(cats, None)
    means = Aggregate(groupby=cats, data=x).sum() / np.maximum(counts, 1)[:, None]
    # nothing writes into `x` on this path, so materialize the output directly
    res = to_dense(x, order="C") if isinstance(x, CSBase) else np.array(x, order="C")
    res = _as_float(res, order="C")
    const_genes = res[:, is_const].copy()
    # subtract per category rather than via `means[codes]`, which would materialize a
    # second float64 copy of the whole matrix
    for cat in range(len(cats.categories)):
        in_cat = cats.codes == cat
        res[in_cat] -= means[cat]
    res[:, is_const] = const_genes
    return res


def regress_out(
    adata: AnnData,
    keys: str | Sequence[str],
    *,
    layer: str | None = None,
    n_jobs: int | None = None,
    copy: bool = False,
) -> AnnData | None:
    """Regress out (mostly) unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R :cite:p:`Satija2015`. Note that this function tends to overcorrect
    in certain circumstances as described in :issue:`526`.

    .. array-support:: pp.regress_out

    Parameters
    ----------
    adata
        The annotated data matrix.
    keys
        Keys for observation annotation on which to regress on.
    layer
        If provided, which element of layers to regress on.
    n_jobs
        Number of Numba threads to use for parallel computation.
        `None` means using :attr:`scanpy.settings.n_jobs`.
        `-1` means using all threads reported by :attr:`numba.config.NUMBA_NUM_THREADS`.
    copy
        Determines whether a copy of `adata` is returned.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an updated `AnnData` object. Sets the following fields:

    `adata.X` | `adata.layers[layer]` : :class:`numpy.ndarray` | :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Corrected count data matrix.

    """
    start = logg.info(f"regressing out {keys}")
    adata = adata.copy() if copy else adata

    sanitize_anndata(adata)

    view_to_actual(adata)

    if isinstance(keys, str):
        keys = [keys]

    x = _get_arr(adata, layer=layer)
    raise_not_implemented_error_if_backed_type(x, "regress_out")
    if isinstance(x, DaskArray):
        # The fit needs every observation at once, so there is nothing to stream.
        # Previously only the rank-deficient path accepted a dask array, and only by
        # materializing it too – now every path does, and says so.
        logg.info("    dask input is computed into memory")
        x = x.compute()

    if isinstance(x, CSBase):
        logg.info("    sparse input is densified and may lead to high memory use")
        if not x.has_canonical_format:
            # An entry stored more than once counts as the sum of its copies, as in
            # `.toarray()`. Canonicalize once so densifying and aggregating agree.
            x = x.copy()
            x.sum_duplicates()

    n_jobs = settings.n_jobs if n_jobs is None else n_jobs

    # every path below is Numba-parallel, so `n_jobs` is simply a thread count
    with _numba_thread_limit(n_jobs) as n_threads:
        logg.debug(f"... using {n_threads} Numba thread(s)")
        res = _regress_out_arr(adata, x, keys)

    _set_obs_rep(adata, res, layer=layer)
    logg.info("    finished", time=start)
    return adata if copy else None


def _regress_out_arr(
    adata: AnnData, x: NDArray[np.floating] | CSBase, keys: Sequence[str]
) -> NDArray[np.floating]:
    """Build the design matrix and return the residuals. Caller sets the thread count."""
    # regress on a single categorical variable
    if keys[0] in adata.obs and isinstance(adata.obs[keys[0]].dtype, CategoricalDtype):
        if len(keys) > 1:
            msg = (
                "If providing categorical variable, "
                "only a single one is allowed. For this one "
                "we regress on the mean for each category."
            )
            raise ValueError(msg)
        logg.debug("... regressing on per-gene means within categories")
        return _regress_out_categorical(x, adata.obs[keys[0]].values)

    # regress on one or several ordinal variables:
    # one design shared by every gene, so solve it once for the whole matrix.
    # `numpy_regress_out` goes through the pseudo-inverse, which handles a
    # rank-deficient design too – no need to branch on rank.
    regressors = adata.obs[keys] if keys else adata.obs.copy()
    # add column of ones at index 0 (first column)
    regressors.insert(0, "ones", 1.0)
    # `bool` and nullable `Int64` columns would otherwise come out as `object`
    regressors = regressors.to_numpy(dtype=np.float64)
    if not np.isfinite(regressors).all():
        msg = f"Cannot regress on `{keys}`: it contains NaN or infinite values."
        raise ValueError(msg)

    x = _as_float(x, order="C")
    x = to_dense(x, order="C") if isinstance(x, CSBase) else x
    # a gene with no variance is left alone, as on the categorical path
    is_const = _zero_variance_genes(x)
    const_genes = x[:, is_const].copy()
    res = numpy_regress_out(x, regressors)
    res[:, is_const] = const_genes
    return res


@overload
def sample(
    data: AnnData,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: Literal[False] = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool] | NDArray[np.floating] | None = None,
) -> None: ...
@overload
def sample(
    data: AnnData,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: Literal[True],
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool] | NDArray[np.floating] | None = None,
) -> AnnData: ...
@overload
def sample[A: np.ndarray | CSBase | DaskArray](
    data: A,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: bool = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool] | NDArray[np.floating] | None = None,
) -> tuple[A, NDArray[np.int64]]: ...


@_doc_params(rng=doc_rng)
def sample(  # noqa: PLR0912
    data: AnnData | np.ndarray | CSBase | DaskArray,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: bool = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool] | NDArray[np.floating] | None = None,
) -> AnnData | tuple[np.ndarray | CSBase | DaskArray, NDArray[np.int64]] | None:
    r"""Sample observations or variables with or without replacement.

    .. array-support:: pp.sample

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    fraction
        Sample to this `fraction` of the number of observations or variables.
        (All of them, even if there are `0`\ s/`False`\ s in `p`.)
        This can be larger than 1.0, if `replace=True`.
        See `axis` and `replace`.
    n
        Sample to this number of observations or variables. See `axis`.
    {rng}
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.
    replace
        If True, samples are drawn with replacement.
    axis
        Sample `obs`\ ervations (axis 0) or `var`\ iables (axis 1).
    p
        Drawing probabilities (floats) or mask (bools).
        Either an `axis`-sized array, or the name of a column.
        If `p` is an array of probabilities, it must sum to 1.

    Returns
    -------
    If `isinstance(data, AnnData)` and `copy=False`,
    this function returns `None`. Otherwise:

    `data[indices, :]` | `data[:, indices]` (depending on `axis`)
        If `data` is array-like or `copy=True`, returns the subset.
    `indices` : numpy.ndarray
        If `data` is array-like, also returns the indices into the original.

    """
    # parameter validation
    if not copy and isinstance(data, AnnData) and data.isbacked:
        msg = "Inplace sampling (`copy=False`) is not implemented for backed objects."
        raise NotImplementedError(msg)
    axis, axis_name = _resolve_axis(axis)
    p = _check_mask(data, p, dim=axis_name, allow_probabilities=True)
    if p is not None and p.dtype == bool:
        p = p.astype(np.float64) / p.sum()
    old_n = data.shape[axis]
    match (fraction, n):
        case (None, None):
            msg = "Either `fraction` or `n` must be set."
            raise TypeError(msg)
        case (None, _):
            pass
        case (_, None):
            if fraction < 0:
                msg = f"`{fraction=}` needs to be nonnegative."
                raise ValueError(msg)
            if not replace and fraction > 1:
                msg = f"If `replace=False`, `{fraction=}` needs to be within [0, 1]."
                raise ValueError(msg)
            n = int(fraction * old_n)
            logg.debug(f"... sampled to {n} {axis_name}")
        case _:
            msg = "Providing both `fraction` and `n` is not allowed."
            raise TypeError(msg)
    del fraction

    # actually do subsampling
    rng = np.random.default_rng(rng)
    indices = rng.choice(old_n, size=n, replace=replace, p=p)

    # overload 1: inplace AnnData subset
    if not copy and isinstance(data, AnnData):
        if axis_name == "obs":
            data._inplace_subset_obs(indices)
        else:
            data._inplace_subset_var(indices)
        return None

    subset = data[indices] if axis_name == "obs" else data[:, indices]

    # overload 2: copy AnnData subset
    if copy and isinstance(data, AnnData):
        assert isinstance(subset, AnnData)
        return subset.to_memory() if data.isbacked else subset.copy()

    # overload 3: return array and indices
    assert isinstance(subset, np.ndarray | CSBase | DaskArray), type(subset)
    if copy:
        subset = subset.copy()
    return subset, indices


@_accepts_legacy_random_state(0)
@_doc_params(rng=doc_rng)
def downsample_counts(
    adata: AnnData,
    counts_per_cell: int | Collection[int] | None = None,
    total_counts: int | None = None,
    *,
    rng: SeedLike | RNGLike | None = None,
    replace: bool = False,
    copy: bool = False,
) -> AnnData | None:
    """Downsample counts from count matrix.

    If `counts_per_cell` is specified, each cell will downsampled.
    If `total_counts` is specified, expression matrix will be downsampled to
    contain at most `total_counts`.

    .. array-support:: pp.downsample_counts

    Parameters
    ----------
    adata
        Annotated data matrix.
    counts_per_cell
        Target total counts per cell. If a cell has more than 'counts_per_cell',
        it will be downsampled to this number. Resulting counts can be specified
        on a per cell basis by passing an array.Should be an integer or integer
        ndarray with same length as number of obs.
    total_counts
        Target total counts. If the count matrix has more than `total_counts`
        it will be downsampled to have this number.
    {rng}
    replace
        Whether to sample the counts with replacement.
    copy
        Determines whether a copy of `adata` is returned.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.X` : :class:`~numpy.ndarray` | :class:`~scipy.sparse.csr_matrix` | :class:`~scipy.sparse.csc_matrix` (dtype `float`)
        Downsampled counts matrix.

    """
    raise_not_implemented_error_if_backed_type(adata.X, "downsample_counts")
    # This logic is all dispatch
    rng = np.random.default_rng(rng)
    if (total_counts is not None) is (counts_per_cell is not None):
        msg = "Must specify exactly one of `total_counts` or `counts_per_cell`."
        raise ValueError(msg)
    if copy:
        adata = adata.copy()
    if total_counts is not None:
        adata.X = _downsample_total_counts(
            adata.X, total_counts, rng=rng, replace=replace
        )
    elif counts_per_cell is not None:
        adata.X = _downsample_per_cell(
            adata.X, counts_per_cell, rng=rng, replace=replace
        )
    return adata if copy else None


def _downsample_per_cell[T: (np.ndarray, CSBase)](
    x: T,
    /,
    counts_per_cell: int | Collection[int],
    *,
    rng: np.random.Generator,
    replace: bool,
) -> T:
    n_obs = x.shape[0]
    counts_per_cell = (
        np.full(n_obs, counts_per_cell)
        if isinstance(counts_per_cell, int)
        else np.asarray(counts_per_cell, np.int_)
    )
    if counts_per_cell.shape != (n_obs,):
        msg = (
            "If provided, 'counts_per_cell' must be either an integer, or "
            "coercible to an `np.ndarray` of length as number of observations"
            " by `np.asarray(counts_per_cell)`."
        )
        raise ValueError(msg)
    with sparse_as_csr(x) as spc:
        x = spc.x  # we only mutate x, so spc.x receives the changes
        rows = np.split(x.data, x.indptr[1:-1]) if isinstance(x, CSRBase) else x
        totals = stats.sum(x, axis=1)
        under_target = np.flatnonzero(totals > counts_per_cell)
        # when parallelizing this, the results will change due to having to spawn `rngs` (e.g. per thread)
        for rowidx in under_target:
            _downsample_array(
                rows[rowidx],
                counts_per_cell[rowidx],
                rng=rng,
                replace=replace,
                inplace=True,
            )
    return spc.x  # use x that was converted back


def _downsample_total_counts[T: (np.ndarray, CSBase)](
    x: T,
    /,
    total_counts: int,
    *,
    rng: np.random.Generator,
    replace: bool,
) -> T:
    total_counts = int(total_counts)
    total = x.sum()
    if total < total_counts:
        return x
    with sparse_as_csr(x) as spc:
        x = spc.x  # we only mutate x, so spc.x receives the changes
        v = x.data if isinstance(x, CSBase) else x.reshape(-1)
        _downsample_array(v, total_counts, rng=rng, replace=replace, inplace=True)
    return spc.x  # use x that was converted back


def _downsample_array(
    col: np.ndarray,
    target: int,
    *,
    rng: np.random.Generator,
    replace: bool = True,
    inplace: bool = False,
) -> np.ndarray:
    """Evenly reduce counts in cell to target amount.

    This is an internal function and has some restrictions:

    * total counts in cell must be less than target
    """
    rng = _if_legacy_apply_global(rng)
    cumcounts = col.cumsum()
    total = np.int_(cumcounts[-1])
    sample = rng.choice(total, target, replace=replace)
    sample.sort()
    return _downsample_array_inner(col, cumcounts, sample, inplace=inplace)


# TODO: can/should this be parallelized?
@numba.njit(cache=True)  # noqa: TID251
def _downsample_array_inner(
    col: np.ndarray, cumcounts: np.ndarray, sample: np.ndarray, *, inplace: bool
) -> np.ndarray:
    if inplace:
        col[:] = 0
    else:
        col = np.zeros_like(col)
    geneptr = 0
    for count in sample:
        while count >= cumcounts[geneptr]:
            geneptr += 1
        col[geneptr] += 1
    return col


@dataclass
class sparse_as_csr[T: (np.ndarray | CSBase)](AbstractContextManager):  # noqa: N801
    """Context manager that converts to CSR while active."""

    x: T
    _original_type: type[T] = field(init=False)

    def __post_init__(self) -> None:
        self._original_type = type(self.x)
        if isinstance(self.x, CSBase) and not isinstance(self.x, CSRBase):
            self.x = self.x.tocsr()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *args: object) -> bool | None:
        if isinstance(self.x, CSBase):
            self.x.eliminate_zeros()
            if not issubclass(self._original_type, CSRBase):
                self.x = self._original_type(self.x)

"""Simple Preprocessing Functions.

Compositions of these functions are found in sc.preprocess.recipes.
"""

from __future__ import annotations

import warnings
from functools import singledispatch
from itertools import repeat
from typing import TYPE_CHECKING, TypeVar, overload

import numba
import numpy as np
from anndata import AnnData
from pandas.api.types import CategoricalDtype
from scipy import sparse
from sklearn.utils import check_array, sparsefuncs

from .. import logging as logg
from .._compat import (
    CSBase,
    CSRBase,
    DaskArray,
    _register_union,
    deprecated,
    njit,
    old_positionals,
)
from .._settings import settings as sett
from .._utils import (
    _check_array_function_arguments,
    _resolve_axis,
    axis_sum,
    is_backed_type,
    raise_not_implemented_error_if_backed_type,
    renamed_arg,
    sanitize_anndata,
    view_to_actual,
)
from ..get import _check_mask, _get_obs_rep, _set_obs_rep
from ._distributed import materialize_as_ndarray
from ._utils import _to_dense

try:
    import dask.array as da
except ImportError:
    da = None

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable, Sequence
    from numbers import Number
    from typing import Literal

    import pandas as pd
    from numpy.typing import NDArray

    from .._utils.random import RNGLike, SeedLike, _LegacyRandom


A = TypeVar("A", bound=np.ndarray | CSBase | DaskArray)


@old_positionals(
    "min_counts", "min_genes", "max_counts", "max_genes", "inplace", "copy"
)
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
        the array stores `n_counts` or `n_cells` per gene.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.krumsiek11()
    UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
        utils.warn_names_duplicates("obs")
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
            "Only provide one of the optional parameters `min_counts`, "
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
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    max_number = max_counts if max_genes is None else max_genes
    number_per_cell = axis_sum(
        X if min_genes is None and max_genes is None else X > 0, axis=1
    )
    if isinstance(number_per_cell, np.matrix):
        number_per_cell = number_per_cell.A1
    if min_number is not None:
        cell_subset = number_per_cell >= min_number
    if max_number is not None:
        cell_subset = number_per_cell <= max_number

    s = axis_sum(~cell_subset)
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


@old_positionals(
    "min_counts", "min_cells", "max_counts", "max_cells", "inplace", "copy"
)
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
            "Only provide one of the optional parameters `min_counts`, "
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

    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_cells is None else min_cells
    max_number = max_counts if max_cells is None else max_cells
    number_per_gene = axis_sum(
        X if min_cells is None and max_cells is None else X > 0, axis=0
    )
    if isinstance(number_per_gene, np.matrix):
        number_per_gene = number_per_gene.A1
    if min_number is not None:
        gene_subset = number_per_gene >= min_number
    if max_number is not None:
        gene_subset = number_per_gene <= max_number

    s = axis_sum(~gene_subset)
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


@renamed_arg("X", "data", pos_0=True)
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
    _check_array_function_arguments(
        chunked=chunked, chunk_size=chunk_size, layer=layer, obsm=obsm
    )
    return log1p_array(data, copy=copy, base=base)


@_register_union(log1p, CSBase)
def log1p_sparse(X: CSBase, *, base: Number | None = None, copy: bool = False):
    X = check_array(
        X, accept_sparse=("csr", "csc"), dtype=(np.float64, np.float32), copy=copy
    )
    X.data = log1p(X.data, copy=False, base=base)
    return X


@log1p.register(np.ndarray)
def log1p_array(X: np.ndarray, *, base: Number | None = None, copy: bool = False):
    # Can force arrays to be np.ndarrays, but would be useful to not
    # X = check_array(X, dtype=(np.float64, np.float32), ensure_2d=False, copy=copy)
    if copy:
        X = X.astype(float) if not np.issubdtype(X.dtype, np.floating) else X.copy()
    elif not (np.issubdtype(X.dtype, np.floating) or np.issubdtype(X.dtype, complex)):
        X = X.astype(float)
    np.log1p(X, out=X)
    if base is not None:
        np.divide(X, np.log(base), out=X)
    return X


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
        X = _get_obs_rep(adata, layer=layer, obsm=obsm)
        if is_backed_type(X):
            msg = f"log1p is not implemented for matrices of type {type(X)}"
            if layer is not None:
                msg = f"{msg} from layers"
                raise NotImplementedError(msg)
            msg = f"{msg} without `chunked=True`"
            raise NotImplementedError(msg)
        X = log1p(X, copy=False, base=base)
        _set_obs_rep(adata, X, layer=layer, obsm=obsm)

    adata.uns["log1p"] = {"base": base}
    if copy:
        return adata


@old_positionals("copy", "chunked", "chunk_size")
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
    X = data  # proceed with data matrix
    return X.sqrt() if isinstance(X, CSBase) else np.sqrt(X)


@deprecated("Use sc.pp.normalize_total instead")
@old_positionals(
    "counts_per_cell_after",
    "counts_per_cell",
    "key_n_counts",
    "copy",
    "layers",
    "use_rep",
    "min_counts",
)
def normalize_per_cell(
    data: AnnData | np.ndarray | CSBase,
    *,
    counts_per_cell_after: float | None = None,
    counts_per_cell: np.ndarray | None = None,
    key_n_counts: str = "n_counts",
    copy: bool = False,
    layers: Literal["all"] | Iterable[str] = (),
    use_rep: Literal["after", "X"] | None = None,
    min_counts: int = 1,
) -> AnnData | np.ndarray | CSBase | None:
    """Normalize total counts per cell.

    .. deprecated:: 1.3.7

       Use :func:`~scanpy.pp.normalize_total` instead.
       The new function is equivalent to the present
       function, except that

       * the new function doesn't filter cells based on `min_counts`,
         use :func:`~scanpy.pp.filter_cells` if filtering is needed.
       * some arguments were renamed
       * `copy` is replaced by `inplace`

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat :cite:p:`Satija2015`, Cell Ranger
    :cite:p:`Zheng2017` or SPRING :cite:p:`Weinreb2017`.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell
        Precomputed counts per cell.
    key_n_counts
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    min_counts
        Cells with counts less than `min_counts` are filtered out during
        normalization.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an updated `AnnData` object. Sets the following fields:

    `adata.X` : :class:`numpy.ndarray` | :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Normalized count data matrix.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = AnnData(np.array([[1, 0], [3, 0], [5, 6]], dtype=np.float32))
    >>> print(adata.X.sum(axis=1))
    [ 1.  3. 11.]
    >>> sc.pp.normalize_per_cell(adata)
    >>> print(adata.obs)
       n_counts
    0       1.0
    1       3.0
    2      11.0
    >>> print(adata.X.sum(axis=1))
    [3. 3. 3.]
    >>> sc.pp.normalize_per_cell(
    ...     adata,
    ...     counts_per_cell_after=1,
    ...     key_n_counts="n_counts2",
    ... )
    >>> print(adata.obs)
       n_counts  n_counts2
    0       1.0        3.0
    1       3.0        3.0
    2      11.0        3.0
    >>> print(adata.X.sum(axis=1))
    [1. 1. 1.]

    """
    if isinstance(data, AnnData):
        start = logg.info("normalizing by total count per cell")
        adata = data.copy() if copy else data
        if counts_per_cell is None:
            cell_subset, counts_per_cell = materialize_as_ndarray(
                filter_cells(adata.X, min_counts=min_counts)
            )
            adata.obs[key_n_counts] = counts_per_cell
            adata._inplace_subset_obs(cell_subset)
            counts_per_cell = counts_per_cell[cell_subset]
        normalize_per_cell(
            adata.X,
            counts_per_cell_after=counts_per_cell_after,
            counts_per_cell=counts_per_cell,
        )

        layers = adata.layers.keys() if layers == "all" else layers
        if use_rep == "after":
            after = counts_per_cell_after
        elif use_rep == "X":
            after = np.median(counts_per_cell[cell_subset])
        elif use_rep is None:
            after = None
        else:
            msg = 'use_rep should be "after", "X" or None'
            raise ValueError(msg)
        for layer in layers:
            _subset, counts = filter_cells(adata.layers[layer], min_counts=min_counts)
            temp = normalize_per_cell(adata.layers[layer], after, counts, copy=True)
            adata.layers[layer] = temp

        logg.info(
            "    finished ({time_passed}): normalized adata.X and added\n"
            f"    {key_n_counts!r}, counts per cell before normalization (adata.obs)",
            time=start,
        )
        return adata if copy else None
    # proceed with data matrix
    X = data.copy() if copy else data
    if counts_per_cell is None:
        if not copy:
            msg = "Can only be run with copy=True"
            raise ValueError(msg)
        cell_subset, counts_per_cell = filter_cells(X, min_counts=min_counts)
        X = X[cell_subset]
        counts_per_cell = counts_per_cell[cell_subset]
    if counts_per_cell_after is None:
        counts_per_cell_after = np.median(counts_per_cell)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        counts_per_cell += counts_per_cell == 0
        counts_per_cell /= counts_per_cell_after
        if not isinstance(X, CSBase):
            X /= counts_per_cell[:, np.newaxis]
        else:
            sparsefuncs.inplace_row_scale(X, 1 / counts_per_cell)
    return X if copy else None


DT = TypeVar("DT")


@njit
def _create_regressor_categorical(
    X: np.ndarray, number_categories: int, cat_array: np.ndarray
) -> np.ndarray:
    # create regressor matrix for categorical variables
    # would be best to use X dtype but this matches old behavior
    regressors = np.zeros(X.shape, dtype=np.float32)
    # iterate over categories
    for category in range(number_categories):
        # iterate over genes and calculate mean expression
        # for each gene per category
        mask = category == cat_array
        for ix in numba.prange(X.T.shape[0]):
            regressors[mask, ix] = X.T[ix, mask].mean()
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
    """Numba kernel for regress out unwanted sorces of variantion.

    Finding coefficient using Linear regression (Linear Least Squares).
    """
    inv_gram_matrix = np.linalg.inv(regressor.T @ regressor)
    coeff = inv_gram_matrix @ (regressor.T @ data)
    data = get_resid(data, regressor, coeff)
    return data


@old_positionals("layer", "n_jobs", "copy")
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

    Parameters
    ----------
    adata
        The annotated data matrix.
    keys
        Keys for observation annotation on which to regress on.
    layer
        If provided, which element of layers to regress on.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy.settings.n_jobs`.
    copy
        Determines whether a copy of `adata` is returned.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an updated `AnnData` object. Sets the following fields:

    `adata.X` | `adata.layers[layer]` : :class:`numpy.ndarray` | :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Corrected count data matrix.

    """
    from joblib import Parallel, delayed

    start = logg.info(f"regressing out {keys}")
    adata = adata.copy() if copy else adata

    sanitize_anndata(adata)

    view_to_actual(adata)

    if isinstance(keys, str):
        keys = [keys]

    X = _get_obs_rep(adata, layer=layer)
    raise_not_implemented_error_if_backed_type(X, "regress_out")

    if isinstance(X, CSBase):
        logg.info("    sparse input is densified and may lead to high memory use")

    n_jobs = sett.n_jobs if n_jobs is None else n_jobs

    # regress on a single categorical variable
    variable_is_categorical = False
    if keys[0] in adata.obs_keys() and isinstance(
        adata.obs[keys[0]].dtype, CategoricalDtype
    ):
        if len(keys) > 1:
            msg = (
                "If providing categorical variable, "
                "only a single one is allowed. For this one "
                "we regress on the mean for each category."
            )
            raise ValueError(msg)
        logg.debug("... regressing on per-gene means within categories")
        # set number of categories to the same dtype as the categories
        cat_array = adata.obs[keys[0]].cat.codes.to_numpy()
        number_categories = cat_array.dtype.type(len(adata.obs[keys[0]].cat.categories))

        X = _to_dense(X, order="F") if isinstance(X, CSBase) else X
        if np.issubdtype(X.dtype, np.integer):
            target_dtype = np.float32 if X.dtype.itemsize <= 4 else np.float64
            X = X.astype(target_dtype)
        regressors = _create_regressor_categorical(X, number_categories, cat_array)
        variable_is_categorical = True
    # regress on one or several ordinal variables
    else:
        # create data frame with selected keys (if given)
        regressors = adata.obs[keys] if keys else adata.obs.copy()

        # add column of ones at index 0 (first column)
        regressors.insert(0, "ones", 1.0)
        regressors = regressors.to_numpy()

    # if the regressors are not categorical and the matrix is not singular
    # use the shortcut numpy_regress_out
    if not variable_is_categorical and np.linalg.det(regressors.T @ regressors) != 0:
        # Because we update `X` in `numpy_regress_out`, it needs to be floating point to match
        # the incoming values.
        if np.issubdtype(X.dtype, np.integer):
            target_dtype = np.float32 if X.dtype.itemsize <= 4 else np.float64
            kwargs = {}
            if isinstance(X, np.ndarray):
                kwargs["order"] = "C"
            X = X.astype(target_dtype, **kwargs)
        X = _to_dense(X, order="C") if isinstance(X, CSBase) else X
        res = numpy_regress_out(X, regressors)

    # for a categorical variable or if the above checks failed,
    # we fall back to the GLM implemetation of regression.
    else:
        # split the adata.X matrix by columns in chunks of size n_chunk
        # (the last chunk could be of smaller size than the others)
        len_chunk = int(np.ceil(min(1000, X.shape[1]) / n_jobs))
        n_chunks = int(np.ceil(X.shape[1] / len_chunk))
        X = _to_dense(X, order="F") if isinstance(X, CSBase) else X
        chunk_list = np.array_split(X, n_chunks, axis=1)
        regressors_chunk = (
            np.array_split(regressors, n_chunks, axis=1)
            if variable_is_categorical
            else repeat(regressors)
        )

        # each task is passed a data chunk (e.g. `adata.X[:, 0:100]``) and the regressors.
        # This data will be passed to each of the jobs.
        # TODO: figure out how to test that this doesn't oversubscribe resources
        res = Parallel(n_jobs=n_jobs)(
            delayed(_regress_out_chunk)(
                data_chunk, regres, variable_is_categorical=variable_is_categorical
            )
            for data_chunk, regres in zip(chunk_list, regressors_chunk, strict=False)
        )

        # res is a list of vectors (each corresponding to a regressed gene column).
        # The transpose is needed to get the matrix in the shape needed
        res = np.vstack(res).T

    _set_obs_rep(adata, res, layer=layer)
    logg.info("    finished", time=start)
    return adata if copy else None


def _regress_out_chunk(
    data_chunk: NDArray[np.floating],
    regressors: pd.DataFrame | NDArray[np.floating],
    *,
    variable_is_categorical: bool,
) -> NDArray[np.floating]:
    import statsmodels.api as sm
    import statsmodels.tools.sm_exceptions as sme

    responses_chunk_list = []
    for col_index in range(data_chunk.shape[1]):
        # if all values are identical, the statsmodel.api.GLM throws an error;
        # but then no regression is necessary anyways...
        if not (data_chunk[:, col_index] != data_chunk[0, col_index]).any():
            responses_chunk_list.append(data_chunk[:, col_index])
            continue

        if variable_is_categorical:
            regres = np.c_[np.ones(regressors.shape[0]), regressors[:, col_index]]
        else:
            regres = regressors

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", sme.PerfectSeparationWarning)
                result = sm.GLM(
                    data_chunk[:, col_index], regres, family=sm.families.Gaussian()
                ).fit()
                new_column = result.resid_response
        except (sme.PerfectSeparationError, sme.PerfectSeparationWarning):
            logg.warning("Encountered perfect separation, setting to 0 as in R.")
            new_column = np.zeros(data_chunk.shape[0])

        responses_chunk_list.append(new_column)

    return np.vstack(responses_chunk_list)


@overload
def sample(
    data: AnnData,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = 0,
    copy: Literal[False] = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool_] | NDArray[np.floating] | None = None,
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
    p: str | NDArray[np.bool_] | NDArray[np.floating] | None = None,
) -> AnnData: ...
@overload
def sample(
    data: A,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: bool = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool_] | NDArray[np.floating] | None = None,
) -> tuple[A, NDArray[np.int64]]: ...
def sample(  # noqa: PLR0912
    data: AnnData | np.ndarray | CSBase | DaskArray,
    fraction: float | None = None,
    *,
    n: int | None = None,
    rng: RNGLike | SeedLike | None = None,
    copy: bool = False,
    replace: bool = False,
    axis: Literal["obs", 0, "var", 1] = "obs",
    p: str | NDArray[np.bool_] | NDArray[np.floating] | None = None,
) -> AnnData | None | tuple[np.ndarray | CSBase | DaskArray, NDArray[np.int64]]:
    r"""Sample observations or variables with or without replacement.

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
    rng
        Random seed to change subsampling.
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


@renamed_arg("target_counts", "counts_per_cell")
def downsample_counts(
    adata: AnnData,
    counts_per_cell: int | Collection[int] | None = None,
    total_counts: int | None = None,
    *,
    random_state: _LegacyRandom = 0,
    replace: bool = False,
    copy: bool = False,
) -> AnnData | None:
    """Downsample counts from count matrix.

    If `counts_per_cell` is specified, each cell will downsampled.
    If `total_counts` is specified, expression matrix will be downsampled to
    contain at most `total_counts`.

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
    random_state
        Random seed for subsampling.
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
    total_counts_call = total_counts is not None
    counts_per_cell_call = counts_per_cell is not None
    if total_counts_call is counts_per_cell_call:
        msg = "Must specify exactly one of `total_counts` or `counts_per_cell`."
        raise ValueError(msg)
    if copy:
        adata = adata.copy()
    if total_counts_call:
        adata.X = _downsample_total_counts(
            adata.X, total_counts, random_state=random_state, replace=replace
        )
    elif counts_per_cell_call:
        adata.X = _downsample_per_cell(
            adata.X, counts_per_cell, random_state=random_state, replace=replace
        )
    if copy:
        return adata


def _downsample_per_cell(
    X: CSBase,
    counts_per_cell: int,
    *,
    random_state: _LegacyRandom,
    replace: bool,
) -> CSBase:
    n_obs = X.shape[0]
    if isinstance(counts_per_cell, int):
        counts_per_cell = np.full(n_obs, counts_per_cell)
    else:
        counts_per_cell = np.asarray(counts_per_cell)
    # np.random.choice needs int arguments in numba code:
    counts_per_cell = counts_per_cell.astype(np.int_, copy=False)
    if not isinstance(counts_per_cell, np.ndarray) or len(counts_per_cell) != n_obs:
        msg = (
            "If provided, 'counts_per_cell' must be either an integer, or "
            "coercible to an `np.ndarray` of length as number of observations"
            " by `np.asarray(counts_per_cell)`."
        )
        raise ValueError(msg)
    if isinstance(X, CSBase):
        original_type = type(X)
        if not isinstance(X, CSRBase):
            X = sparse.csr_matrix(X)  # noqa: TID251
        totals = np.ravel(axis_sum(X, axis=1))  # Faster for csr matrix
        under_target = np.nonzero(totals > counts_per_cell)[0]
        rows = np.split(X.data, X.indptr[1:-1])
        for rowidx in under_target:
            row = rows[rowidx]
            _downsample_array(
                row,
                counts_per_cell[rowidx],
                random_state=random_state,
                replace=replace,
                inplace=True,
            )
        X.eliminate_zeros()
        if not issubclass(original_type, CSRBase):  # Put it back
            X = original_type(X)
    else:
        totals = np.ravel(axis_sum(X, axis=1))
        under_target = np.nonzero(totals > counts_per_cell)[0]
        for rowidx in under_target:
            row = X[rowidx, :]
            _downsample_array(
                row,
                counts_per_cell[rowidx],
                random_state=random_state,
                replace=replace,
                inplace=True,
            )
    return X


def _downsample_total_counts(
    X: CSBase,
    total_counts: int,
    *,
    random_state: _LegacyRandom,
    replace: bool,
) -> CSBase:
    total_counts = int(total_counts)
    total = X.sum()
    if total < total_counts:
        return X
    if isinstance(X, CSBase):
        original_type = type(X)
        if not isinstance(X, CSRBase):
            X = sparse.csr_matrix(X)  # noqa: TID251
        _downsample_array(
            X.data,
            total_counts,
            random_state=random_state,
            replace=replace,
            inplace=True,
        )
        X.eliminate_zeros()
        if not issubclass(original_type, CSRBase):
            X = original_type(X)
    else:
        v = X.reshape(np.multiply(*X.shape))
        _downsample_array(
            v, total_counts, random_state=random_state, replace=replace, inplace=True
        )
    return X


# TODO: can/should this be parallelized?
@numba.njit(cache=True)  # noqa: TID251
def _downsample_array(
    col: np.ndarray,
    target: int,
    *,
    random_state: _LegacyRandom = 0,
    replace: bool = True,
    inplace: bool = False,
):
    """Evenly reduce counts in cell to target amount.

    This is an internal function and has some restrictions:

    * total counts in cell must be less than target
    """
    np.random.seed(random_state)
    cumcounts = col.cumsum()
    if inplace:
        col[:] = 0
    else:
        col = np.zeros_like(col)
    total = np.int_(cumcounts[-1])
    sample = np.random.choice(total, target, replace=replace)
    sample.sort()
    geneptr = 0
    for count in sample:
        while count >= cumcounts[geneptr]:
            geneptr += 1
        col[geneptr] += 1
    return col

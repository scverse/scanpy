"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""
from functools import singledispatch
from numbers import Number
import warnings
from typing import Union, Optional, Tuple, Collection, Sequence, Iterable

import numba
import numpy as np
import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from sklearn.utils import sparsefuncs, check_array
from pandas.api.types import is_categorical_dtype
from anndata import AnnData

from .. import logging as logg
from .._settings import settings as sett
from .._utils import (
    sanitize_anndata,
    deprecated_arg_names,
    view_to_actual,
    AnyRandom,
    _check_array_function_arguments,
)
from .._compat import Literal
from ..get import _get_obs_rep, _set_obs_rep
from ._distributed import materialize_as_ndarray
from ._utils import _get_mean_var

# install dask if available
try:
    import dask.array as da
except ImportError:
    da = None

# backwards compat
from ._deprecated.highly_variable_genes import filter_genes_dispersion


def filter_cells(
    data: AnnData,
    min_counts: Optional[int] = None,
    min_genes: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_genes: Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """\
    Filter cell outliers based on counts and numbers of genes expressed.

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
        Depending on what was tresholded (`counts` or `genes`),
        the array stores `n_counts` or `n_cells` per gene.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.krumsiek11()
    >>> adata.n_obs
    640
    >>> adata.var_names
    ['Gata2' 'Gata1' 'Fog1' 'EKLF' 'Fli1' 'SCL' 'Cebpa'
     'Pu.1' 'cJun' 'EgrNab' 'Gfi1']
    >>> # add some true zeros
    >>> adata.X[adata.X < 0.3] = 0
    >>> # simply compute the number of genes per cell
    >>> sc.pp.filter_cells(adata, min_genes=0)
    >>> adata.n_obs
    640
    >>> adata.obs['n_genes'].min()
    1
    >>> # filter manually
    >>> adata_copy = adata[adata.obs['n_genes'] >= 3]
    >>> adata_copy.obs['n_genes'].min()
    >>> adata.n_obs
    554
    >>> adata.obs['n_genes'].min()
    3
    >>> # actually do some filtering
    >>> sc.pp.filter_cells(adata, min_genes=3)
    >>> adata.n_obs
    554
    >>> adata.obs['n_genes'].min()
    3
    """
    if copy:
        logg.warning('`copy` is deprecated, use `inplace` instead.')
    n_given_options = sum(
        option is not None for option in [min_genes, min_counts, max_genes, max_counts]
    )
    if n_given_options != 1:
        raise ValueError(
            'Only provide one of the optional parameters `min_counts`, '
            '`min_genes`, `max_counts`, `max_genes` per call.'
        )
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        cell_subset, number = materialize_as_ndarray(
            filter_cells(adata.X, min_counts, min_genes, max_counts, max_genes)
        )
        if not inplace:
            return cell_subset, number
        if min_genes is None and max_genes is None:
            adata.obs['n_counts'] = number
        else:
            adata.obs['n_genes'] = number
        adata._inplace_subset_obs(cell_subset)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    max_number = max_counts if max_genes is None else max_genes
    number_per_cell = np.sum(
        X if min_genes is None and max_genes is None else X > 0, axis=1
    )
    if issparse(X):
        number_per_cell = number_per_cell.A1
    if min_number is not None:
        cell_subset = number_per_cell >= min_number
    if max_number is not None:
        cell_subset = number_per_cell <= max_number

    s = np.sum(~cell_subset)
    if s > 0:
        msg = f'filtered out {s} cells that have '
        if min_genes is not None or min_counts is not None:
            msg += 'less than '
            msg += (
                f'{min_genes} genes expressed'
                if min_counts is None
                else f'{min_counts} counts'
            )
        if max_genes is not None or max_counts is not None:
            msg += 'more than '
            msg += (
                f'{max_genes} genes expressed'
                if max_counts is None
                else f'{max_counts} counts'
            )
        logg.info(msg)
    return cell_subset, number_per_cell


def filter_genes(
    data: AnnData,
    min_counts: Optional[int] = None,
    min_cells: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_cells: Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Union[AnnData, None, Tuple[np.ndarray, np.ndarray]]:
    """\
    Filter genes based on number of cells or counts.

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
        Depending on what was tresholded (`counts` or `cells`), the array stores
        `n_counts` or `n_cells` per gene.
    """
    if copy:
        logg.warning('`copy` is deprecated, use `inplace` instead.')
    n_given_options = sum(
        option is not None for option in [min_cells, min_counts, max_cells, max_counts]
    )
    if n_given_options != 1:
        raise ValueError(
            'Only provide one of the optional parameters `min_counts`, '
            '`min_cells`, `max_counts`, `max_cells` per call.'
        )

    if isinstance(data, AnnData):
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
            adata.var['n_counts'] = number
        else:
            adata.var['n_cells'] = number
        adata._inplace_subset_var(gene_subset)
        return adata if copy else None

    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_cells is None else min_cells
    max_number = max_counts if max_cells is None else max_cells
    number_per_gene = np.sum(
        X if min_cells is None and max_cells is None else X > 0, axis=0
    )
    if issparse(X):
        number_per_gene = number_per_gene.A1
    if min_number is not None:
        gene_subset = number_per_gene >= min_number
    if max_number is not None:
        gene_subset = number_per_gene <= max_number

    s = np.sum(~gene_subset)
    if s > 0:
        msg = f'filtered out {s} genes that are detected '
        if min_cells is not None or min_counts is not None:
            msg += 'in less than '
            msg += (
                f'{min_cells} cells' if min_counts is None else f'{min_counts} counts'
            )
        if max_cells is not None or max_counts is not None:
            msg += 'in more than '
            msg += (
                f'{max_cells} cells' if max_counts is None else f'{max_counts} counts'
            )
        logg.info(msg)
    return gene_subset, number_per_gene


@singledispatch
def log1p(
    X: Union[AnnData, np.ndarray, spmatrix],
    *,
    base: Optional[Number] = None,
    copy: bool = False,
    chunked: bool = None,
    chunk_size: Optional[int] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
):
    """\
    Logarithmize the data matrix.

    Computes :math:`X = \\log(X + 1)`,
    where :math:`log` denotes the natural logarithm unless a different base is given.

    Parameters
    ----------
    X
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
        Entry of layers to tranform.
    obsm
        Entry of obsm to transform.

    Returns
    -------
    Returns or updates `data`, depending on `copy`.
    """
    _check_array_function_arguments(
        chunked=chunked, chunk_size=chunk_size, layer=layer, obsm=obsm
    )
    return log1p_array(X, copy=copy, base=base)


@log1p.register(spmatrix)
def log1p_sparse(X, *, base: Optional[Number] = None, copy: bool = False):
    X = check_array(
        X, accept_sparse=("csr", "csc"), dtype=(np.float64, np.float32), copy=copy
    )
    X.data = log1p(X.data, copy=False, base=base)
    return X


@log1p.register(np.ndarray)
def log1p_array(X, *, base: Optional[Number] = None, copy: bool = False):
    # Can force arrays to be np.ndarrays, but would be useful to not
    # X = check_array(X, dtype=(np.float64, np.float32), ensure_2d=False, copy=copy)
    if copy:
        if not np.issubdtype(X.dtype, np.floating):
            X = X.astype(float)
        else:
            X = X.copy()
    elif not (np.issubdtype(X.dtype, np.floating) or np.issubdtype(X.dtype, complex)):
        X = X.astype(float)
    np.log1p(X, out=X)
    if base is not None:
        np.divide(X, np.log(base), out=X)
    return X


@log1p.register(AnnData)
def log1p_anndata(
    adata,
    *,
    base: Optional[Number] = None,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
) -> Optional[AnnData]:
    if "log1p" in adata.uns_keys():
        logg.warning("adata.X seems to be already log-transformed.")

    adata = adata.copy() if copy else adata
    view_to_actual(adata)

    if chunked:
        if (layer is not None) or (obsm is not None):
            raise NotImplementedError(
                "Currently cannot perform chunked operations on arrays not stored in X."
            )
        for chunk, start, end in adata.chunked_X(chunk_size):
            adata.X[start:end] = log1p(chunk, base=base, copy=False)
    else:
        X = _get_obs_rep(adata, layer=layer, obsm=obsm)
        X = log1p(X, copy=False, base=base)
        _set_obs_rep(adata, X, layer=layer, obsm=obsm)

    adata.uns["log1p"] = {"base": base}
    if copy:
        return adata


def sqrt(
    data: AnnData,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    """\
    Square root the data matrix.

    Computes :math:`X = \\sqrt(X)`.

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
    if not issparse(X):
        return np.sqrt(X)
    else:
        return X.sqrt()


def normalize_per_cell(
    data: Union[AnnData, np.ndarray, spmatrix],
    counts_per_cell_after: Optional[float] = None,
    counts_per_cell: Optional[np.ndarray] = None,
    key_n_counts: str = 'n_counts',
    copy: bool = False,
    layers: Union[Literal['all'], Iterable[str]] = (),
    use_rep: Optional[Literal['after', 'X']] = None,
    min_counts: int = 1,
) -> Optional[AnnData]:
    """\
    Normalize total counts per cell.

    .. warning::
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

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

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
    Returns or updates `adata` with normalized version of the original
    `adata.X`, depending on `copy`.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = AnnData(np.array([[1, 0], [3, 0], [5, 6]]))
    >>> print(adata.X.sum(axis=1))
    [  1.   3.  11.]
    >>> sc.pp.normalize_per_cell(adata)
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts
    0       1.0
    1       3.0
    2      11.0
    [ 3.  3.  3.]
    >>> sc.pp.normalize_per_cell(
    >>>     adata, counts_per_cell_after=1,
    >>>     key_n_counts='n_counts2',
    >>> )
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts  n_counts2
    0       1.0        3.0
    1       3.0        3.0
    2      11.0        3.0
    [ 1.  1.  1.]
    """
    if isinstance(data, AnnData):
        start = logg.info('normalizing by total count per cell')
        adata = data.copy() if copy else data
        if counts_per_cell is None:
            cell_subset, counts_per_cell = materialize_as_ndarray(
                filter_cells(adata.X, min_counts=min_counts)
            )
            adata.obs[key_n_counts] = counts_per_cell
            adata._inplace_subset_obs(cell_subset)
            counts_per_cell = counts_per_cell[cell_subset]
        normalize_per_cell(adata.X, counts_per_cell_after, counts_per_cell)

        layers = adata.layers.keys() if layers == 'all' else layers
        if use_rep == 'after':
            after = counts_per_cell_after
        elif use_rep == 'X':
            after = np.median(counts_per_cell[cell_subset])
        elif use_rep is None:
            after = None
        else:
            raise ValueError('use_rep should be "after", "X" or None')
        for layer in layers:
            subset, counts = filter_cells(adata.layers[layer], min_counts=min_counts)
            temp = normalize_per_cell(adata.layers[layer], after, counts, copy=True)
            adata.layers[layer] = temp

        logg.info(
            '    finished ({time_passed}): normalized adata.X and added'
            f'    {key_n_counts!r}, counts per cell before normalization (adata.obs)',
            time=start,
        )
        return adata if copy else None
    # proceed with data matrix
    X = data.copy() if copy else data
    if counts_per_cell is None:
        if not copy:
            raise ValueError('Can only be run with copy=True')
        cell_subset, counts_per_cell = filter_cells(X, min_counts=min_counts)
        X = X[cell_subset]
        counts_per_cell = counts_per_cell[cell_subset]
    if counts_per_cell_after is None:
        counts_per_cell_after = np.median(counts_per_cell)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        counts_per_cell += counts_per_cell == 0
        counts_per_cell /= counts_per_cell_after
        if not issparse(X):
            X /= materialize_as_ndarray(counts_per_cell[:, np.newaxis])
        else:
            sparsefuncs.inplace_row_scale(X, 1 / counts_per_cell)
    return X if copy else None


def regress_out(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    n_jobs: Optional[int] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Regress out (mostly) unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R [Satija15]. Note that this function tends to overcorrect
    in certain circumstances as described in :issue:`526`.

    Parameters
    ----------
    adata
        The annotated data matrix.
    keys
        Keys for observation annotation on which to regress on.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    copy
        Determines whether a copy of `adata` is returned.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with the corrected data matrix.
    """
    start = logg.info(f'regressing out {keys}')
    if issparse(adata.X):
        logg.info('    sparse input is densified and may ' 'lead to high memory use')
    adata = adata.copy() if copy else adata

    sanitize_anndata(adata)

    # TODO: This should throw an implicit modification warning
    if adata.is_view:
        adata._init_as_actual(adata.copy())

    if isinstance(keys, str):
        keys = [keys]

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    n_jobs = sett.n_jobs if n_jobs is None else n_jobs

    # regress on a single categorical variable
    variable_is_categorical = False
    if keys[0] in adata.obs_keys() and is_categorical_dtype(adata.obs[keys[0]]):
        if len(keys) > 1:
            raise ValueError(
                'If providing categorical variable, '
                'only a single one is allowed. For this one '
                'we regress on the mean for each category.'
            )
        logg.debug('... regressing on per-gene means within categories')
        regressors = np.zeros(adata.X.shape, dtype='float32')
        for category in adata.obs[keys[0]].cat.categories:
            mask = (category == adata.obs[keys[0]]).values
            for ix, x in enumerate(adata.X.T):
                regressors[mask, ix] = x[mask].mean()
        variable_is_categorical = True
    # regress on one or several ordinal variables
    else:
        # create data frame with selected keys (if given)
        if keys:
            regressors = adata.obs[keys]
        else:
            regressors = adata.obs.copy()

        # add column of ones at index 0 (first column)
        regressors.insert(0, 'ones', 1.0)

    len_chunk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_chunks = np.ceil(adata.X.shape[1] / len_chunk).astype(int)

    tasks = []
    # split the adata.X matrix by columns in chunks of size n_chunk
    # (the last chunk could be of smaller size than the others)
    chunk_list = np.array_split(adata.X, n_chunks, axis=1)
    if variable_is_categorical:
        regressors_chunk = np.array_split(regressors, n_chunks, axis=1)
    for idx, data_chunk in enumerate(chunk_list):
        # each task is a tuple of a data_chunk eg. (adata.X[:,0:100]) and
        # the regressors. This data will be passed to each of the jobs.
        if variable_is_categorical:
            regres = regressors_chunk[idx]
        else:
            regres = regressors
        tasks.append(tuple((data_chunk, regres, variable_is_categorical)))

    from joblib import Parallel, delayed

    # TODO: figure out how to test that this doesn't oversubscribe resources
    res = Parallel(n_jobs=n_jobs)(delayed(_regress_out_chunk)(task) for task in tasks)

    # res is a list of vectors (each corresponding to a regressed gene column).
    # The transpose is needed to get the matrix in the shape needed
    adata.X = np.vstack(res).T.astype(adata.X.dtype)
    logg.info('    finished', time=start)
    return adata if copy else None


def _regress_out_chunk(data):
    # data is a tuple containing the selected columns from adata.X
    # and the regressors dataFrame
    data_chunk = data[0]
    regressors = data[1]
    variable_is_categorical = data[2]

    responses_chunk_list = []
    import statsmodels.api as sm
    from statsmodels.tools.sm_exceptions import PerfectSeparationError

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
            result = sm.GLM(
                data_chunk[:, col_index], regres, family=sm.families.Gaussian()
            ).fit()
            new_column = result.resid_response
        except PerfectSeparationError:  # this emulates R's behavior
            logg.warning('Encountered PerfectSeparationError, setting to 0 as in R.')
            new_column = np.zeros(data_chunk.shape[0])

        responses_chunk_list.append(new_column)

    return np.vstack(responses_chunk_list)


@singledispatch
def scale(
    X: Union[AnnData, spmatrix, np.ndarray],
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
):
    """\
    Scale data to unit variance and zero mean.

    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and (for zero_center==True) set to 0
        during this operation. In the future, they might be set to NaNs.

    Parameters
    ----------
    X
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        Whether this function should be performed inplace. If an AnnData object
        is passed, this also determines if a copy is returned.
    layer
        If provided, which element of layers to scale.
    obsm
        If provided, which element of obsm to scale.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with a scaled `adata.X`,
    annotated with `'mean'` and `'std'` in `adata.var`.
    """
    _check_array_function_arguments(layer=layer, obsm=obsm)
    if layer is not None:
        raise ValueError(f"`layer` argument inappropriate for value of type {type(X)}")
    if obsm is not None:
        raise ValueError(f"`obsm` argument inappropriate for value of type {type(X)}")
    return scale_array(X, zero_center=zero_center, max_value=max_value, copy=copy)


@scale.register(np.ndarray)
def scale_array(
    X,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    return_mean_std: bool = False,
):
    if copy:
        X = X.copy()
    if not zero_center and max_value is not None:
        logg.info(  # Be careful of what? This should be more specific
            "... be careful when using `max_value` " "without `zero_center`."
        )

    if np.issubdtype(X.dtype, np.integer):
        logg.info(
            '... as scaling leads to float results, integer '
            'input is cast to float, returning copy.'
        )
        X = X.astype(float)

    mean, var = _get_mean_var(X)
    std = np.sqrt(var)
    std[std == 0] = 1
    if issparse(X):
        if zero_center:
            raise ValueError("Cannot zero-center sparse matrix.")
        sparsefuncs.inplace_column_scale(X, 1 / std)
    else:
        if zero_center:
            X -= mean
        X /= std

    # do the clipping
    if max_value is not None:
        logg.debug(f"... clipping at max_value {max_value}")
        X[X > max_value] = max_value

    if return_mean_std:
        return X, mean, std
    else:
        return X


@scale.register(spmatrix)
def scale_sparse(
    X,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    return_mean_std: bool = False,
):
    # need to add the following here to make inplace logic work
    if zero_center:
        logg.info(
            "... as `zero_center=True`, sparse input is "
            "densified and may lead to large memory consumption"
        )
        X = X.toarray()
        copy = False  # Since the data has been copied
    return scale_array(
        X,
        zero_center=zero_center,
        copy=copy,
        max_value=max_value,
        return_mean_std=return_mean_std,
    )


@scale.register(AnnData)
def scale_anndata(
    adata: AnnData,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata
    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer, obsm=obsm)
    X, adata.var["mean"], adata.var["std"] = scale(
        X,
        zero_center=zero_center,
        max_value=max_value,
        copy=False,  # because a copy has already been made, if it were to be made
        return_mean_std=True,
    )
    _set_obs_rep(adata, X, layer=layer, obsm=obsm)
    if copy:
        return adata


def subsample(
    data: Union[AnnData, np.ndarray, spmatrix],
    fraction: Optional[float] = None,
    n_obs: Optional[int] = None,
    random_state: AnyRandom = 0,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Subsample to a fraction of the number of observations.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    fraction
        Subsample to this `fraction` of the number of observations.
    n_obs
        Subsample to this number of observations.
    random_state
        Random seed to change subsampling.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.

    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~anndata.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """
    np.random.seed(random_state)
    old_n_obs = data.n_obs if isinstance(data, AnnData) else data.shape[0]
    if n_obs is not None:
        new_n_obs = n_obs
    elif fraction is not None:
        if fraction > 1 or fraction < 0:
            raise ValueError(f'`fraction` needs to be within [0, 1], not {fraction}')
        new_n_obs = int(fraction * old_n_obs)
        logg.debug(f'... subsampled to {new_n_obs} data points')
    else:
        raise ValueError('Either pass `n_obs` or `fraction`.')
    obs_indices = np.random.choice(old_n_obs, size=new_n_obs, replace=False)
    if isinstance(data, AnnData):
        if copy:
            return data[obs_indices].copy()
        else:
            data._inplace_subset_obs(obs_indices)
    else:
        X = data
        return X[obs_indices], obs_indices


@deprecated_arg_names({"target_counts": "counts_per_cell"})
def downsample_counts(
    adata: AnnData,
    counts_per_cell: Optional[Union[int, Collection[int]]] = None,
    total_counts: Optional[int] = None,
    *,
    random_state: AnyRandom = 0,
    replace: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Downsample counts from count matrix.

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
    Depending on `copy` returns or updates an `adata` with downsampled `.X`.
    """
    # This logic is all dispatch
    total_counts_call = total_counts is not None
    counts_per_cell_call = counts_per_cell is not None
    if total_counts_call is counts_per_cell_call:
        raise ValueError(
            "Must specify exactly one of `total_counts` or `counts_per_cell`."
        )
    if copy:
        adata = adata.copy()
    if total_counts_call:
        adata.X = _downsample_total_counts(adata.X, total_counts, random_state, replace)
    elif counts_per_cell_call:
        adata.X = _downsample_per_cell(adata.X, counts_per_cell, random_state, replace)
    if copy:
        return adata


def _downsample_per_cell(X, counts_per_cell, random_state, replace):
    n_obs = X.shape[0]
    if isinstance(counts_per_cell, int):
        counts_per_cell = np.full(n_obs, counts_per_cell)
    else:
        counts_per_cell = np.asarray(counts_per_cell)
    # np.random.choice needs int arguments in numba code:
    counts_per_cell = counts_per_cell.astype(np.int_, copy=False)
    if not isinstance(counts_per_cell, np.ndarray) or len(counts_per_cell) != n_obs:
        raise ValueError(
            "If provided, 'counts_per_cell' must be either an integer, or "
            "coercible to an `np.ndarray` of length as number of observations"
            " by `np.asarray(counts_per_cell)`."
        )
    if issparse(X):
        original_type = type(X)
        if not isspmatrix_csr(X):
            X = csr_matrix(X)
        totals = np.ravel(X.sum(axis=1))  # Faster for csr matrix
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
        if original_type is not csr_matrix:  # Put it back
            X = original_type(X)
    else:
        totals = np.ravel(X.sum(axis=1))
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


def _downsample_total_counts(X, total_counts, random_state, replace):
    total_counts = int(total_counts)
    total = X.sum()
    if total < total_counts:
        return X
    if issparse(X):
        original_type = type(X)
        if not isspmatrix_csr(X):
            X = csr_matrix(X)
        _downsample_array(
            X.data,
            total_counts,
            random_state=random_state,
            replace=replace,
            inplace=True,
        )
        X.eliminate_zeros()
        if original_type is not csr_matrix:
            X = original_type(X)
    else:
        v = X.reshape(np.multiply(*X.shape))
        _downsample_array(v, total_counts, random_state, replace=replace, inplace=True)
    return X


@numba.njit(cache=True)
def _downsample_array(
    col: np.ndarray,
    target: int,
    random_state: AnyRandom = 0,
    replace: bool = True,
    inplace: bool = False,
):
    """\
    Evenly reduce counts in cell to target amount.

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


# --------------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------------


def _pca_fallback(data, n_comps=2):
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    C = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since C is symmetric,
    # the performance gain is substantial
    # evals, evecs = np.linalg.eigh(C)
    evals, evecs = sp.sparse.linalg.eigsh(C, k=n_comps)
    # sort eigenvalues in decreasing order
    idcs = np.argsort(evals)[::-1]
    evecs = evecs[:, idcs]
    evals = evals[idcs]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or n_comps)
    evecs = evecs[:, :n_comps]
    # project data points on eigenvectors
    return np.dot(evecs.T, data.T).T

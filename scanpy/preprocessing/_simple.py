"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""
import warnings
from typing import Union, Optional, Tuple, Collection

import numba
import numpy as np
import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from sklearn.utils import sparsefuncs
from pandas.api.types import is_categorical_dtype
from anndata import AnnData

from .._settings import settings as sett
from .. import logging as logg
from ..utils import sanitize_anndata, deprecated_arg_names
from ._distributed import materialize_as_ndarray
from ._utils import _get_mean_var

# install dask if available
try:
    import dask.array as da
except ImportError:
    da = None

N_PCS = 50  # default number of PCs

# backwards compat
from ._deprecated.highly_variable_genes import filter_genes_dispersion


def filter_cells(
    data: AnnData,
    min_counts: Optional[int] = None,
    min_genes:  Optional[int] = None,
    max_counts: Optional[int] = None,
    max_genes:  Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters ``min_counts``, ``min_genes``,
    ``max_counts``, ``max_genes`` per call.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
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
    Depending on ``inplace``, returns the following arrays or directly subsets
    and annotates the data matrix:

    cells_subset : numpy.ndarray
        Boolean index mask that does filtering. ``True`` means that the
        cell is kept. ``False`` means the cell is removed.
    number_per_cell : numpy.ndarray
        Depending on what was tresholded (``counts`` or ``genes``), the array stores
        ``n_counts`` or ``n_cells`` per gene.

    Examples
    --------
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
       logg.warn('`copy` is deprecated, use `inplace` instead.')
    n_given_options = sum(
        option is not None for option in
        [min_genes, min_counts, max_genes, max_counts])
    if n_given_options != 1:
        raise ValueError(
            'Only provide one of the optional parameters `min_counts`,'
            '`min_genes`, `max_counts`, `max_genes` per call.')
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        cell_subset, number = materialize_as_ndarray(filter_cells(adata.X, min_counts, min_genes, max_counts, max_genes))
        if not inplace:
            return cell_subset, number
        if min_genes is None and max_genes is None: adata.obs['n_counts'] = number
        else: adata.obs['n_genes'] = number
        adata._inplace_subset_obs(cell_subset)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    max_number = max_counts if max_genes is None else max_genes
    number_per_cell = np.sum(X if min_genes is None and max_genes is None
                             else X > 0, axis=1)
    if issparse(X): number_per_cell = number_per_cell.A1
    if min_number is not None:
        cell_subset = number_per_cell >= min_number
    if max_number is not None:
        cell_subset = number_per_cell <= max_number

    s = np.sum(~cell_subset)
    if s > 0:
        logg.info('filtered out {} cells that have'.format(s), end=' ')
        if min_genes is not None or min_counts is not None:
            logg.info('less than',
                   str(min_genes) + ' genes expressed'
                   if min_counts is None else str(min_counts) + ' counts', no_indent=True)
        if max_genes is not None or max_counts is not None:
            logg.info('more than ',
                   str(max_genes) + ' genes expressed'
                   if max_counts is None else str(max_counts) + ' counts', no_indent=True)
    return cell_subset, number_per_cell


def filter_genes(
    data: AnnData,
    min_counts: Optional[int] = None,
    min_cells:  Optional[int] = None,
    max_counts: Optional[int] = None,
    max_cells:  Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Union[AnnData, None, Tuple[np.ndarray, np.ndarray]]:
    """Filter genes based on number of cells or counts.

    Keep genes that have at least ``min_counts`` counts or are expressed in at
    least ``min_cells`` cells or have at most ``max_counts`` counts or are expressed
    in at most ``max_cells`` cells.

    Only provide one of the optional parameters ``min_counts``, ``min_cells``,
    ``max_counts``, ``max_cells`` per call.

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

    gene_subset : numpy.ndarray
        Boolean index mask that does filtering. `True` means that the
        gene is kept. `False` means the gene is removed.
    number_per_gene : numpy.ndarray
        Depending on what was tresholded (`counts` or `cells`), the array stores
        `n_counts` or `n_cells` per gene.
    """
    if copy:
       logg.warn('`copy` is deprecated, use `inplace` instead.')
    n_given_options = sum(
        option is not None for option in
        [min_cells, min_counts, max_cells, max_counts])
    if n_given_options != 1:
        raise ValueError(
            'Only provide one of the optional parameters `min_counts`,'
            '`min_cells`, `max_counts`, `max_cells` per call.')

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        gene_subset, number = materialize_as_ndarray(
            filter_genes(adata.X, min_cells=min_cells,
                         min_counts=min_counts, max_cells=max_cells,
                         max_counts=max_counts))
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
    number_per_gene = np.sum(X if min_cells is None and max_cells is None
                             else X > 0, axis=0)
    if issparse(X):
        number_per_gene = number_per_gene.A1
    if min_number is not None:
        gene_subset = number_per_gene >= min_number
    if max_number is not None:
        gene_subset = number_per_gene <= max_number

    s = np.sum(~gene_subset)
    if s > 0:
        logg.info('filtered out {} genes that are detected'.format(s), end=' ')
        if min_cells is not None or min_counts is not None:
            logg.info('in less than',
                   str(min_cells) + ' cells'
                   if min_counts is None else str(min_counts) + ' counts', no_indent=True)
        if max_cells is not None or max_counts is not None:
            logg.info('in more than ',
                   str(max_cells) + ' cells'
                   if max_counts is None else str(max_counts) + ' counts', no_indent=True)
    return gene_subset, number_per_gene


def log1p(
    data: Union[AnnData, np.ndarray, spmatrix],
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    """Logarithmize the data matrix.

    Computes :math:`X = \\log(X + 1)`, where :math:`log` denotes the natural logarithm.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    chunked
        Process the data matrix in chunks, which will save memory.
        Applies only to :class:`~anndata.AnnData`.
    chunk_size
        ``n_obs`` of the chunks to process the data in.

    Returns
    -------
    Returns or updates ``data``, depending on ``copy``.
    """
    if copy:
        if not isinstance(data, AnnData):
            data = data.astype(np.floating)
        else:
            data = data.copy()
    elif not isinstance(data, AnnData) and np.issubdtype(data.dtype, np.integer):
        raise TypeError("Cannot perform inplace log1p on integer array")

    def _log1p(X):
        if issparse(X):
            np.log1p(X.data, out=X.data)
        else:
            np.log1p(X, out=X)
        return X

    if isinstance(data, AnnData):
        if not np.issubdtype(data.X.dtype, np.floating):
            data.X = data.X.astype(np.float32)
        if chunked:
            for chunk, start, end in data.chunked_X(chunk_size):
                 data.X[start:end] = _log1p(chunk)
        else:
            _log1p(data.X)
    else:
        _log1p(data)

    return data if copy else None


def sqrt(
    data: AnnData,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    """Square root the data matrix.

    Computes :math:`X = \\sqrt(X)`.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    copy
        If an :class:`~scanpy.api.AnnData` is passed,
        determines whether a copy is returned.
    chunked
        Process the data matrix in chunks, which will save memory.
        Applies only to :class:`~anndata.AnnData`.
    chunk_size
        ``n_obs`` of the chunks to process the data in.

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


def pca(
    data: Union[AnnData, np.ndarray, spmatrix],
    n_comps: int = N_PCS,
    zero_center: Optional[bool] = True,
    svd_solver: str = 'auto',
    random_state: int = 0,
    return_info: bool = False,
    use_highly_variable: Optional[bool] = None,
    dtype: str = 'float32',
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Union[AnnData, np.ndarray, spmatrix]:
    """Principal component analysis [Pedregosa11]_.

    Computes PCA coordinates, loadings and variance decomposition. Uses the
    implementation of *scikit-learn* [Pedregosa11]_.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute.
    zero_center
        If `True`, compute standard PCA from covariance matrix.
        If ``False``, omit zero-centering variables
        (uses :class:`~sklearn.decomposition.TruncatedSVD`),
        which allows to handle sparse input efficiently.
        Passing ``None`` decides automatically based on sparseness of the data.
    svd_solver
        SVD solver to use:

        ``'arpack'``
          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)

        ``'randomized'``
          for the randomized algorithm due to Halko (2009).

        ``'auto'`` (the default)
          chooses automatically depending on the size of the problem.

    random_state
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        ``.var['highly_variable']``.
        By default uses them if they have been determined beforehand.
    dtype
        Numpy data type string to which to convert the result.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    chunked
        If ``True``, perform an incremental PCA on segments of ``chunk_size``.
        The incremental PCA automatically zero centers and ignores settings of
        ``random_seed`` and ``svd_solver``. If ``False``, perform a full PCA.
    chunk_size
        Number of observations to include in each chunk.
        Required if ``chunked=True`` was passed.

    Returns
    -------
    X_pca : :class:`scipy.sparse.spmatrix` or :class:`numpy.ndarray`
        If `data` is array-like and ``return_info=False`` was passed,
        this function only returns `X_pca`…
    adata : anndata.AnnData
        …otherwise if ``copy=True`` it returns or else adds fields to ``adata``:

        ``.obsm['X_pca']``
             PCA representation of data.

        ``.varm['PCs']``
             The principal components containing the loadings.

        ``.uns['pca']['variance_ratio']``)
             Ratio of explained variance.

        ``.uns['pca']['variance']``
             Explained variance, equivalent to the eigenvalues of the covariance matrix.
    """
    # chunked calculation is not randomized, anyways
    if svd_solver in {'auto', 'randomized'} and not chunked:
        logg.info(
            'Note that scikit-learn\'s randomized PCA might not be exactly '
            'reproducible across different computational platforms. For exact '
            'reproducibility, choose `svd_solver=\'arpack\'.` This will likely '
            'become the Scanpy default in the future.')

    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data)

    logg.msg('computing PCA with n_comps =', n_comps, r=True, v=4)

    if adata.n_vars < n_comps:
        n_comps = adata.n_vars - 1
        logg.msg('reducing number of computed PCs to',
               n_comps, 'as dim of data is only', adata.n_vars, v=4)

    if use_highly_variable is True and 'highly_variable' not in adata.var.keys():
        raise ValueError('Did not find adata.var[\'highly_variable\']. '
                         'Either your data already only consists of highly-variable genes '
                         'or consider running `pp.filter_genes_dispersion` first.')
    if use_highly_variable is None:
        use_highly_variable = True if 'highly_variable' in adata.var.keys() else False
    if use_highly_variable:
        logg.info('computing PCA on highly variable genes')
    adata_comp = adata[:, adata.var['highly_variable']] if use_highly_variable else adata

    if chunked:
        if not zero_center or random_state or svd_solver != 'auto':
            logg.msg('Ignoring zero_center, random_state, svd_solver', v=4)

        from sklearn.decomposition import IncrementalPCA

        X_pca = np.zeros((adata_comp.X.shape[0], n_comps), adata_comp.X.dtype)

        pca_ = IncrementalPCA(n_components=n_comps)

        for chunk, _, _ in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            pca_.partial_fit(chunk)

        for chunk, start, end in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            X_pca[start:end] = pca_.transform(chunk)
    else:
        if zero_center is None:
            zero_center = not issparse(adata_comp.X)
        if zero_center:
            from sklearn.decomposition import PCA
            if issparse(adata_comp.X):
                logg.msg('    as `zero_center=True`, '
                       'sparse input is densified and may '
                       'lead to huge memory consumption', v=4)
                X = adata_comp.X.toarray()  # Copying the whole adata_comp.X here, could cause memory problems
            else:
                X = adata_comp.X
            pca_ = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
        else:
            from sklearn.decomposition import TruncatedSVD
            logg.msg('    without zero-centering: \n'
                   '    the explained variance does not correspond to the exact statistical defintion\n'
                   '    the first component, e.g., might be heavily influenced by different means\n'
                   '    the following components often resemble the exact PCA very closely', v=4)
            pca_ = TruncatedSVD(n_components=n_comps, random_state=random_state)
            X = adata_comp.X
        X_pca = pca_.fit_transform(X)

    if X_pca.dtype.descr != np.dtype(dtype).descr: X_pca = X_pca.astype(dtype)

    if data_is_AnnData:
        adata.obsm['X_pca'] = X_pca
        if use_highly_variable:
            adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm['PCs'][adata.var['highly_variable']] = pca_.components_.T
        else:
            adata.varm['PCs'] = pca_.components_.T
        adata.uns['pca'] = {}
        adata.uns['pca']['variance'] = pca_.explained_variance_
        adata.uns['pca']['variance_ratio'] = pca_.explained_variance_ratio_
        logg.msg('    finished', t=True, end=' ', v=4)
        logg.msg('and added\n'
                 '    \'X_pca\', the PCA coordinates (adata.obs)\n'
                 '    \'PC1\', \'PC2\', ..., the loadings (adata.var)\n'
                 '    \'pca_variance\', the variance / eigenvalues (adata.uns)\n'
                 '    \'pca_variance_ratio\', the variance ratio (adata.uns)', v=4)
        return adata if copy else None
    else:
        if return_info:
            return X_pca, pca_.components_, pca_.explained_variance_ratio_, pca_.explained_variance_
        else:
            return X_pca


def normalize_per_cell(
    data,
    counts_per_cell_after=None,
    counts_per_cell=None,
    key_n_counts=None,
    copy=False,
    layers=[],
    use_rep=None,
    min_counts=1,
) -> Optional[AnnData]:
    """Normalize total counts per cell.

    .. warning::
        .. deprecated:: 1.3.7
            Use :func:`~scanpy.api.pp.normalize_total` instead.
            The new function is equivalent to the present
            function, except that

            * the new function doesn't filter cells based on `min_counts`,
              use :func:`~scanpy.api.pp.filter_cells` if filtering is needed.
            * some arguments were renamed
            * `copy` is replaced by `inplace`

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    key_n_counts : `str`, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    min_counts : `int`, optional (default: 1)
        Cells with counts less than `min_counts` are filtered out during
        normalization.

    Returns
    -------
    Returns or updates `adata` with normalized version of the original
    `adata.X`, depending on `copy`.

    Examples
    --------
    >>> adata = AnnData(
    >>>     data=np.array([[1, 0], [3, 0], [5, 6]]))
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
    >>> sc.pp.normalize_per_cell(adata, counts_per_cell_after=1,
    >>>                          key_n_counts='n_counts2')
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts  n_counts2
    0       1.0        3.0
    1       3.0        3.0
    2      11.0        3.0
    [ 1.  1.  1.]
    """
    if key_n_counts is None: key_n_counts = 'n_counts'
    if isinstance(data, AnnData):
        logg.msg('normalizing by total count per cell', r=True)
        adata = data.copy() if copy else data
        if counts_per_cell is None:
            cell_subset, counts_per_cell = materialize_as_ndarray(
                        filter_cells(adata.X, min_counts=min_counts))
            adata.obs[key_n_counts] = counts_per_cell
            adata._inplace_subset_obs(cell_subset)
            counts_per_cell=counts_per_cell[cell_subset]
        normalize_per_cell(adata.X, counts_per_cell_after, counts_per_cell)

        layers = adata.layers.keys() if layers == 'all' else layers
        if use_rep == 'after':
            after = counts_per_cell_after
        elif use_rep == 'X':
            after = np.median(counts_per_cell[cell_subset])
        elif use_rep is None:
            after = None
        else: raise ValueError('use_rep should be "after", "X" or None')
        for layer in layers:
            subset, counts = filter_cells(adata.layers[layer],
                    min_counts=min_counts)
            temp = normalize_per_cell(adata.layers[layer], after, counts, copy=True)
            adata.layers[layer] = temp

        logg.msg('    finished', t=True, end=': ')
        logg.msg('normalized adata.X and added', no_indent=True)
        logg.msg('    \'{}\', counts per cell before normalization (adata.obs)'
            .format(key_n_counts))
        return adata if copy else None
    # proceed with data matrix
    X = data.copy() if copy else data
    if counts_per_cell is None:
        if copy == False:
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
        if not issparse(X): X /= materialize_as_ndarray(counts_per_cell[:, np.newaxis])
        else: sparsefuncs.inplace_row_scale(X, 1/counts_per_cell)
    return X if copy else None


def normalize_per_cell_weinreb16_deprecated(
    X,
    max_fraction=1,
    mult_with_mean=False,
) -> np.ndarray:
    """Normalize each cell [Weinreb17]_.

    This is a deprecated version. See `normalize_per_cell` instead.

    Normalize each cell by UMI count, so that every cell has the same total
    count.

    Parameters
    ----------
    X : np.ndarray
        Expression matrix. Rows correspond to cells and columns to genes.
    max_fraction : float, optional
        Only use genes that make up more than max_fraction of the total
        reads in every cell.
    mult_with_mean: bool, optional
        Multiply the result with the mean of total counts.

    Returns
    -------
    Normalized version of the original expression matrix.
    """
    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('Choose max_fraction between 0 and 1.')

    counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
    gene_subset = np.all(X <= counts_per_cell[:, None] * max_fraction, axis=0)
    if issparse(X): gene_subset = gene_subset.A1
    tc_include = X[:, gene_subset].sum(1).A1 if issparse(X) else X[:, gene_subset].sum(1)

    X_norm = X.multiply(csr_matrix(1/tc_include[:, None])) if issparse(X) else X / tc_include[:, None]
    if mult_with_mean:
        X_norm *= np.mean(counts_per_cell)

    return X_norm


def regress_out(adata, keys, n_jobs=None, copy=False) -> Optional[AnnData]:
    """Regress out unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R [Satija15].

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    keys : `str` or list of `str`
        Keys for observation annotation on which to regress on.
    n_jobs : `int` or `None`, optional. If None is given, then the n_jobs seting is used (default: `None`)
        Number of jobs for parallel computation.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with the corrected data matrix.
    """
    logg.info('regressing out', keys, r=True)
    if issparse(adata.X):
        logg.info('    sparse input is densified and may '
                  'lead to high memory use')
    adata = adata.copy() if copy else adata
    if isinstance(keys, str):
        keys = [keys]

    if issparse(adata.X):
        adata.X = adata.X.toarray()

    n_jobs = sett.n_jobs if n_jobs is None else n_jobs

    # regress on a single categorical variable
    sanitize_anndata(adata)
    variable_is_categorical = False
    if keys[0] in adata.obs_keys() and is_categorical_dtype(adata.obs[keys[0]]):
        if len(keys) > 1:
            raise ValueError(
                'If providing categorical variable, '
                'only a single one is allowed. For this one '
                'we regress on the mean for each category.')
        logg.msg('... regressing on per-gene means within categories')
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
    # split the adata.X matrix by columns in chunks of size n_chunk (the last chunk could be of smaller
    # size than the others)
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

    if n_jobs > 1 and n_chunks > 1:
        import multiprocessing
        pool = multiprocessing.Pool(n_jobs)
        res = pool.map_async(_regress_out_chunk, tasks).get(9999999)
        pool.close()

    else:
        res = list(map(_regress_out_chunk, tasks))

    # res is a list of vectors (each corresponding to a regressed gene column).
    # The transpose is needed to get the matrix in the shape needed
    adata.X = np.vstack(res).T.astype(adata.X.dtype)
    logg.info('    finished', t=True)
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
        if variable_is_categorical:
            regres = np.c_[np.ones(regressors.shape[0]), regressors[:, col_index]]
        else:
            regres = regressors
        try:
            result = sm.GLM(data_chunk[:, col_index], regres, family=sm.families.Gaussian()).fit()
            new_column = result.resid_response
        except PerfectSeparationError:  # this emulates R's behavior
            logg.warn('Encountered PerfectSeparationError, setting to 0 as in R.')
            new_column = np.zeros(data_chunk.shape[0])

        responses_chunk_list.append(new_column)

    return np.vstack(responses_chunk_list)


def scale(data, zero_center=True, max_value=None, copy=False) -> Optional[AnnData]:
    """Scale data to unit variance and zero mean.

    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and set to 0 during this operation. In
        the future, they might be set to NaNs.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    zero_center : `bool`, optional (default: `True`)
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value : `float` or `None`, optional (default: `None`)
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with a scaled `adata.X`.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        # need to add the following here to make inplace logic work
        if zero_center and issparse(adata.X):
            logg.msg(
                '... scale_data: as `zero_center=True`, sparse input is '
                'densified and may lead to large memory consumption')
            adata.X = adata.X.toarray()
        scale(adata.X, zero_center=zero_center, max_value=max_value, copy=False)
        return adata if copy else None
    X = data.copy() if copy else data  # proceed with the data matrix
    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
    if not zero_center and max_value is not None:
        logg.msg(
            '... scale_data: be careful when using `max_value` without `zero_center`',
            v=4)
    if max_value is not None:
        logg.msg('... clipping at max_value', max_value)
    if zero_center and issparse(X):
        logg.msg('... scale_data: as `zero_center=True`, sparse input is '
                 'densified and may lead to large memory consumption, returning copy',
                 v=4)
        X = X.toarray()
        copy = True
    _scale(X, zero_center)
    if max_value is not None: X[X > max_value] = max_value
    return X if copy else None


def subsample(data, fraction=None, n_obs=None, random_state=0, copy=False) -> Optional[AnnData]:
    """Subsample to a fraction of the number of observations.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    fraction : `float` in [0, 1] or `None`, optional (default: `None`)
        Subsample to this `fraction` of the number of observations.
    n_obs : `int` or `None`, optional (default: `None`)
        Subsample to this number of observations.
    random_state : `int` or `None`, optional (default: 0)
        Random seed to change subsampling.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

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
            raise ValueError('`fraction` needs to be within [0, 1], not {}'
                             .format(fraction))
        new_n_obs = int(fraction * old_n_obs)
        logg.msg('... subsampled to {} data points'.format(new_n_obs))
    else:
        raise ValueError('Either pass `n_obs` or `fraction`.')
    obs_indices = np.random.choice(old_n_obs, size=new_n_obs, replace=False)
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata._inplace_subset_obs(obs_indices)
        return adata if copy else None
    else:
        X = data
        return X[obs_indices], obs_indices


@deprecated_arg_names({"target_counts": "counts_per_cell"})
def downsample_counts(
    adata: AnnData,
    counts_per_cell: Optional[Union[int, Collection[int]]] = None,
    total_counts: Optional[int] = None,
    random_state: Optional[int] = 0,
    replace: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """Downsample counts from count matrix.

    If `counts_per_cell` is specified, each cell will downsampled. If
    `total_counts` is specified, expression matrix will be downsampled to
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
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates an `adata` with downsampled `.X`.
    """
    # This logic is all dispatch
    total_counts_call = total_counts is not None
    counts_per_cell_call = counts_per_cell is not None
    if total_counts_call is counts_per_cell_call:
        raise ValueError("Must specify exactly one of `total_counts` or `counts_per_cell`.")
    if copy:
        adata = adata.copy()
    adata.X = adata.X.astype(np.integer)  # Numba doesn't want floats
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
        rows = np.split(X.data.view(), X.indptr[1:-1])
        for rowidx in under_target:
            row = rows[rowidx]
            _downsample_array(row, counts_per_cell[rowidx], random_state=random_state,
                              replace=replace, inplace=True)
        X.eliminate_zeros()
        if original_type is not csr_matrix:  # Put it back
            X = original_type(X)
    else:
        totals = np.ravel(X.sum(axis=1))
        under_target = np.nonzero(totals > counts_per_cell)[0]
        for rowidx in under_target:
            row = X[rowidx, :].view()
            _downsample_array(row, counts_per_cell[rowidx], random_state=random_state,
                              replace=replace, inplace=True)
    return X


def _downsample_total_counts(X, total_counts, random_state, replace):
    total = X.sum()
    if total < total_counts:
        return X
    if issparse(X):
        original_type = type(X)
        if not isspmatrix_csr(X):
            X = csr_matrix(X)
        _downsample_array(X.data, total_counts, random_state=random_state,
                          replace=replace, inplace=True)
        X.eliminate_zeros()
        if original_type is not csr_matrix:
            X = original_type(X)
    else:
        v = X.view().reshape(np.multiply(*X.shape))
        _downsample_array(v, total_counts, random_state, replace=replace,
                          inplace=True)
    return X


@numba.njit(cache=True)
def _downsample_array(col: np.array, target: int, random_state: int=0,
                      replace: bool = True, inplace: bool=False):
    """
    Evenly reduce counts in cell to target amount.

    This is an internal function and has some restrictions:

    * `dtype` of col must be an integer (i.e. satisfy issubclass(col.dtype.type, np.integer))
    * total counts in cell must be less than target
    """
    np.random.seed(random_state)
    cumcounts = col.cumsum()
    if inplace:
        col[:] = 0
    else:
        col = np.zeros_like(col)
    total = cumcounts[-1]
    sample = np.random.choice(total, target, replace=replace)
    sample.sort()
    geneptr = 0
    for count in sample:
        while count >= cumcounts[geneptr]:
            geneptr += 1
        col[geneptr] += 1
    return col


def zscore_deprecated(X: np.ndarray) -> np.ndarray:
    """Z-score standardize each variable/gene in X.

    Use `scale` instead.

    Reference: Weinreb et al. (2017).

    Parameters
    ----------
    X
        Data matrix. Rows correspond to cells and columns to genes.

    Returns
    -------
    Z-score standardized version of the data matrix.
    """
    means = np.tile(np.mean(X, axis=0)[None, :], (X.shape[0], 1))
    stds = np.tile(np.std(X, axis=0)[None, :], (X.shape[0], 1))
    return (X - means) / (stds + .0001)


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


def _scale(X, zero_center=True):
    # - using sklearn.StandardScaler throws an error related to
    #   int to long trafo for very large matrices
    # - using X.multiply is slower
    #   the result differs very slightly, why?
    if True:
        mean, var = _get_mean_var(X)
        scale = np.sqrt(var)
        if issparse(X):
            if zero_center: raise ValueError('Cannot zero-center sparse matrix.')
            sparsefuncs.inplace_column_scale(X, 1/scale)
        else:
            X -= mean
            scale[scale == 0] = 1e-12
            X /= scale
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=zero_center, copy=False).partial_fit(X)
        # user R convention (unbiased estimator)
        scaler.scale_ *= np.sqrt(X.shape[0]/(X.shape[0]-1))
        scaler.transform(X)

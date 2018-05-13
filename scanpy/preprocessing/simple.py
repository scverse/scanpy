"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""

import numpy as np
import scipy as sp
import warnings
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from pandas.api.types import is_categorical_dtype
from anndata import AnnData
from .. import settings as sett
from .. import logging as logg
from ..utils import sanitize_anndata

N_PCS = 50  # default number of PCs


def filter_cells(data, min_counts=None, min_genes=None, max_counts=None,
                 max_genes=None, copy=False):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers, i.e.,
    "unreliable" observations.

    Only provide one of the optional parameters `min_counts`, `min_genes`,
    `max_counts`, `max_genes` per call.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts : `int`, optional (default: `None`)
        Minimum number of counts required for a cell to pass filtering.
    min_genes : `int`, optional (default: `None`)
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts : `int`, optional (default: `None`)
        Maximum number of counts required for a cell to pass filtering.
    max_genes : `int`, optional (default: `None`)
        Maximum number of genes expressed required for a cell to pass filtering.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    If `data` is an :class:`~scanpy.api.AnnData`, filters the object and adds\
    either `n_genes` or `n_counts` to `adata.obs`. Otherwise a tuple

    cell_subset : `np.ndarray`
        Boolean index mask that does filtering. `True` means that the cell is
        kept. `False` means the cell is removed.
    number_per_cell : `np.ndarray`
        Either `n_counts` or `n_genes` per cell.

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
    if min_genes is not None and min_counts is not None:
        raise ValueError('Either provide min_counts or min_genes, but not both.')
    if min_genes is not None and max_genes is not None:
        raise ValueError('Either provide min_genes or max_genes, but not both.')
    if min_counts is not None and max_counts is not None:
        raise ValueError('Either provide min_counts or max_counts, but not both.')
    if min_genes is None and min_counts is None and max_genes is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        cell_subset, number = filter_cells(adata.X, min_counts, min_genes, max_counts, max_genes)
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
    logg.msg('filtered out {} cells that have'.format(s), end=' ', v=4)
    if min_genes is not None or min_counts is not None:
        logg.msg('less than',
               str(min_genes) + ' genes expressed'
               if min_counts is None else str(min_counts) + ' counts', v=4, no_indent=True)
    if max_genes is not None or max_counts is not None:
        logg.msg('more than ',
               str(max_genes) + ' genes expressed'
               if max_counts is None else str(max_counts) + ' counts', v=4, no_indent=True)
    return cell_subset, number_per_cell


def filter_genes(data, min_counts=None, min_cells=None, max_counts=None,
                 max_cells=None, copy=False):
    """Filter genes based on number of cells or counts.

    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells or have at most `max_counts` counts or are expressed
    in at most `max_cells` cells.

    Only provide one of the optional parameters `min_counts`, `min_cells`,
    `max_counts`, `max_cells` per call.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts : `int`, optional (default: `None`)
        Minimum number of counts required for a cell to pass filtering.
    min_cells : `int`, optional (default: `None`)
        Minimum number of cells expressed required for a cell to pass filtering.
    max_counts : `int`, optional (default: `None`)
        Maximum number of counts required for a cell to pass filtering.
    max_cells : `int`, optional (default: `None`)
        Maximum number of cells expressed required for a cell to pass filtering.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    If `data` is an :class:`~scanpy.api.AnnData`, filters the object and adds\
    either `n_cells` or `n_counts` to `adata.var`. Otherwise a tuple

    gene_subset : `np.ndarray`
        Boolean index mask that does filtering. `True` means that the gene is
        kept. `False` means the gene is removed.
    number_per_cell : `np.ndarray`
        Either `n_counts` or `n_cells` per cell.
    """
    n_given_options = sum(
        option is not None for option in
        [min_cells, min_counts, max_cells, max_counts])
    if n_given_options != 1:
        raise ValueError(
            'Only provide one of the optional parameters `min_counts`,'
            '`min_cells`, `max_counts`, `max_cells` per call.')

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        gene_subset, number = filter_genes(adata.X, min_cells=min_cells,
                                           min_counts=min_counts, max_cells=max_cells,
                                           max_counts=max_counts)
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
    logg.msg('filtered out {} genes that are detected'.format(s), end=' ', v=4)
    if min_cells is not None or min_counts is not None:
        logg.msg('in less than',
               str(min_cells) + ' cells'
               if min_counts is None else str(min_counts) + ' counts', v=4, no_indent=True)
    if max_cells is not None or max_counts is not None:
        logg.msg('in more than ',
               str(max_cells) + ' cells'
               if max_counts is None else str(max_counts) + ' counts', v=4, no_indent=True)
    return gene_subset, number_per_gene


def filter_genes_dispersion(data,
                            flavor='seurat',
                            min_disp=None, max_disp=None,
                            min_mean=None, max_mean=None,
                            n_bins=20,
                            n_top_genes=None,
                            log=True,
                            copy=False):
    """Extract highly variable genes [Satija15]_ [Zheng17]_.

    If trying out parameters, pass the data matrix instead of AnnData.

    Depending on `flavor`, this reproduces the R-implementations of Seurat
    [Satija15]_ and Cell Ranger [Zheng17]_.

    The normalized dispersion is obtained by scaling with the mean and standard
    deviation of the dispersions for genes falling into a given bin for mean
    expression of genes. This means that for each bin of mean expression, highly
    variable genes are selected.

    Use `flavor='cell_ranger'` with care and in the same way as in
    :func:`~scanpy.api.pp.recipe_zheng17`.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    flavor : {'seurat', 'cell_ranger'}, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion. If choosing
        'seurat', this expects non-logarithmized data - the logarithm of mean
        and dispersion is taken internally when `log` is at its default value
        `True`. For 'cell_ranger', this is usually called for logarithmized data
        - in this case you should set `log` to `False`. In their default
        workflows, Seurat passes the cutoffs whereas Cell Ranger passes
        `n_top_genes`.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=`None` : `float`, optional
        If `n_top_genes` unequals `None`, these cutoffs for the means and the
        normalized dispersions are ignored.
    n_bins : `int` (default: 20)
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    n_top_genes : `int` or `None` (default: `None`)
        Number of highly-variable genes to keep.
    log : `bool`, optional (default: `True`)
        Use the logarithm of the mean to variance ratio.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    If an AnnData `adata` is passed, returns or updates `adata` depending on \
    `copy`. It filters the `adata` and adds the annotations

    means : adata.var
        Means per gene. Logarithmized when `log` is `True`.
    dispersions : adata.var
        Dispersions per gene. Logarithmized when `log` is `True`.
    dispersions_norm : adata.var
        Normalized dispersions per gene. Logarithmized when `log` is `True`.

    If a data matrix `X` is passed, the annotation is returned as `np.recarray` \
    with the same information stored in fields: `gene_subset`, `means`, `dispersions`, `dispersion_norm`.
    """
    if n_top_genes is not None and not all([
            min_disp is None, max_disp is None, min_mean is None, max_mean is None]):
        logg.warn('If you pass `n_top_genes`, all cutoffs are ignored.')
    if min_disp is None: min_disp = 0.5
    if min_mean is None: min_mean = 0.0125
    if max_mean is None: max_mean = 3
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        result = filter_genes_dispersion(adata.X, log=log,
                                         min_disp=min_disp, max_disp=max_disp,
                                         min_mean=min_mean, max_mean=max_mean,
                                         n_top_genes=n_top_genes,
                                         flavor=flavor)
        adata.var['means'] = result['means']
        adata.var['dispersions'] = result['dispersions']
        adata.var['dispersions_norm'] = result['dispersions_norm']
        adata._inplace_subset_var(result['gene_subset'])
        return adata if copy else None
    logg.msg('extracting highly variable genes',
              r=True, v=4)
    X = data  # no copy necessary, X remains unchanged in the following
    mean, var = _get_mean_var(X)
    # now actually compute the dispersion
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    if log:  # logarithmized mean as in Seurat
        dispersion[dispersion == 0] = np.nan
        dispersion = np.log(dispersion)
        mean = np.log1p(mean)
    # all of the following quantities are "per-gene" here
    import pandas as pd
    df = pd.DataFrame()
    df['mean'] = mean
    df['dispersion'] = dispersion
    if flavor == 'seurat':
        df['mean_bin'] = pd.cut(df['mean'], bins=n_bins)
        disp_grouped = df.groupby('mean_bin')['dispersion']
        disp_mean_bin = disp_grouped.mean()
        disp_std_bin = disp_grouped.std(ddof=1)
        # retrieve those genes that have nan std, these are the ones where
        # only a single gene fell in the bin and implicitly set them to have
        # a normalized disperion of 1
        one_gene_per_bin = disp_std_bin.isnull()
        gen_indices = np.where(one_gene_per_bin[df['mean_bin']])[0].tolist()
        if len(gen_indices) > 0:
            logg.msg(
                'Gene indices {} fell into a single bin: their '
                'normalized dispersion was set to 1.\n    '
                'Decreasing `n_bins` will likely avoid this effect.'
                .format(gen_indices), v=4)
        disp_std_bin[one_gene_per_bin] = disp_mean_bin[one_gene_per_bin]
        disp_mean_bin[one_gene_per_bin] = 0
        # actually do the normalization
        df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
                                 - disp_mean_bin[df['mean_bin']].values) \
                                 / disp_std_bin[df['mean_bin']].values
    elif flavor == 'cell_ranger':
        from statsmodels import robust
        df['mean_bin'] = pd.cut(df['mean'], np.r_[-np.inf,
            np.percentile(df['mean'], np.arange(10, 105, 5)), np.inf])
        disp_grouped = df.groupby('mean_bin')['dispersion']
        disp_median_bin = disp_grouped.median()
        # the next line raises the warning: "Mean of empty slice"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            disp_mad_bin = disp_grouped.apply(robust.mad)
        df['dispersion_norm'] = np.abs((df['dispersion'].values
                                 - disp_median_bin[df['mean_bin']].values)) \
                                / disp_mad_bin[df['mean_bin']].values
    else:
        raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
    dispersion_norm = df['dispersion_norm'].values.astype('float32')
    if n_top_genes is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_genes-1]
        gene_subset = df['dispersion_norm'].values >= disp_cut_off
        logg.msg('the {} top genes correspond to a normalized dispersion cutoff of'
                 .format(n_top_genes, disp_cut_off), v=5)
    else:
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce((mean > min_mean, mean < max_mean,
                                             dispersion_norm > min_disp,
                                             dispersion_norm < max_disp))
    logg.msg('    finished', time=True, v=4)
    return np.rec.fromarrays((gene_subset,
                              df['mean'].values,
                              df['dispersion'].values,
                              df['dispersion_norm'].values.astype('float32', copy=False)),
                              dtype=[('gene_subset', bool),
                                     ('means', 'float32'),
                                     ('dispersions', 'float32'),
                                     ('dispersions_norm', 'float32')])


def filter_genes_cv_deprecated(X, Ecutoff, cvFilter):
    """Filter genes by coefficient of variance and mean.

    See `filter_genes_dispersion`.

    Reference: Weinreb et al. (2017).
    """
    if issparse(X):
        raise ValueError('Not defined for sparse input. See `filter_genes_dispersion`.')
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.std(X, axis=0) / (np.mean(X, axis=0) + .0001) > cvFilter
    gene_subset = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_subset


def filter_genes_fano_deprecated(X, Ecutoff, Vcutoff):
    """Filter genes by fano factor and mean.

    See `filter_genes_dispersion`.

    Reference: Weinreb et al. (2017).
    """
    if issparse(X):
        raise ValueError('Not defined for sparse input. See `filter_genes_dispersion`.')
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.var(X, axis=0) / (np.mean(X, axis=0) + .0001) > Vcutoff
    gene_subset = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_subset


def log1p(data, copy=False, chunked=False, chunk_size=None):
    """Logarithmize the data matrix.

    Computes `X = log(X + 1)`, where `log` denotes the natural logrithm.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Returns or updates `data`, depending on `copy`.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        if chunked:
            for chunk, start, end in adata.chunked_X(chunk_size):
                adata.X[start:end] = log1p(chunk)
        else:
            adata.X = log1p(data.X)
        return adata if copy else None
    X = data  # proceed with data matrix
    if not issparse(X):
        return np.log1p(X)
    else:
        return X.log1p()


def pca(data, n_comps=None, zero_center=True, svd_solver='auto', random_state=0,
        return_info=False, dtype='float32', copy=False, chunked=False, chunk_size=None):
    """Principal component analysis [Pedregosa11]_.

    Computes PCA coordinates, loadings and variance decomposition. Uses the
    implementation of *scikit-learn* [Pedregosa11]_.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    n_comps : `int`, optional (default: 50)
        Number of principal components to compute.
    zero_center : `bool` or `None`, optional (default: `True`)
        If `True`, compute standard PCA from covariance matrix. If `False`, omit
        zero-centering variables (uses *TruncatedSVD* from scikit-learn), which
        allows to handle sparse input efficiently.
    svd_solver : `str`, optional (default: 'auto')
        SVD solver to use. Either 'arpack' for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or 'randomized' for the randomized algorithm
        due to Halko (2009). "auto" chooses automatically depending on the size
        of the problem.
    random_state : `int`, optional (default: 0)
        Change to use different intial states for the optimization.
    return_info : `bool`, optional (default: `False`)
        Only relevant when not passing an :class:`~scanpy.api.AnnData`: see
        "Returns".
    dtype : `str` (default: 'float32')
        Numpy data type string to which to convert the result.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    If `data` is array-like and `return_info == False`, only returns `X_pca`,\
    otherwise returns or adds to `adata`:
    X_pca : `.obsm`
         PCA representation of data.
    PCs : `.varm`
         The principal components containing the loadings.
    variance_ratio : `.uns['pca']`
         Ratio of explained variance.
    variance : `.uns['pca']`
         Explained variance, equivalent to the eigenvalues of the covariance matrix.
    """

    if n_comps is None: n_comps = N_PCS

    if isinstance(data, AnnData):
        data_is_AnnData = True
        adata = data.copy() if copy else data
    else:
        data_is_AnnData = False
        adata = AnnData(data)

    logg.msg('computing PCA with n_comps =', n_comps, r=True, v=4)

    if adata.n_vars < n_comps:
        n_comps = adata.n_vars - 1
        logg.msg('reducing number of computed PCs to',
               n_comps, 'as dim of data is only', adata.n_vars, v=4)

    if chunked:
        if not zero_center or random_state or svd_solver != 'auto':
            logg.msg('Ignoring zero_center, random_state, svd_solver', v=4)

        from sklearn.decomposition import IncrementalPCA

        X_pca = np.zeros((adata.X.shape[0], n_comps), adata.X.dtype)

        pca_ = IncrementalPCA(n_components=n_comps)

        for chunk, _, _ in adata.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            pca_.partial_fit(chunk)

        for chunk, start, end in adata.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            X_pca[start:end] = pca_.transform(chunk)
    else:
        zero_center = zero_center if zero_center is not None else False if issparse(adata.X) else True
        if zero_center:
            from sklearn.decomposition import PCA
            if issparse(adata.X):
                logg.msg('    as `zero_center=True`, '
                       'sparse input is densified and may '
                       'lead to huge memory consumption', v=4)
                X = adata.X.toarray()  # Copying the whole adata.X here, could cause memory problems
            else:
                X = adata.X
            pca_ = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
        else:
            from sklearn.decomposition import TruncatedSVD
            logg.msg('    without zero-centering: \n'
                   '    the explained variance does not correspond to the exact statistical defintion\n'
                   '    the first component, e.g., might be heavily influenced by different means\n'
                   '    the following components often resemble the exact PCA very closely', v=4)
            pca_ = TruncatedSVD(n_components=n_comps, random_state=random_state)
            X = adata.X
        X_pca = pca_.fit_transform(X)

    if X_pca.dtype.descr != np.dtype(dtype).descr: X_pca = X_pca.astype(dtype)

    if data_is_AnnData:
        adata.obsm['X_pca'] = X_pca
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


def normalize_per_cell(data, counts_per_cell_after=None, counts_per_cell=None,
                       key_n_counts=None, copy=False):
    """Normalize total counts per cell.

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
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
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

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
        cell_subset, counts_per_cell = filter_cells(adata.X, min_counts=1)
        adata.obs[key_n_counts] = counts_per_cell
        adata._inplace_subset_obs(cell_subset)
        normalize_per_cell(adata.X, counts_per_cell_after,
                           counts_per_cell=counts_per_cell[cell_subset])
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
        cell_subset, counts_per_cell = filter_cells(X, min_counts=1)
        X = X[cell_subset]
        counts_per_cell = counts_per_cell[cell_subset]
    if counts_per_cell_after is None:
        counts_per_cell_after = np.median(counts_per_cell)
    counts_per_cell /= counts_per_cell_after
    if not issparse(X): X /= counts_per_cell[:, np.newaxis]
    else: sparsefuncs.inplace_row_scale(X, 1/counts_per_cell)
    return X if copy else None


def normalize_per_cell_weinreb16_deprecated(X, max_fraction=1,
                                            mult_with_mean=False):
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
    X_norm : np.ndarray
        Normalized version of the original expression matrix.
    """
    if issparse(X):
        raise ValueError('Sparse input not allowed. '
                         'Consider `sc.pp.normalize_per_cell` instead.')
    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('Choose max_fraction between 0 and 1.')
    counts_per_cell = np.sum(X, axis=1)
    if max_fraction == 1:
        X_norm = X / counts_per_cell[:, np.newaxis]
        return X_norm
    # restrict computation of counts to genes that make up less than
    # constrain_theshold of the total reads
    tc_tiled = np.tile(counts_per_cell[:, np.newaxis], (1, X.shape[1]))
    included = np.all(X <= tc_tiled * max_fraction, axis=0)
    tc_include = np.sum(X[:, included], axis=1)
    tc_tiled = np.tile(tc_include[:, np.newaxis], (1, X.shape[1])) + 1e-6
    X_norm = X / tc_tiled
    if mult_with_mean:
        X_norm *= np.mean(counts_per_cell)
    return X_norm


def regress_out(adata, keys, n_jobs=None, copy=False):
    """Regress out unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R [Satija15].

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    keys : `str` or list of `str`
        Keys for observation annotation on which to regress on.
    n_jobs : `int` or `None`, optional (default: `None`)
        Number of jobs for parallel computation. Currently has no effect.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depening on `copy` returns or updates `adata` with the corrected data matrix.
    """
    logg.info('regressing out', keys, r=True)
    if issparse(adata.X):
        logg.info('    sparse input is densified and may '
                  'lead to high memory use')
    adata = adata.copy() if copy else adata
    if isinstance(keys, str): keys = [keys]
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    if n_jobs is not None:
        logg.warn('Parallelization is currently broke, will be restored soon. Running on 1 core.')
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    # regress on a single categorical variable
    sanitize_anndata(adata)
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
    # regress on one or several ordinal variables
    else:
        regressors = np.array(
            [adata.obs[key].values if key in adata.obs_keys()
             else adata[:, key].X for key in keys]).T
    regressors = np.c_[np.ones(adata.X.shape[0]), regressors]
    len_chunk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_chunks = np.ceil(adata.X.shape[1] / len_chunk).astype(int)
    chunks = [np.arange(start, min(start + len_chunk, adata.X.shape[1]))
              for start in range(0, n_chunks * len_chunk, len_chunk)]

    import statsmodels.api as sm
    from statsmodels.tools.sm_exceptions import PerfectSeparationError

    def _regress_out(col_index, responses, regressors):
        try:
            if regressors.shape[1] - 1 == responses.shape[1]:
                regressors_view = np.c_[regressors[:, 0], regressors[:, col_index + 1]]
            else:
                regressors_view = regressors
            result = sm.GLM(responses[:, col_index],
                            regressors_view, family=sm.families.Gaussian()).fit()
            new_column = result.resid_response
        except PerfectSeparationError:  # this emulates R's behavior
            logg.warn('Encountered PerfectSeparationError, setting to 0 as in R.')
            new_column = np.zeros(responses.shape[0])
        return new_column

    def _regress_out_chunk(chunk, responses, regressors):
        chunk_array = np.zeros((responses.shape[0], chunk.size),
                               dtype=responses.dtype)
        for i, col_index in enumerate(chunk):
            chunk_array[:, i] = _regress_out(col_index, responses, regressors)
        return chunk_array

    for chunk in chunks:
        # why did this break after migrating to dataframes?
        # result_lst = Parallel(n_jobs=n_jobs)(
        #     delayed(_regress_out)(
        #         col_index, adata.X, regressors) for col_index in chunk)
        result_lst = [_regress_out(
            col_index, adata.X, regressors) for col_index in chunk]
        for i_column, column in enumerate(chunk):
            adata.X[:, column] = result_lst[i_column]
    logg.info('    finished', t=True)
    return adata if copy else None


def scale(data, zero_center=True, max_value=None, copy=False):
    """Scale data to unit variance and zero mean.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    zero_center : `bool`, optional (default: `True`)
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value : `float` or `None`, optional (default: `None`)
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
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


def subsample(data, fraction, random_state=0, copy=False):
    """Subsample to a fraction of the number of observations.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    fraction : float in [0, 1]
        Subsample to this `fraction` of the number of observations.
    random_state : `int` or `None`, optional (default: 0)
        Random seed to change subsampling.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~scanpy.api.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """
    if fraction > 1 or fraction < 0:
        raise ValueError('`fraction` needs to be within [0, 1], not {}'
                         .format(fraction))
    np.random.seed(random_state)
    n_obs = data.n_obs if isinstance(data, AnnData) else data.shape[0]
    new_n_obs = int(fraction * n_obs)
    obs_indices = np.random.choice(n_obs, size=new_n_obs, replace=False)
    logg.msg('... subsampled to {} data points'.format(new_n_obs))
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata._inplace_subset_obs(obs_indices)
        return adata if copy else None
    else:
        X = data
        return X[obs_indices], obs_indices


def downsample_counts(adata, target_counts=20000, random_state=0, copy=False):
    """Downsample counts so that each cell has no more than `target_counts`.

    Cells with fewer counts than `target_counts` are unaffected by this. This
    has been implemented by M. D. Luecken.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    target_counts : `int` (default: 20,000)
        Target number of counts for downsampling. Cells with more counts than
        'target_counts' will be downsampled to have 'target_counts' counts.
    random_state : `int` or `None`, optional (default: 0)
        Random seed to change subsampling.
    copy : `bool`, optional (default: `False`)
        If an :class:`~scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates an `adata` with downsampled `.X`.
    """
    if target_counts < 1:
        raise ValueError('`target_counts` must be a positive integer'
                         .format(target_counts))
    if not isinstance(adata, AnnData):
        raise ValueError('`adata` must be an `AnnData` object'.format(adata))
    logg.msg('downsampling to {} counts'.format(target_counts), r=True)
    adata = adata.copy() if copy else adata
    np.random.seed(random_state)
    counts = adata.X.sum(axis=1)
    adata.obs['n_counts'] = counts
    for icell, _ in enumerate(adata.obs_names):
        if counts[icell] > target_counts:
            idx_vec = []
            for ix, i in enumerate(adata.X[icell].astype(int)):
                idx_vec.extend([ix]*i)
            # idx_vec = np.array(idx_vec)
            downsamp = np.random.choice(idx_vec, target_counts)
            cell_profile = np.zeros(adata.n_vars)
            indices, values = np.unique(downsamp, return_counts=True)
            for i in range(len(indices)):
                cell_profile[indices[i]] = values[i]
            adata.X[icell] = cell_profile
    logg.msg('finished', t=True)
    return adata if copy else None


def zscore_deprecated(X):
    """Z-score standardize each variable/gene in X.

    Use `scale` instead.

    Reference: Weinreb et al. (2017).

    Parameters
    ----------
    X : np.ndarray
        Data matrix. Rows correspond to cells and columns to genes.

    Returns
    -------
    XZ : np.ndarray
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


def _get_mean_var(X):
    # - using sklearn.StandardScaler throws an error related to
    #   int to long trafo for very large matrices
    # - using X.multiply is slower
    if True:
        mean = X.mean(axis=0)
        if issparse(X):
            mean_sq = X.multiply(X).mean(axis=0)
            mean = mean.A1
            mean_sq = mean_sq.A1
        else:
            mean_sq = np.multiply(X, X).mean(axis=0)
        # enforece R convention (unbiased estimator) for variance
        var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=False).partial_fit(X)
        mean = scaler.mean_
        # enforce R convention (unbiased estimator)
        var = scaler.var_ * (X.shape[0]/(X.shape[0]-1))
    return mean, var


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
            X /= scale
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=zero_center, copy=False).partial_fit(X)
        # user R convention (unbiased estimator)
        scaler.scale_ *= np.sqrt(X.shape[0]/(X.shape[0]-1))
        scaler.transform(X)

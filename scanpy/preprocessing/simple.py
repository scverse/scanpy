"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""

import numpy as np
import scipy as sp
import warnings
from joblib import Parallel, delayed
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from pandas.api.types import is_categorical_dtype
from anndata import AnnData
from .. import settings as sett
from .. import logging as logg


def filter_cells(data, min_counts=None, min_genes=None, max_counts=None,
                 max_genes=None, copy=False):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers, i.e.,
    "unreliable" observations.

    Only provide one of the optional arguments (`min_counts`, `min_genes`,
    `max_counts`, `max_genes`) per call.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.spmatrix`
        Data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and
        columns to genes.
    min_counts : `int`, optional (default: `None`)
        Minimum number of counts required for a cell to pass filtering.
    min_genes : `int`, optional (default: `None`)
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts : `int`, optional (default: `None`)
        Maximum number of counts required for a cell to pass filtering.
    max_genes : `int`, optional (default: `None`)
        Maximum number of genes expressed required for a cell to pass filtering.
    copy : `bool`, optional (default: `False`)
        If an :class:`scanpy.api.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    If `data` is an :class:`~scanpy.api.AnnData`, filters the object and adds
    either `n_genes` or `n_counts` to `adata.obs`. Otherwise a tuple

    cell_subset : `np.ndarray`
        Boolean index mask that does filtering. `True` means that the cell is
        kept. `False` means the cell is removed.
    number_per_cell: `np.ndarray`
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
    logg.m('filtered out {} cells that have'.format(s), end=' ', v=4)
    if min_genes is not None or min_counts is not None:
        logg.m('less than',
               str(min_genes) + ' genes expressed'
               if min_counts is None else str(min_counts) + ' counts', v=4, no_indent=True)
    if max_genes is not None or max_counts is not None:
        logg.m('more than ',
               str(max_genes) + ' genes expressed'
               if max_counts is None else str(max_counts) + ' counts', v=4, no_indent=True)
    return cell_subset, number_per_cell


def filter_genes(data, min_cells=None, min_counts=None, copy=False):
    """Filter genes based on minimal number of cells or counts.

    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells.

    See :func:`~scanpy.api.pp.filter_cells`.
    """
    if min_cells is not None and min_counts is not None:
        raise ValueError('Either specify min_counts or min_cells, but not both.')
    if min_cells is None and min_counts is None:
        raise ValueError('Provide one of min_counts or min_cells.')
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        gene_subset, number = filter_genes(adata.X, min_cells=min_cells,
                                           min_counts=min_counts)
        if min_cells is None:
            adata.var['n_counts'] = number
        else:
            adata.var['n_cells'] = number
        adata._inplace_subset_var(gene_subset)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    number_per_gene = np.sum(X if min_cells is None else X > 0, axis=0)
    min_number = min_counts if min_cells is None else min_cells
    if issparse(X):
        number_per_gene = number_per_gene.A1
    gene_subset = number_per_gene >= min_number
    logg.m('filtered out', np.sum(~gene_subset),
           'genes that are detected',
           'in less than ' + str(min_cells) + ' cells' if min_counts is None
           else 'with less than ' + str(min_counts) + ' counts', v=4)
    return gene_subset, number_per_gene


def filter_genes_dispersion(data,
                            flavor='seurat',
                            min_disp=None, max_disp=None,
                            min_mean=None, max_mean=None,
                            n_top_genes=None,
                            log=True,
                            copy=False):
    """Filter genes based on dispersion: extract highly variable genes.

    If trying out parameters, pass the data matrix instead of AnnData.

    Similar functions are used, for example, by Cell Ranger [Zheng17]_
    and Seurat [Satija15]_.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.sparse`
        Data matrix.
    flavor : {'seurat', 'cell_ranger'}, optional (default: 'seurat')
        Choose method for computing normalized dispersion. If choosing 'Seurat',
        this expects non-logarithmized data, you can change this by setting
        `log` to `False`. Note that Seurat passes the cutoffs whereas Cell
        Ranger passes `n_top_genes`.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=`None` : `float`, optional
        If `n_top_genes` is not `None`, these cutoffs for the normalized gene
        expression are ignored.
    n_top_genes : `int` or `None` (default: `None`)
        Number of highly-variable genes to keep.
    log : `bool`, optional (default: True)
        Use the logarithm of mean and variance.
    copy : `bool`, optional (default: `False`)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    If an AnnData `adata` is passed, returns or updates `adata` depending on
    `copy`. It filters the adata object and adds the annotations

    means : pd.Series (adata.var)
        Means per gene.
    dispersions : pd.Series (adata.var)
        Dispersions per gene.
    dispersions_norm : pd.Series (adata.var)
        Normalized dispersions per gene.

    If a data matrix `X` is passed, the annotation is returned as `np.recarray`
    with the columns:
        gene_subset, means, dispersions, dispersion_norm
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
    logg.info('filter highly variable genes by dispersion and mean',
              r=True, end=' ')
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
        df['mean_bin'] = pd.cut(df['mean'], bins=20)
        disp_grouped = df.groupby('mean_bin')['dispersion']
        disp_mean_bin = disp_grouped.mean()
        disp_std_bin = disp_grouped.std(ddof=1)
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
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_genes-1]
        gene_subset = df['dispersion_norm'].values >= disp_cut_off
        logg.m(t=True)
        logg.m('the', n_top_genes,
               'top genes correspond to a normalized dispersion cutoff of',
               disp_cut_off, v=4)
    else:
        logg.m(t=True, no_indent=True)
        logg.m('using `min_disp={}`, `max_disp={}`, `min_mean={}` and `max_mean={}`'
               .format(min_disp, max_disp, min_mean, max_mean), v=4)
        logg.hint('set `n_top_genes` to simply select top-scoring genes instead')
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce((mean > min_mean, mean < max_mean,
                                             dispersion_norm > min_disp,
                                             dispersion_norm < max_disp))
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


def log1p(data, copy=False):
    """Logarithmize the data matrix.

    Computes `X = log(X + 1)`, where `log` denotes the natural logrithm.

    Parameters
    ----------
    data : array-like or AnnData
        The data matrix.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    Returns or updates data, depending on `copy`.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata.X = log1p(data.X)
        return adata if copy else None
    X = data  # proceed with data matrix
    if not issparse(X):
        return np.log1p(X)
    else:
        return X.log1p()


def pca(data, n_comps=50, zero_center=True, svd_solver='auto', random_state=0,
        recompute=True, mute=False, return_info=None, copy=False,
        dtype='float32'):
    """Principal component analysis [Pedregosa11]_.

    Computes PCA coordinates, loadings and variance decomposition. Uses the
    implementation of *scikit-learn* [Pedregosa11]_.

    Parameters
    ----------
    data : AnnData, array-like
        Data matrix of shape n_obs × n_vars.
    n_comps : int, optional (default: 10)
        Number of principal components to compute.
    zero_center : bool or None, optional (default: None)
        If True, compute standard PCA from Covariance matrix. If False, omit
        zero-centering variables, which allows to handle sparse input
        efficiently. If None, defaults to True for dense and to False for sparse
        input.
    svd_solver : str, optional (default: 'auto')
        SVD solver to use. Either 'arpack' for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or 'randomized' for the randomized algorithm
        due to Halko (2009). "auto" chooses automatically depending on the size
        of the problem.
    random_state : int, optional (default: 0)
        Change to use different intial states for the optimization.
    recompute : bool, optional (default: True)
        Use the result of previous calculation, if possible.
    return_info : bool or None, optional (default: None)
        If providing an array, this defaults to False, if providing an AnnData,
        defaults to true.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    dtype : str (default: 'float32')
        Numpy data type string to which to convert the result.

    Returns
    -------
    If X is array-like and ``return_info == True``, only returns ``X_pca``, otherwise adds to ``adata``:
    X_pca : np.ndarray (adata.obs)
         PCA representation of the data with shape n_variables × n_comps.
    components / PC1, PC2, PC3, ... : np.ndarray (adata.var)
         The PCs containing the loadings as shape n_comps × n_vars.
    variance_ratio : np.ndarray (adata.uns)
         Ratio of explained variance.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        from .. import settings as sett  # why is this necessary?
        if ('X_pca' in adata.obs
            and adata.obsm['X_pca'].shape[1] >= n_comps
            and not recompute
            and (sett.recompute == 'none' or sett.recompute == 'pp')):
            logg.m('    not recomputing PCA, using "X_pca" contained '
                   'in `adata.obs` (set `recompute=True` to avoid this)', v=4)
            return adata
        else:
            logg.m('compute PCA with n_comps =', n_comps, r=True, v=4)
            result = pca(adata.X, n_comps=n_comps, zero_center=zero_center,
                         svd_solver=svd_solver, random_state=random_state,
                         recompute=recompute, mute=mute, return_info=True)
            X_pca, components, pca_variance_ratio = result
            adata.obsm['X_pca'] = X_pca
            adata.varm['PCs'] = components.T
            adata.uns['pca_variance_ratio'] = pca_variance_ratio
            logg.m('    finished', t=True, end=' ', v=4)
            logg.m('and added\n'
                      '    "X_pca", the PCA coordinates (adata.obs)\n'
                      '    "PC1", "PC2", ..., the loadings (adata.var)\n'
                      '    "pca_variance_ratio", the variance ratio (adata.uns)', v=4)
        return adata if copy else None
    X = data  # proceed with data matrix
    from .. import settings as sett
    if X.shape[1] < n_comps:
        n_comps = X.shape[1] - 1
        logg.m('reducing number of computed PCs to',
               n_comps, 'as dim of data is only', X.shape[1], v=4)
    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
    from sklearn.decomposition import PCA, TruncatedSVD
    verbosity_level = np.inf if mute else 0
    if zero_center:
        if issparse(X):
            logg.m('    as `zero_center=True`, '
                   'sparse input is densified and may '
                   'lead to huge memory consumption', v=4)
            X = X.toarray()
        pca_ = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
    else:
        logg.m('    without zero-centering: \n'
               '    the explained variance does not correspond to the exact statistical defintion\n'
               '    the first component, e.g., might be heavily influenced by different means\n'
               '    the following components often resemble the exact PCA very closely', v=4)
        pca_ = TruncatedSVD(n_components=n_comps, random_state=random_state)
    X_pca = pca_.fit_transform(X)
    if X_pca.dtype.descr != np.dtype(dtype).descr: X_pca = X_pca.astype(dtype)
    if False if return_info is None else return_info:
        return X_pca, pca_.components_, pca_.explained_variance_ratio_
    else:
        return X_pca


def normalize_per_cell(data, counts_per_cell_after=None, counts_per_cell=None,
                       key_n_counts=None, copy=False):
    """Normalize each cell.

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~scanpy.api.AnnData`, `np.ndarray`, `sp.spmatrix`
        Data matrix. Rows correspond to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    key_n_counts : str, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy : `bool` (default: `False`)
        Determines whether function operates inplace (default) or a copy is
        returned.

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
        logg.info('normalizing by total count per cell', r=True)
        adata = data.copy() if copy else data
        cell_subset, counts_per_cell = filter_cells(adata.X, min_counts=1)
        adata.obs[key_n_counts] = counts_per_cell
        adata._inplace_subset_obs(cell_subset)
        normalize_per_cell(adata.X, counts_per_cell_after,
                           counts_per_cell=counts_per_cell[cell_subset])
        logg.info('    finished', t=True, end=': ')
        logg.info('normalized adata.X and added', no_indent=True)
        logg.info('    \'{}\', counts per cell before normalization (adata.obs)'
               .format(key_n_counts),
               no_indent=True)
        return adata if copy else None
    # proceed with data matrix
    X = data.copy() if copy else data
    if counts_per_cell is None:
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
    """Normalize each cell.

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
    adata : AnnData
        The annotated data matrix.
    keys : str or list of strings
        Keys for sample annotation on which to regress on.
    n_jobs : int
        Number of jobs for parallel computation.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    Depening on `copy` returns or updates `adata` with the corrected data matrix.
    """
    logg.info('regressing out', keys, r=True)
    if issparse(adata.X):
        logg.info('... sparse input is densified and may '
                  'lead to huge memory consumption')
    adata = adata.copy() if copy else adata
    if isinstance(keys, str): keys = [keys]
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    if n_jobs is not None:
        logg.warn('Parallelization is currently broke, will be restored soon. Running on 1 core.')
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    # regress on a single categorical variable
    if keys[0] in adata.obs_keys() and is_categorical_dtype(adata.obs[keys[0]]):
        if len(keys) > 1:
            raise ValueError(
                'If providing categorical variable, '
                'only a single one is allowed. For this one '
                'the mean is computed for each variable/gene.')
        logg.m('... regressing on per-gene means within categories', v=4)
        unique_categories = np.unique(adata.obs[keys[0]].values)
        regressors = np.zeros(adata.X.shape, dtype='float32')
        for category in unique_categories:
            mask = category == adata.obs[keys[0]].values
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
    logg.info('finished', t=True)
    logg.hint('after `sc.pp.regress_out`, consider rescaling the adata using `sc.pp.scale`')
    return adata if copy else None


def scale(data, zero_center=True, max_value=None, copy=False):
    """Scale data to unit variance and zero mean.

    Parameters
    ----------
    zero_center : `bool`, optional (default: `True`)
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value : `float` or `None`, optional (default: `None`)
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy : `bool` (default: `False`)
        Perfrom operation inplace if `False`.

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
                'densified and may lead to large memory consumption', v=4)
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
        logg.msg('... clipping at max_value', max_value, v=4)
    if zero_center and issparse(X):
        logg.msg('... scale_data: as `zero_center=True`, sparse input is '
                 'densified and may lead to large memory consumption, returning copy',
                 v=4)
        X = X.toarray()
        copy = True
    _scale(X, zero_center)
    if max_value is not None: X[X > max_value] = max_value
    return X if copy else None


def subsample(data, fraction, seed=0, simply_skip_samples=False, copy=False):
    """Subsample to a fraction of the number of samples.

    Parameters
    ----------
    data : AnnData or array-like
        Annotated data matrix.
    fraction : float in [0, 1]
        Subsample to this `fraction` of the number of samples.
    seed : int, optional (default: 0)
        Random seed to change subsampling.
    simply_skip_samples : bool, optional (default: False)
        Simply skip samples instead of true sampling.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    Updates or returns the subsampled data, depending on `copy`. Returns
    ``X, obs_indices`` if data is array-like, otherwise subsamples the passed
    `AnnData` (``copy == False``) or a copy of it (``copy == True``).
    """
    if fraction > 1 or fraction < 0:
        raise ValueError('`fraction` needs to be within [0, 1], not {}'
                         .format(fraction))
    np.random.seed(seed)
    n_obs = data.n_obs if isinstance(data, AnnData) else data.shape[0]
    new_n_obs = int(fraction * n_obs)
    if simply_skip_samples:
        obs_indices = np.arange(0, n_obs, int(1./fraction), dtype=int)
    else:
        obs_indices = np.random.choice(n_obs, size=new_n_obs, replace=False)
    logg.m('... subsampled to {} data points'.format(new_n_obs), v=4)
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata._inplace_subset_obs(obs_indices)
        return adata if copy else None
    else:
        X = data
        return X[obs_indices], obs_indices


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

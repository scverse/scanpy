# Author: F. Alex Wolf (http://falexwolf.de)
"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""

import numpy as np
import scipy as sp
import warnings
from joblib import Parallel, delayed
from scipy.sparse import issparse
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from sklearn.utils import sparsefuncs
from ..data_structs import AnnData
from .. import settings as sett
from .. import logging as logg


def filter_cells(data, min_counts=None, min_genes=None, max_counts=None,
                 max_genes=None, copy=False):
    """Filter outliers based on counts and number of genes expressed.

    For instance, only keep cells with at least `min_counts` UMI counts or
    `min_genes` genes expressed.

    Only provide one of the optional arguments per call.

    This is to filter measurement outliers, i.e., "unreliable" samples.

    Paramaters
    ----------
    data : np.ndarray or AnnData
        Data matrix of shape n_sample x n_variables. Rows correspond to cells
        and columns to genes.
    min_counts : int
        Minimum number of counts required for a cell to pass filtering.
    min_genes : int
        Minimum number of genes expressed required for a cell to pass filtering.
    min_counts : int
        Maximum number of counts required for a cell to pass filtering.
    max_genes : int
        Maximum number of genes expressed required for a cell to pass filtering.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    If data is a data matrix X, the following to arrays are returned
        cell_subset : np.ndarray
            Boolean index mask that does filtering. True means that the cell is
            kept. False means the cell is removed.
        number_per_cell: np.ndarray
            Either n_counts or n_genes per cell.
    otherwise, if copy == False, the adata object is updated, otherwise a copy is returned
        adata : AnnData
            The filtered adata object, with the count info stored in adata.smp.
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
        if min_genes is None and max_genes is None: adata.smp['n_counts'] = number
        else: adata.smp['n_genes'] = number
        adata.inplace_subset_smp(cell_subset)
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
    logg.m('... filtered out {} cells that have'.format(s), end=' ', v=4)
    if min_genes is not None or min_counts is not None:
        logg.m('less than',
               str(min_genes) + ' genes expressed'
               if min_counts is None else str(min_counts) + ' counts', v=4)
    if max_genes is not None or max_counts is not None:
        logg.m('more than ',
               str(max_genes) + ' genes expressed'
               if max_counts is None else str(max_counts) + ' counts', v=4)
    return cell_subset, number_per_cell


def filter_genes(data, min_cells=None, min_counts=None, copy=False):
    """Filter genes.

    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells / samples. The latter means that genes are filtered
    whose expression is not supported by sufficient statistical evidence.

    Parameters
    ----------
    See filter_cells.
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
        adata.inplace_subset_var(gene_subset)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    number_per_gene = np.sum(X if min_cells is None else X > 0, axis=0)
    min_number = min_counts if min_cells is None else min_cells
    if issparse(X):
        number_per_gene = number_per_gene.A1
    gene_subset = number_per_gene >= min_number
    logg.m('... filtered out', np.sum(~gene_subset),
           'genes that are detected',
           'in less than ' + str(min_cells) + ' cells' if min_counts is None
           else 'with less than ' + str(min_counts) + ' counts', v=4)
    return gene_subset, number_per_gene


def filter_genes_dispersion(data,
                            flavor='seurat',
                            min_disp=0.5, max_disp=None,
                            min_mean=0.0125, max_mean=3,
                            n_top_genes=None,
                            log=True,
                            copy=False):
    """Extract highly variable genes.

    If trying out parameters, pass the data matrix instead of AnnData.

    Similar functions are used, for example, by Cell Ranger (Zheng et al., 2017)
    and Seurat (Macosko et al., 2015).

    Parameters
    ----------
    X : AnnData or array-like
        Data matrix storing unlogarithmized data.
    flavor : {'seurat', 'cell_ranger'}
        Choose method for computing normalized dispersion. Note that Seurat
        passes the cutoffs whereas Cell Ranger passes `n_top_genes`.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=None : float
        Cutoffs for the gene expression, used if n_top_genes is None.
    n_top_genes : int or None (default: None)
        Number of highly-variable genes to keep.
    log : bool
        Use the logarithm of mean and variance.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    If an AnnData `adata` is passed, returns or updates `adata` depending on
    `copy`. It filters the adata object and adds the annotations
        "means",  means per gene (adata.var)
        "dispersions", dispersions per gene (adata.var)
        "dispersions_norm", dispersions per gene (adata.var)
    If a data matrix `X` is passed, the annotation is returned as np.recarray
    with the columns:
        gene_subset, means, dispersions, dispersion_norm
    """
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
        adata.inplace_subset_var(result['gene_subset'])
        return adata if copy else None
    logg.m('... filter highly varying genes by dispersion and mean', r=True, end=' ', v=4)
    X = data  # proceed with data matrix
    mean, var = _get_mean_var(X)
    # now actually compute the dispersion
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
        logg.m(t=True, v=4)
        logg.m('    the', n_top_genes,
               'top genes correspond to a normalized dispersion cutoff of',
               disp_cut_off, v=4)
    else:
        logg.m(t=True, v=4)
        logg.m('    using `min_disp={}`, `max_disp={}`, `min_mean={}` and `max_mean={}`'
               .format(min_disp, max_disp, min_mean, max_mean), v=4)
        logg.m('--> set `n_top_genes` to simply select top-scoring genes instead', v=4)
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
    """Apply logarithm to count data "+1", i.e., `adata.X+1`.

    Parameters
    ----------
    data : array-like or AnnData
        The data matrix.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.
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


def pca(data, n_comps=50, zero_center=True, svd_solver='auto',
        random_state=0, recompute=True, mute=False, return_info=None, copy=False, dtype='float32'):
    """Embed data using PCA.

    Parameters
    ----------
    data : AnnData or array_like
        X : np.ndarray
            Data matrix of shape n_samples x n_variables.
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
    dtype : str
        Numpy data type string to which to convert the result.

    Notes
    -----
    If X is array-like and return_info == True, returns, otherwise adds to (a copy
    of) AnnData:
    X_pca : np.ndarray
         PCA representation of the data with shape n_variables x n_comps.
         Depending on whether an AnnData or a data matrix has been
         provided, the array is written to AnnData or returned directly.
    components / PC1, PC2, PC3, ... : np.ndarray
         The PCs containing the loadings as shape n_comps x n_vars. Written to
         adata.var if an AnnData object is provided.
    variance_ratio : np.ndarray
         Ratio of explained variance. Written as unstructured annotation to
         adata, if provided.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        from .. import settings as sett  # why is this necessary?
        if ('X_pca' in adata.smp
            and adata.smp['X_pca'].shape[1] >= n_comps
            and not recompute
            and (sett.recompute == 'none' or sett.recompute == 'pp')):
            logg.m('    not recomputing PCA, using "X_pca" contained '
                   'in `adata.smp` (set `recompute=True` to avoid this)', v=4)
            return adata
        else:
            logg.m('compute PCA with n_comps =', n_comps, r=True, v=4)
            result = pca(adata.X, n_comps=n_comps, zero_center=zero_center,
                         svd_solver=svd_solver, random_state=random_state,
                         recompute=recompute, mute=mute, return_info=True)
            X_pca, components, pca_variance_ratio = result
            adata.smp['X_pca'] = X_pca  # this is multicolumn-sample annotation
            for icomp, comp in enumerate(components):
                adata.var['PC' + str(icomp+1)] = comp
            adata.add['pca_variance_ratio'] = pca_variance_ratio
            logg.m('    finished', t=True, end=' ', v=4)
            logg.m('and added\n'
                      '    "X_pca", the PCA coordinates (adata.smp)\n'
                      '    "PC1", "PC2", ..., the loadings (adata.var)\n'
                      '    "pca_variance_ratio", the variance ratio (adata.add)', v=4)
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


def normalize_per_cell(data, counts_per_cell_after=None, copy=False,
                       counts_per_cell=None):
    """Normalize each cell.

    Normalize each cell by UMI count, so that every cell has the same total
    count.

    Similar functions are used, for example, by Cell Ranger (Zheng et al.,
    2017), Seurat (Macosko et al., 2015), Haghverdi et al. (2016) or Weinreb et
    al. (2016).

    Parameters
    ----------
    data : array_like, sparse or AnnData
        Data matrix. Rows correspond to cells and columns to genes.
    counts_per_cell_after : float or None (default: None)
        If None, after normalization, each cell has a total count equal
        to the median of the counts_per_cell before normalization.
    counts_per_cell : array (default: None)
        Precomputed counts per cell.
    copy : bool (default: False)
        Determines whether function operates inplace (default) or a copy is
        returned.

    Returns
    -------
    Returns or updates adata with normalized version of the original adata.X,
    depending on `copy`.
    """
    if isinstance(data, AnnData):
        logg.m('... normalizing by total count per cell', r=True, v=4)
        adata = data.copy() if copy else data
        cell_subset, counts_per_cell = filter_cells(adata.X, min_counts=1)
        adata.smp['n_counts'] = counts_per_cell
        adata.inplace_subset_smp(cell_subset)
        normalize_per_cell(adata.X, counts_per_cell_after, copy,
                           counts_per_cell=counts_per_cell[cell_subset])
        logg.m('    finished', t=True, v=4)
        logg.m('    normalized adata.X and added, '
                  '"n_counts", counts per cell before normalization (adata.smp)', v=4)
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

    Yields a dense matrix.

    This is inspired by Seurat's `regressOut` function in R.

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
    """
    logg.m('regressing out', keys, r=True, v=4)
    if issparse(adata.X):
        logg.m('... sparse input is densified and may '
                  'lead to huge memory consumption', v=4)
    if not copy:
        logg.m('--> note that this is an inplace computation '
               'and will return None: set `copy=True` if you want a copy', v=4)
    adata = adata.copy() if copy else adata
    if isinstance(keys, str): keys = [keys]
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    # regress on categorical variable
    if adata.smp[keys[0]].dtype.char == np.dtype('U').char:
        if len(keys) > 1:
            raise ValueError(
                'If providing categorical variable, '
                'only a single one is allowed. For this one '
                'the mean is computed for each variable/gene.')
        logg.m('... regressing on per-gene means within categories', v=4)
        unique_categories = np.unique(adata.smp[keys[0]])
        regressors = np.zeros(adata.X.shape, dtype='float32')
        for category in unique_categories:
            mask = category == adata.smp[keys[0]]
            for ix, x in enumerate(adata.X.T):
                regressors[mask, ix] = x[mask].mean()
    # regress on one or several ordinal variables
    else:
        regressors = np.array([adata.smp[key] for key in keys]).T
    regressors = np.c_[np.ones(adata.X.shape[0]), regressors]
    len_chunk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_chunks = np.ceil(adata.X.shape[1] / len_chunk).astype(int)
    chunks = [np.arange(start, min(start + len_chunk, adata.X.shape[1]))
              for start in range(0, n_chunks * len_chunk, len_chunk)]
    if sett.is_run_from_ipython:
        from tqdm import tqdm_notebook as tqdm
        # does not work in Rodeo, should be solved sometime soon
        # tqdm = lambda x: x
        # from tqdm import tqdm
        # sett.m(1, 'TODO: get nice waitbars also in interactive mode')
        # logg.m('... nicer progress bars as on command line come soon')
    else:
        from tqdm import tqdm
    for chunk in tqdm(chunks):
        result_lst = Parallel(n_jobs=n_jobs)(
            delayed(_regress_out)(
                col_index, adata.X, regressors) for col_index in chunk)
        for i_column, column in enumerate(chunk):
            adata.X[:, column] = result_lst[i_column]
    logg.m('finished', t=True, v=4)
    logg.m('--> after `sc.pp.regress_out`, consider rescaling the adata using `sc.pp.scale`', v=4)
    return adata if copy else None


def scale(data, zero_center=True, max_value=None, copy=False):
    """Scale data to unit variance and zero mean (`if zero_center`).

    Parameters
    ----------
    zero_center : bool or None, optional (default: None)
        If False, omit zero-centering variables, which allows to handle sparse
        input efficiently. If None, for dense input, defaults to True, for
        sparse input, defaults to False.
    max_value : None or float, optional (default: None)
        Clip to this value after scaling. If None, do not clip.
    copy : bool (default: False)
        Perfrom operation inplace if False.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        # need to add the following here to make inplace logic work
        if zero_center and issparse(adata.X):
            logg.m('... scale_data: as `zero_center=True`, sparse input is '
                   'densified and may lead to large memory consumption', v=4)
            adata.X = adata.X.toarray()
        scale(adata.X, zero_center=zero_center, max_value=max_value, copy=copy)
        return adata if copy else None
    X = data.copy() if copy else data  # proceed with the data matrix
    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
    if not zero_center and max_value is not None:
        logg.m('... scale_data: be very careful to use `max_value` without `zero_center`', v=4)
    if max_value is not None:
        logg.m('... clipping at max_value', max_value, v=4)
    if zero_center and issparse(X):
        logg.m('... scale_data: as `zero_center=True`, sparse input is '
               'densified and may lead to large memory consumption, returning copy', v=4)
        X = X.toarray()
        copy = True
    _scale(X, zero_center)
    if max_value is not None: X[X > max_value] = max_value
    return X if copy else None


def subsample(data, fraction, seed=0, simply_skip_samples=False, copy=False):
    """Subsample to `fraction` of the data.

    Parameters
    ----------
    data : AnnData or array-like
        Annotated data matrix.
    fraction : float in [0, 1]
        Subsample to a fraction the number of samples.
    seed : int
        Random seed to change subsampling.
    simply_skip_samples : bool, optional (default: False)
        Simply skip samples instead of true sampling.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    Updates or returns the subsampled data object, depending on `copy`.

    Notes
    -----
    Returns X, smp_indices if data is array-like, otherwise subsamples the passed
    AnnData (copy == False) or a copy of it (copy == True).
    """
    if fraction > 1 or fraction < 0:
        raise ValueError('`fraction` needs to be within [0, 1], not {}'
                         .format(fraction))
    np.random.seed(seed)
    n_smps = data.n_smps if isinstance(data, AnnData) else data.shape[0]
    new_n_smps = int(fraction * n_smps)
    if simply_skip_samples:
        smp_indices = np.arange(0, n_smps, int(1./fraction), dtype=int)
    else:
        smp_indices = np.random.choice(n_smps, size=new_n_smps, replace=False)
    logg.m('... subsampled to {} data points'.format(new_n_smps), v=4)
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata.inplace_subset_smp(smp_indices)
        return adata if copy else None
    else:
        X = data
        return X[smp_indices], smp_indices


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
        logg.warn('Encountered PerfectSeparationError, setting to zero as in R.')
        new_column = np.zeros(responses.shape[0])
    return new_column


def _regress_out_chunk(chunk, responses, regressors):
    chunk_array = np.zeros((responses.shape[0], chunk.size),
                           dtype=responses.dtype)
    for i, col_index in enumerate(chunk):
        chunk_array[:, i] = _regress_out(col_index, responses, regressors)
    return chunk_array


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

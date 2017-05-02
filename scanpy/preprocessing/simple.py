# Author: F. Alex Wolf (http://falexwolf.de)
"""Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""

import numpy as np
import scipy as sp
from joblib import Parallel, delayed
from scipy.sparse import issparse
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from ..data_structs import AnnData
from .. import sett


def filter_cells(data, min_counts=None, min_genes=None, copy=False):
    """
    Keep cells that have at least `min_counts` UMI counts or `min_genes`
    genes expressed.

    This is to filter measurement outliers, i.e., "unreliable" samples.

    Paramaters
    ----------
    data : np.ndarray or AnnData
        Data matrix of shape n_sample x n_variables. Rows correspond to cells
        and columns to genes.
    min_counts : int
        Minimum number of counts required for a cell to pass filtering.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Notes
    -----
    If data is a data matrix X, the following to arrays are returned
        cell_filter : np.ndarray
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
    if min_genes is None and min_counts is None:
        raise ValueError('Provide one of min_counts or min_genes.')
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        cell_filter, number = filter_cells(adata.X, min_counts, min_genes)
        if min_genes is None:
            adata.smp['n_counts'] = number
        else:
            adata.smp['n_genes'] = number
        adata.filter_smp(cell_filter)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    number_per_cell = np.sum(X if min_genes is None else X > 0, axis=1)
    if issparse(X):
        number_per_cell = number_per_cell.A1
    cell_filter = number_per_cell >= min_number
    sett.m(0, '... filtered out', np.sum(~cell_filter), 'outlier cells')
    return cell_filter, number_per_cell


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
        gene_filter, number = filter_genes(adata.X, min_cells=min_cells,
                                           min_counts=min_counts)
        if min_cells is None:
            adata.var['n_counts'] = number
        else:
            adata.var['n_cells'] = number
        adata.filter_var(gene_filter)
        return adata if copy else None
    X = data  # proceed with processing the data matrix
    number_per_gene = np.sum(X if min_cells is None else X > 0, axis=0)
    min_number = min_counts if min_cells is None else min_cells
    if issparse(X):
        number_per_gene = number_per_gene.A1
    gene_filter = number_per_gene >= min_number
    sett.m(0, '... filtered out', np.sum(~gene_filter),
           'genes that are detected',
           'in less than ' + str(min_cells) + ' cells' if min_counts is None
           else 'with less than ' + str(min_counts) + ' counts')
    return gene_filter, number_per_gene


def filter_genes_dispersion(data, log=True,
                            min_disp=0.5, max_disp=None,
                            min_mean=0.0125, max_mean=3,
                            n_top_genes=None,
                            flavor='seurat',
                            plot=False, copy=False):
    """Extract highly variable genes.

    Similar functions are used, for example, by Cell Ranger (Zheng et al., 2017)
    and Seurat (Macosko et al., 2015).

    Parameters
    ----------
    X : AnnData or array-like
        Data matrix storing unlogarithmized data.
    log : bool
        Use the logarithm of mean and variance.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=None : float
        Cutoffs for the gene expression, used if n_top_genes is None.
    n_top_genes : int or None (default: None)
        Number of highly-variable genes to keep.
    flavor : {'seurat', 'cell_ranger'}
        Choose method for computing normalized dispersion. Note that Seurat
        passes the cutoffs whereas Cell Ranger passes `n_top_genes`.
    plot : bool (default: False)
        Plot the result.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Notes
    -----
    If an AnnData is passed and copy == True, the following is returned, otherwise, adata is updated.
    adata : AnnData
        Filtered AnnData object.
    with the following fields to adata.var:
        means : np.ndarray of shape n_genes
            Means per gene.
        dispersions : np.ndarray of shape n_genes
            Dispersions per gene.
        dispersions_norm : np.ndarray of shape n_genes
            Dispersions per gene.
    If a data matrix is passed, the information is returned as tuple.
        gene_filter, means, dispersions, dispersion_norm
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        result = filter_genes_dispersion(adata.X, log=log,
                                         min_disp=min_disp, max_disp=max_disp,
                                         min_mean=min_mean, max_mean=max_mean,
                                         n_top_genes=n_top_genes,
                                         flavor=flavor, plot=plot)
        gene_filter, means, dispersions, dispersions_norm = result
        adata.var['means'] = means
        adata.var['dispersions'] = dispersions
        adata.var['dispersions_norm'] = dispersions_norm
        if plot:
            plot_filter_genes_dispersion(adata, gene_filter=gene_filter, log=not log)
        adata.filter_var(gene_filter)
        return adata if copy else None
    sett.m(0, '... filter highly varying genes by dispersion and mean')
    X = data  # proceed with data matrix
    if False:  # the following is less efficient and has no support for sparse matrices
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0, ddof=1)  # use R convention
        var = np.var(X, axis=0, ddof=1)
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=False).partial_fit(X)
        mean = scaler.mean_
        var = scaler.var_ * (X.shape[0]/(X.shape[0]-1))  # user R convention (unbiased estimator)
        dispersion = var / (mean + 1e-12)
        if log:  # consider logarithmized mean as in Seurat
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
        df['mean_bin'] = pd.cut(df['mean'],
                                np.r_[-np.inf, np.percentile(df['mean'],
                                                             np.arange(10, 105, 5)),
                                      np.inf])
        var_by_bin = pd.DataFrame()
        from statsmodels import robust
        import warnings  # this raises a warning we do not want to display
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            var_by_bin['bin_disp_median'] = df.groupby('mean_bin').apply(
                lambda group: np.median(group['dispersion']))
            var_by_bin['bin_disp_mad'] = df.groupby('mean_bin').apply(
                lambda group: robust.mad(group['dispersion']))
        df = df.merge(var_by_bin, left_on='mean_bin', right_index=True)
        df['dispersion_norm'] = np.abs(df['dispersion']
                                       - df['bin_disp_median']) \
                                       / df['bin_disp_mad']
    else:
        raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
    dispersion_norm = np.array(df['dispersion_norm'].values)
    if n_top_genes is not None:
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_genes-1]
        gene_filter = df['dispersion_norm'].values >= disp_cut_off
        sett.m(0, 'dispersion cutoff', disp_cut_off)
    else:
        sett.m(0, '    using `min_disp`, `max_disp`, `min_mean` and `max_mean`')
        sett.m(0, '--> set `n_top_genes` to simply select top-scoring genes instead')
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_filter = np.logical_and.reduce((mean > min_mean, mean < max_mean,
                                             dispersion_norm > min_disp,
                                             dispersion_norm < max_disp))
    return gene_filter, df['mean'].values, df['dispersion'].values, df['dispersion_norm'].values


def filter_genes_cv_deprecated(X, Ecutoff, cvFilter):
    """Filter genes by coefficient of variance and mean.

    See `filter_genes_dispersion`.
    """
    if issparse(X):
        raise ValueError('Not defined for sparse input. See `filter_genes_dispersion`.')
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.std(X, axis=0) / (np.mean(X, axis=0) + .0001) > cvFilter
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def filter_genes_fano_deprecated(X, Ecutoff, Vcutoff):
    """Filter genes by fano factor and mean.

    See `filter_genes_dispersion`.
    """
    if issparse(X):
        raise ValueError('Not defined for sparse input. See `filter_genes_dispersion`.')
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.var(X, axis=0) / (np.mean(X, axis=0) + .0001) > Vcutoff
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def log1p(data, copy=False):
    """Apply logarithm to count data "plus 1".

    Parameters
    ----------
    data : array-like or AnnData
        The data matrix.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        data.X = log1p(data.X)
        return adata if copy else None
    X = data  # proceed with data matrix
    if not issparse(X):
        return np.log1p(X)
    else:
        return X.log1p()


def pca(data, n_comps=10, zero_center=None, svd_solver='auto',
        random_state=None, recompute=True, mute=False, return_info=None, copy=False, dtype='float32'):
    """Embed data using PCA.

    Parameters
    ----------
    data : AnnData or array_like
        X : np.ndarray
            Data matrix of shape n_samples x n_variables.
    n_comps : int, optional (default: 10)
        Number of principal components to compute.
    zero_center : bool or None, optional (default: None)
        If True, compute standard PCA from Covariance matrix. Default for dense
        input. If False, omit zero-centering variables, which allows to handle
        sparse input efficiently. Defaults for sparse input.
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
            sett.m(0, '... not recomputing, using X_pca contained '
                   'in adata (set `recompute` to avoid this)')
            return adata
        else:
            sett.mt(0, 'compute PCA with n_comps =', n_comps, start=True)
            result = pca(adata.X, n_comps=n_comps, zero_center=zero_center,
                         svd_solver=svd_solver, random_state=random_state,
                         recompute=recompute, mute=mute, return_info=True)
            X_pca, components, pca_variance_ratio = result
            adata.smp['X_pca'] = X_pca  # this is multicolumn-sample annotation
            for icomp, comp in enumerate(components):
                adata.var['PC' + str(icomp+1)] = comp
            adata.add['pca_variance_ratio'] = pca_variance_ratio
            sett.mt(0, 'finished, added\n'
                    '    the data representation "X_pca" (adata.smp)\n'
                    '    the loadings "PC1", "PC2", ... (adata.var)\n'
                    '    and "pca_variance_ratio" (adata.add)')
        return adata if copy else None
    X = data  # proceed with data matrix
    from .. import settings as sett
    if X.shape[1] < n_comps:
        n_comps = X.shape[1] - 1
        sett.m(0, 'reducing number of computed PCs to',
               n_comps, 'as dim of data is only', X.shape[1])
    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
    from sklearn.decomposition import PCA, TruncatedSVD
    verbosity_level = np.inf if mute else 0
    if zero_center:
        if issparse(X):
            sett.m(0, 'pca: as `zero_center=True`, '
                   'sparse input is densified and may '
                   'lead to huge memory consumption')
            X = X.toarray()
        pca_ = PCA(n_components=n_comps, svd_solver=svd_solver)
    else:
        sett.m(verbosity_level, '... without zero-centering')
        pca_ = TruncatedSVD(n_components=n_comps)
    X_pca = pca_.fit_transform(X)
    if X_pca.dtype.descr != np.dtype(dtype).descr:
        X_pca = X_pca.astype(dtype)
    if False if return_info is None else return_info:
        return X_pca, pca_.components_, pca_.explained_variance_ratio_
    else:
        return X_pca


def normalize_per_cell(data, scale_factor=None, copy=False):
    """Normalize each cell.

    Normalize each cell by UMI count, so that every cell has the same total
    count.

    Similar functions are used, for example, by Cell Ranger (Zheng et al.,
    2017), Seurat (Macosko et al., 2015), Haghverdi et al. (2016) or Weinreb et
    al. (2016).

    Parameters
    ----------
    data : np.ndarray or AnnData
        Data matrix. Rows correspond to cells and columns to genes.
    scale_factor : float or None (default: None)
        If None, multiply by median.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Returns
    -------
    X_norm : np.ndarray
        Normalized version of the original expression matrix.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata.X = normalize_per_cell(adata.X, scale_factor)
        return adata if copy else None
    X = data  # proceed with the data matrix
    counts_per_cell = np.sum(X, axis=1)
    if issparse(X):
        counts_per_cell = counts_per_cell.A1
    if scale_factor is None:
        scale_factor = np.median(counts_per_cell)
    if not issparse(X):
        X = X / (counts_per_cell[:, np.newaxis] + 1e-6) * scale_factor
    else:
        Norm = sp.sparse.diags(scale_factor / (counts_per_cell + 1e-6))
        X = Norm.dot(X.tobsr()).tocsr()
    return X


def normalize_per_cell_weinreb16(X, max_fraction=1, mult_with_mean=False):
    """Normalize each cell.

    This is a legacy version. See `normalize_per_cell` instead.

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


def regress_out(adata, smp_keys, n_jobs=None, copy=False):
    """Regress out unwanted sources of variation.

    Yields a dense matrix.

    Parameters
    ----------
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    """
    sett.mt(0, 'regress out', smp_keys)
    if issparse(adata.X):
        sett.m(0, '... sparse input is densified and may '
               'lead to huge memory consumption')
    if not copy:
        sett.m(0, '... note that this is an inplace computation '
               'and will return None, set copy true if you want a copy')
    adata = adata.copy() if copy else adata
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    # the code here can still be much optimized
    regressors = np.array([adata.smp[key] for key in smp_keys]).T
    regressors = np.c_[np.ones(adata.X.shape[0]), regressors]
    len_chunk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_chunks = np.ceil(adata.X.shape[1] / len_chunk).astype(int)
    chunks = [np.arange(start, min(start + len_chunk, adata.X.shape[1]))
              for start in range(0, n_chunks * len_chunk, len_chunk)]
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    if sett._is_interactive:
        # from tqdm import tqdm_notebook as tqdm
        # does not work in Rodeo, should be solved sometime soon
        # tqdm = lambda x: x
        from tqdm import tqdm
        # sett.m(1, 'TODO: get nice waitbars also in interactive mode')
        sett.m(0, '... nicer progress bars as on command line come soon')
    else:
        from tqdm import tqdm
    for chunk in tqdm(chunks):
        result_lst = Parallel(n_jobs=n_jobs)(
            delayed(_regress_out)(
                col_index, adata.X, regressors) for col_index in chunk)
        for i_column, column in enumerate(chunk):
            adata.X[:, column] = result_lst[i_column]
    sett.mt(0, 'finished')
    return adata if copy else None


def scale(data, zero_center=None, max_value=None, copy=False):
    """Scale data to unit variance and zero mean (`if zero_center`).

    Parameters
    ----------
    zero_center : bool or None, optional (default: None)
        If False, omit zero-centering variables, which allows to handle sparse
        input efficiently. For dense input, defaults to True, for sparse input,
        defaults to False.
    max_value : float or None, optional (default: None)
        Clip to this value after scaling. Defaults to 10 if zero_center is True,
        otherwise np.inf (no cutoff).
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    """
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        adata.X = scale(adata.X, zero_center)
        return adata if copy else None
    X = data  # proceed with the data matrix
    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
    if zero_center and max_value:
        sett.m(0, 'scale_data: be very careful to use `max_value` without `zero_center`')
    if zero_center and max_value is None:
        max_value = 10
    if not zero_center and max_value is None:
        max_value = np.inf
    if zero_center and issparse(X):
        sett.m(0, 'scale_data: as `zero_center=True`, '
               'sparse input is densified and may '
               'lead to huge memory consumption')
        X = X.toarray()
    from sklearn.preprocessing import StandardScaler
    # the following doesn't use the unbiased estimator for variance
    # hence the result differs slightly from R's result
    scaler = StandardScaler(with_mean=zero_center).partial_fit(X)
    # user R convention (unbiased estimator)
    scaler.scale_ = scaler.scale_ * np.sqrt(X.shape[0]/(X.shape[0]-1))
    X_scaled = scaler.transform(X)
    X_scaled[X_scaled > max_value] = max_value  # np.clip not implementd for sparse matrices?
    return X_scaled


def subsample(data, subsample, seed=0, copy=False):
    """Subsample.

    Parameters
    ----------
    data : AnnData or array-like
        Annotated data matrix.
    subsample : int
        Subsample to a fraction of 1/subsample of the data.
    seed : int
        Random seed to change subsampling.
    copy : bool (default: False)
        If an AnnData is passed, determines whether a copy is returned.

    Notes
    -----
    Returns X, smp_indices if data is array-like, otherwise subsamples the passed
    AnnData (copy == False) or a copy of it (copy == True).
    """
    from .. import utils
    if not isinstance(data, AnnData):
        X = data
        return utils.subsample(X, subsample, seed)
    adata = data.copy() if copy else data
    _, smp_indices = utils.subsample(adata.X, subsample, seed)
    adata.filter_smp(smp_indices)
    for k in adata.smp_keys():
        # TODO: this should also be taken into account when slicing
        if k + '_masks' in adata.add:
            adata.add[k + '_masks'] = adata[k + '_masks'][:, smp_indices]
    return adata if copy else None


def zscore_deprecated(X):
    """Z-score standardize each variable/gene in X.

    Use `scale` instead.

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
# Plot result of preprocessing functions
# --------------------------------------------------------------------------------


def plot_filter_genes_dispersion(adata, gene_filter, log=True):
    """Plot dispersions vs. means for genes.

    Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() of Seurat.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    gene_filter : np.ndarray of shape n_top_genes of dtype bool
        Boolean index array.
    log : bool
        Plot on logarithmic axes.
    """
    from matplotlib import pyplot as pl
    for d in ['dispersions_norm', 'dispersions']:
        means, dispersions = adata.var['means'], adata.var[d]
        pl.figure(figsize=(4, 4))
        for label, color, mask in zip(['highly variable genes', 'other genes'],
                                      ['black', 'grey'],
                                      [gene_filter, ~gene_filter]):
            pl.scatter(means[mask], dispersions[mask],
                       label=label, c=color, s=1)
        if log:  # there's a bug in autoscale
            pl.xscale('log')
            pl.yscale('log')
            min_dispersion = np.min(dispersions)
            y_min = 0.95*min_dispersion if min_dispersion > 0 else 1e-3
            pl.xlim(0.95*np.min(means), 1.05*np.max(means))
            pl.ylim(y_min, 1.05*np.max(dispersions))
        pl.legend()
        pl.xlabel('mean expression of gene')
        pl.ylabel('dispersion of gene' + (' (normalized)' if 'norm' in d else ' (not normalized)'))
    from .. import plotting as plott
    plott.savefig_or_show('high_var_genes')


# --------------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------------


def _regress_out(col_index, responses, regressors):
    try:
        result = sm.GLM(responses[:, col_index],
                        regressors, family=sm.families.Gaussian()).fit()
        new_column = result.resid_response
    except PerfectSeparationError:  # this emulates R's behavior
        sett.m(0, 'warning: encountered PerfectSeparationError, setting to zero')
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

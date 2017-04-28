# Author: F. Alex Wolf (http://falexwolf.de)
"""
Simple Preprocessing Functions

Compositions of these functions are found in sc.preprocess.recipes.
"""

import sys
import numpy as np
import scipy as sp
from joblib import Parallel, delayed
from collections import OrderedDict
from scipy.sparse import issparse
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from ..classes.ann_data import AnnData
from .. import settings as sett


def filter_cells(data, min_counts=None, min_genes=None):
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

    Returns
    -------
    If data is a data matrix X
        cell_filter : np.ndarray
            Boolean index mask that does filtering. True means that the cell is
            kept. False means the cell is removed.
        number_per_cell: np.ndarray
            Either n_counts or n_genes per cell.
    otherwise:
        adata : AnnData
            The filtered adata object, with the count info stored in adata.smp.
    """
    if min_genes is not None and min_counts is not None:
        raise ValueError('Either provide min_counts or min_genes, but not both.')
    if min_genes is None and min_counts is None:
        raise ValueError('Provide one of min_counts or min_genes.')
    if isinstance(data, AnnData):
        adata = data
        cell_filter, number = filter_cells(adata.X, min_counts, min_genes)
        if min_genes is None:
            adata.smp['n_counts'] = number
        else:
            adata.smp['n_genes'] = number
        return adata[cell_filter]
    X = data  # proceed with processing the data matrix
    min_number = min_counts if min_genes is None else min_genes
    number_per_cell = np.sum(X if min_genes is None else X > 0, axis=1)
    if issparse(X):
       number_per_cell = number_per_cell.A1
    cell_filter = number_per_cell >= min_number
    sett.m(0, '... filtered out', np.sum(~cell_filter), 'outlier cells')
    return cell_filter, number_per_cell


def filter_genes(data, min_cells=None, min_counts=None):
    """
    Keep genes that have at least `min_counts` counts or are expressed in
    at least `min_cells` cells / samples.

    Filter out genes whose expression is not supported by sufficient statistical
    evidence.
    """
    if min_cells is not None and min_counts is not None:
        raise ValueError('Either specify min_counts or min_cells, but not both.')
    if min_cells is None and min_counts is None:
        raise ValueError('Provide one of min_counts or min_cells.')
    if isinstance(data, AnnData):
        adata = data
        gene_filter, number = filter_genes(adata.X, min_cells=min_cells,
                                           min_counts=min_counts)
        if min_cells is None:
            adata.var['n_counts'] = number
        else:
            adata.var['n_cells'] = number
        return adata[:, gene_filter]
    X = data  # proceed with processing the data matrix
    number_per_gene = np.sum(X if min_cells is None else X > 0, axis=0)
    min_number = min_counts if min_cells is None else min_cells
    if issparse(X):
       number_per_gene = number_per_gene.A1
    gene_filter = number_per_gene >= min_number
    sett.m(0, '... filtered out', np.sum(~gene_filter),
           'genes that are detected',
           'in less than ' + str(min_cells)  + ' cells' if min_counts is None
           else 'with less than ' + str(min_counts) + ' counts')
    return gene_filter, number_per_gene


def filter_genes_dispersion(data, log=True,
                            min_disp=0.5, max_disp=None,
                            min_mean=0.0125, max_mean=3,
                            n_top_genes=None,
                            norm_method='seurat',
                            plot=False):
    """
    Extract highly variable genes.

    Similar functions are used, for example, by Cell Ranger (Zheng et al., 2017)
    and Seurat (Macosko et al., 2015).

    Parameters
    ----------
    X : AnnData or array-like
        Data matrix storing unlogarithmized data.
    log : bool
        Use the logarithm of mean and variance.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=None : float
        Cutoffs for the gene expression.
    n_top_genes : int or None (default: None)
        Number of highly-variable genes to keep.
    plot : bool (default: False)
        Plot the result.

    Returns
    -------
    adata : AnnData
        Filtered AnnData object.
    Writes the following fields to adata.var:
        means : np.ndarray of shape n_genes
            Means per gene.
        dispersions : np.ndarray of shape n_genes
            Dispersions per gene.
    """
    if isinstance(data, AnnData):
        adata = data
        result = filter_genes_dispersion(adata.X, log=log,
                                         min_disp=min_disp, max_disp=max_disp,
                                         min_mean=min_mean, max_mean=max_mean,
                                         n_top_genes=n_top_genes,
                                         norm_method=norm_method, plot=plot)
        gene_filter, means, dispersions, dispersions_norm = result
        adata.var['means'] = means
        adata.var['dispersions'] = dispersions
        adata.var['dispersions_norm'] = dispersions_norm
        if plot:
            plot_filter_genes_dispersion(adata, gene_filter=gene_filter, log=not log)
        return adata[:, gene_filter]
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
    if norm_method == 'seurat':
        df['mean_bin'] = pd.cut(df['mean'], bins=20)
        disp_grouped = df.groupby('mean_bin')['dispersion']
        disp_mean_bin = disp_grouped.mean()
        disp_std_bin = disp_grouped.std(ddof=1)
        df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
                                 - disp_mean_bin[df['mean_bin']].values) \
                                 / disp_std_bin[df['mean_bin']].values
    elif norm_method == 'cell_ranger':
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
        raise ValueError('norm_method needs to be `seurat` or `cell_ranger`')
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


def filter_genes_cv(X, Ecutoff, cvFilter):
    """
    Filter genes by coefficient of variance and mean.
    """
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.std(X, axis=0) / (np.mean(X, axis=0) + .0001) > cvFilter
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def filter_genes_fano(X, Ecutoff, Vcutoff):
    """
    Filter genes by fano factor and mean.
    """
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.var(X, axis=0) / (np.mean(X, axis=0) + .0001) > Vcutoff
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def log1p(data):
    """
    Apply logarithm to count data "plus 1".
    """
    if isinstance(data, AnnData):
        return AnnData(log1p(data.X), data.smp, data.var, **data.add)
    X = data  # proceed with data matrix
    if not issparse(X):
        X = np.log1p(X)
    else:
        X = X.log1p()
    return X


def pca(data, n_comps=10, zero_center=None, svd_solver='auto',
        random_state=None, recompute=True, mute=False, return_info=None):
    """
    Embed data using PCA.

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

    Returns
    -------
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
        adata = data
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
                adata.var['PC' + str(icomp)] = comp
            adata.add['pca_variance_ratio'] = pca_variance_ratio
            sett.mt(0, 'finished, added\n'
                    '    "X_pca" to adata.smp, "PC1", "PC2", ... to adata.var\n'
                    '    and "pca_variance_ratio" to adata.add')
        return adata
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
    if False if return_info is None else return_info:
        return X_pca.astype(np.float32), pca_.components_, pca_.explained_variance_ratio_
    else:
        return X_pca.astype(np.float32)


def normalize_per_cell(data, scale_factor=None):
    """
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

    Returns
    -------
    X_norm : np.ndarray
        Normalized version of the original expression matrix.
    """
    if isinstance(data, AnnData):
        adata = data
        X = normalize_per_cell(adata.X, scale_factor)
        return AnnData(X, adata.smp, adata.var, **adata.add)
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
    """
    Normalize each cell by UMI count, so that every cell has the same total
    count.

    See normalize_per_cell(...).

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


def regress_out(adata, smp_keys, n_jobs=2):
    """
    Regress out unwanted sources of variation.

    Yields a dense matrix.
    """
    sett.mt(0, 'regress out', smp_keys)
    if issparse(adata.X):
        sett.m(0, '... sparse input is densified and may '
               'lead to huge memory consumption')
    # the code here can still be much optimized
    # ensuring a homogeneous data type seems to be necessary for GLM
    regressors = np.array([adata.smp[key].astype(float)
                           for key in smp_keys]).T
    regressors = np.c_[np.ones(adata.X.shape[0]), regressors]
    adata_corrected = adata.copy()
    len_junk = np.ceil(min(1000, adata.X.shape[1]) / n_jobs).astype(int)
    n_junks = np.ceil(adata.X.shape[1] / len_junk).astype(int)
    junks = [np.arange(start, min(start + len_junk, adata.X.shape[1]))
             for start in range(0, n_junks * len_junk, len_junk)]
    if issparse(adata_corrected.X):
        adata_corrected.X = adata_corrected.X.toarray()
    if sett._is_interactive:
        # from tqdm import tqdm_notebook as tqdm
        # does not work in Rodeo, should be solved sometime soon
        # tqdm = lambda x: x
        from tqdm import tqdm
        # sett.m(1, 'TODO: get nice waitbars also in interactive mode')
        sett.m(0, '... nicer progress bars as on command line come soon')
    else:
        from tqdm import tqdm
    for junk in tqdm(junks):
        result_lst = Parallel(n_jobs=n_jobs)(
                              delayed(_regress_out)(
                                      col_index, adata.X, regressors)
                                      for col_index in junk)
        for i_column, column in enumerate(junk):
            adata_corrected.X[:, column] = result_lst[i_column]
    sett.mt(0, 'finished')
    return adata_corrected


def scale(data, zero_center=None, max_value=None):
    """
    Scale data to unit variance and zero mean (`if zero_center`).

    Parameters
    ----------
    zero_center : bool or None, optional (default: None)
        If False, omit zero-centering variables, which allows to handle sparse
        input efficiently. For dense input, defaults to True, for sparse input,
        defaults to False.
    max_value : float or None, optional (default: None)
        Clip to this value after scaling. Defaults to 10 if zero_center is True,
        otherwise np.inf (no cutoff).
    """
    if isinstance(data, AnnData):
        adata = data
        X = scale(adata.X, zero_center)
        return AnnData(X, adata.smp, adata.var, **adata.add)
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


def subsample(data, subsample, seed=0):
    """
    Subsample.

    Parameters
    ----------
    data : AnnData or array-like
        Annotated data matrix.
    subsample : int
        Subsample to a fraction of 1/subsample of the data.
    seed : int
        Root to change subsampling.

    Returns
    -------
    adata : dict containing modified entries
        'row_names', 'expindices', 'explabels', 'expcolors'
    """
    from .. import utils
    if not isinstance(X, AnnData):
        X = data
        return utils.subsample(X, subsample, seed)
    adata = data
    _, smp_indices = utils.subsample(adata.X, subsample, seed)
    adata = adata[smp_indices, ]
    for k in adata.smp_keys():
        if k + '_masks' in adata.adata:  # TODO: this should also be taken into account when slicing
            adata[k + '_masks'] = adata[k + '_masks'][:, smp_indices]
    adata.add['subsampled'] = True
    return adata


def zscore(X):
    """
    Z-score standardize each variable/gene in X.

    Parameters
    ----------
    X : np.ndarray
        Data matrix. Rows correspond to cells and columns to genes.

    Returns
    -------
    XZ : np.ndarray
        Z-score standardized version of the data matrix.
    """
    means = np.tile(np.mean(X,axis=0)[None,:],(X.shape[0],1))
    stds = np.tile(np.std(X,axis=0)[None,:],(X.shape[0],1))
    return (X - means) / (stds + .0001)


#--------------------------------------------------------------------------------
# Plot result of preprocessing functions
#--------------------------------------------------------------------------------


def plot_filter_genes_dispersion(adata, gene_filter, log=True):
    """
    Plot dispersions vs. means for genes.

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
    for d in [ 'dispersions_norm', 'dispersions']:
        means, dispersions = adata.var['means'], adata.var[d]
        pl.figure()
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


#--------------------------------------------------------------------------------
# Helper Functions
#--------------------------------------------------------------------------------


def _regress_out(col_index, responses, regressors):
    try:
        result = sm.GLM(responses[:, col_index].todense().A1,
                        regressors, family=sm.families.Gaussian()).fit()
        new_column = result.resid_response
    except PerfectSeparationError:  # this emulates R's behavior
        sett.m(0, 'warning: encountered PerfectSeparationError, setting to zero')
        new_column = np.zeros(responses.shape[0])
    return new_column


def _regress_out_junk(junk, responses, regressors):
    junk_array = np.zeros((responses.shape[0], junk.size),
                           dtype=responses.dtype)
    for i, col_index in enumerate(junk):
        junk_array[:, i] = _regress_out(col_index,  responses, regressors)
    return junk_array


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

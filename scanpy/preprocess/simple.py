# Author: F. Alex Wolf (http://falexwolf.de)
"""
Preprocessing functions that take X as an argument
-> "simple" preprocessing functions
"""

import sys
import numpy as np
import scipy as sp
import pandas as pd
from joblib import Parallel, delayed
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from collections import OrderedDict
from scipy.sparse import issparse
from ..classes.ann_data import AnnData
from .. import settings as sett


def filter_cells(data, min_counts=None, min_genes=None):
    """
    Keep cells that have at least `min_counts` UMI counts or `min_genes`
    genes expressed.

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
    return cell_filter, number_per_cell


def filter_genes(data, min_counts=None, min_cells=None):
    """
    Keep genes that have at least `min_counts` counts or are expressed in
    `min_cells` cells.
    """
    if min_cells is not None and min_counts is not None:
        raise ValueError('Either specify min_counts or min_cells, but not both.')
    if min_cells is None and min_counts is None:
        raise ValueError('Provide one of min_counts or min_cells.')
    if isinstance(data, AnnData):
        adata = data
        gene_filter, number = filter_genes(adata.X, min_counts, min_cells)
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
    return gene_filter, number_per_gene


def filter_genes_zheng17(X, n_genes=1000, return_info=False):
    """
    Highly variable genes as in Zheng et al. (2017).

    Parameters
    ----------
    X : np.ndarrary, sp.sparse.matrix of shape n_samples x n_variables
        Data matrix.
    n_genes : int
        Number of highly-variable genes to keep.
    return_info : bool (default: False)
        Return means and dispersions for each gene.

    Returns
    -------
    gene_filter : np.ndarray of shape n_genes of dtype bool
        Boolean index array.
    If return_info, in addition:
        means : np.ndarray of shape n_genes
            Means per gene.
        dispersions : np.ndarray of shape n_genes
            Dispersions per gene.
    """
    from statsmodels import robust
    if False:  # the following is less efficient and has no support for sparse matrices
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0, ddof=1)  # use R convention
        var = np.var(X, axis=0, ddof=1)
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=False).partial_fit(X)
        mean = scaler.mean_
        var = scaler.var_ * (X.shape[0]/(X.shape[0]-1)) # user R convention (unbiased estimator)
        std = np.sqrt(var)
    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['mean'] = mean
    df['cv'] = std / (df['mean'] + 1e-6)
    df['var'] = var
    df['dispersion'] = df['var'] / (df['mean'] + 1e-6)
    df['mean_bin'] = pd.cut(df['mean'],
                            np.r_[-np.inf, np.percentile(df['mean'], np.arange(10, 105, 5)), np.inf])
    var_by_bin = pd.DataFrame()
    import warnings  # this raises a warning we do not want to display
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        var_by_bin['bin_disp_median'] = df.groupby('mean_bin').apply(
            lambda group: np.median(group['dispersion']))
        var_by_bin['bin_disp_mad'] = df.groupby('mean_bin').apply(
            lambda group: robust.mad(group['dispersion']))
    df = df.merge(var_by_bin, left_on='mean_bin', right_index=True)
    df['dispersion_norm'] = np.abs(df['dispersion'] - df['bin_disp_median']) / df['bin_disp_mad']
    dispersion_norm = np.array(df['dispersion_norm'].values)
    dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
    disp_cut_off = dispersion_norm[n_genes-1]
    sett.m(0, 'dispersion cutoff', disp_cut_off)
    # boolean index array for highly-variable genes
    gene_filter = df['dispersion_norm'].values >= disp_cut_off
    if return_info:
        return gene_filter, df['mean'].values, df['dispersion'].values
    else:
        return gene_filter


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


def pca(X, n_comps=50, zero_center=None, svd_solver='randomized',
        random_state=0, mute=False):
    """
    Return PCA representation of data.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    n_comps : int, optional (default: 50)
        Number of PCs to compute.
    zero_center : bool or None, optional (default: None)
        If True, compute standard PCA from Covariance matrix. If False, omit
        zero-centering variables, which allows to handle sparse input efficiently.
        For sparse intput, automatically defaults to False.
    svd_solver : str, optional (default: 'randomized')
        SVD solver to use. Either "arpac" for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or "randomized" for the randomized algorithm
        due to Halko (2009).

    Returns
    -------
    X_pca : np.ndarray
        Data projected on n_comps PCs.
    """
    if isinstance(X, AnnData):
        sys.exit('Use sc.pca(adata, ...) instead. This is function is only defined for a data matrix.')
    from .. import settings as sett
    if X.shape[1] < n_comps:
        n_comps = X.shape[1] - 1
        sett.m(0, 'reducing number of computed PCs to',
               n_comps, 'as dim of data is only', X.shape[1])
    try:
        from scipy.sparse import issparse
        zero_center = zero_center if zero_center is not None else False if issparse(X) else True
        from sklearn.decomposition import PCA, TruncatedSVD
        sett.mt(0 if not mute else 10, 'compute PCA with n_comps =', n_comps)
        if zero_center:
            if issparse(X):
                X = X.toarray()
            X_pca = PCA(n_components=n_comps, svd_solver=svd_solver).fit_transform(X)
        else:
            sett.m(0 if not mute else 10, '... without zero-centering')
            X_pca = TruncatedSVD(n_components=n_comps).fit_transform(X)
        sett.mt(0 if not mute else 10, 'finished')
    except ImportError:
        X_pca = _pca_fallback(X, n_comps=n_comps)
        sett.mt(0 if not mute else 10, 'preprocess: computed PCA using fallback code\n',
                '--> can be sped up by installing package scikit-learn\n',
                '    or by setting the option exact=False')
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


def cell_norm_weinreb16(X, max_fraction=1, mult_with_mean=False):
    """
    Normalize by UMI, so that every cell has the same total read count.

    Used, for example, by Haghverdi et al. (2016), Weinreb et al. (2016) or
    Zheng et al. (2016).

    Using squared Euclidian distance after this normalization will yield the
    same result as using cosine distance.
        sq_eucl_dist(Y1, Y2) = (Y1-Y2)^2
        = Y1^2 +Y2^2 - 2Y1*Y2 = 1 + 1 - 2 Y1*Y2 = 2*(1-(Y1*Y2))
        = 2*(1-(X1*X2)/(|X1|*|X2|)) = 2*cosine_dist(X1,X2)

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


def _regress_out(col_index, responses, regressors):
    try:
        result = sm.GLM(responses[:, col_index].todense().A1,
                        regressors, family=sm.families.Gaussian()).fit()
        new_column = result.resid_response
    except PerfectSeparationError:  # this emulates R's behavior
        new_column = np.zeros(responses.shape[0])
    return new_column


def _regress_out_junk(junk, responses, regressors):
    junk_array = np.zeros((responses.shape[0], junk.size),
                           dtype=responses.dtype)
    for i, col_index in enumerate(junk):
        junk_array[:, i] = _regress_out(col_index,  responses, regressors)
    return junk_array

def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

def regress_out(adata, smp_keys, n_jobs=2):
    """
    Regress out unwanted sources of variation.

    Yields a dense matrix.
    """
    sett.mt(0, 'regress out', smp_keys)
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
    adata_corrected.X = adata_corrected.X.toarray()
    if is_interactive():
        # from tqdm import tqdm_notebook as tqdm  # does not work in Rodeo, should be solved sometime soon
        sett.m('TODO: nice waitbars also in interactive mode')
    else:
        from tqdm import tqdm
    for junk in tqdm(junks):
        result_lst = Parallel(n_jobs=n_jobs)(
                              delayed(_regress_out)(
                                      col_index, adata.X, regressors)
                                      for col_index in junk)
        for i_column, column in enumerate(junk):
            adata_corrected.X[:, column] = result_lst[i_column]
    sett.mt(0, 'finished regress out')
    return adata_corrected


def subsample(adata, subsample, seed=0):
    """
    Subsample.

    Parameters
    ----------
    adata : AnnData
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
    _, smp_indices = utils.subsample(adata.X, subsample, seed)
    adata = adata[smp_indices, ]
    for key in ['X_pca']:
        if key in adata:
            adata[key] = adata[key][smp_indices]
    for k in adata.smp_keys():
        if k + '_masks' in adata:
            adata[k + '_masks'] = adata[k + '_masks'][:, smp_indices]
    adata['subsample'] = True
    return adata


def zscore(X):
    """
    Z-score standardize each column of X.

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


def plot_high_var_genes_zheng17(gene_filter, means, dispersions):
    """
    Plot dispersions vs. means for genes.

    Produces Supp. Fig. 5c of Zheng et al. (2017).

    Takes as parameters the result of high_var_genes_zheng17().

    Parameters
    ----------
    gene_filter : np.ndarray of shape n_genes of dtype bool
        Boolean index array.
    means : np.ndarray of shape n_genes
        Means per gene.
    dispersions : np.ndarray of shape n_genes
        Dispersion per gene.
    """
    from matplotlib import pyplot as pl
    from .. import plotting as plott
    for label, color, mask in zip(['highly variable genes', 'other genes'],
                                  ['black', 'grey'],
                                  [gene_filter, ~gene_filter]):
        pl.scatter(means[mask], dispersions[mask],
                   label=label, c=color, s=1)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlim(0.95*np.min(means), 1.05*np.max(means)) # there's a bug in autoscale
    pl.ylim(0.95*np.min(dispersions), 1.05*np.max(dispersions))
    pl.legend()
    pl.xlabel('means')
    pl.ylabel('dispersions')
    plott.savefig_or_show('high_var_genes')


#--------------------------------------------------------------------------------
# Helper Functions
#--------------------------------------------------------------------------------


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

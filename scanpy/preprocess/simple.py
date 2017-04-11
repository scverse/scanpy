# Author: F. Alex Wolf (http://falexwolf.de)
"""
Preprocessing functions that take X as an argument
-> "simple" preprocessing functions
"""

from ..classes.ann_data import AnnData
from .. import settings as sett
import numpy as np
import scipy as sp
import pandas as pd
from collections import OrderedDict


def filter_cells(X, min_reads):
    """
    Filter out cells with total UMI count < min_reads.

    Paramaters
    ----------
    X : np.ndarray
        Data matrix. Rows correspond to cells and columns to genes.
    min_reads : int
        Minimum number of reads required for a cell to survive filtering.

    Returns
    -------
    X : np.ndarray
        Filtered data matrix.
    cell_filter : np.ndarray
        Boolean mask that reports filtering. True means that the cell is
        kept. False means the cell is removed.
    """
    total_counts = np.sum(X, axis=1)
    cell_filter = total_counts >= min_reads
    return X[cell_filter], cell_filter


def high_var_genes_zheng17(X, n_genes=1000):
    """
    Highly variable genes as in Zheng et al. (2017).

    Parameters
    ----------
    X : np.ndarrary, sp.sparse.matrix of shape n_samples x n_variables
        Data matrix.
    n_genes : int
        Number of highly-variable genes to keep.

    Returns
    -------
    gene_filter : np.ndarray of shape n_genes of dtype bool
        Boolean index array.
    means : np.ndarray of shape n_genes
        Means per gene.
    dispersions : np.ndarray of shape n_genes
        Dispersions per gene.
    """
    from statsmodels import robust
    if False: # the following is less efficient and has no support for sparse matrices
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0, ddof=1) # use R convention
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
    return gene_filter, df['mean'].values, df['dispersion'].values


def gene_filter_cv(X, Ecutoff, cvFilter):
    """
    Filter genes by coefficient of variance and mean.
    """
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.std(X, axis=0) / (np.mean(X, axis=0) + .0001) > cvFilter
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def gene_filter_fano(X, Ecutoff, Vcutoff):
    """
    Filter genes by fano factor and mean.
    """
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.var(X, axis=0) / (np.mean(X, axis=0) + .0001) > Vcutoff
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter


def log1p(X):
    """
    Apply logarithm to count data "plus 1".
    """
    if isinstance(X, AnnData):
        return AnnData(log1p(X.X), X.smp, X.var, **X.add)
    if not sp.sparse.issparse(X):
        X = np.log1p(X)
    else:
        X = X.log1p()
    return X


log = log1p
""" Same as log1p. For backwards compatibility. """


def pca(X, n_comps=50, zero_center=None, svd_solver='randomized', random_state=0):
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
    from .. import settings as sett
    if X.shape[1] < n_comps:
        n_comps = X.shape[1]-1
        sett.m(0, 'reducing number of computed PCs to',
               n_comps, 'as dim of data is only', X.shape[1])
    try:
        from scipy.sparse import issparse
        zero_center = zero_center if zero_center is not None else False if issparse(X) else True
        from sklearn.decomposition import PCA, TruncatedSVD
        sett.mt(0, 'compute PCA with n_comps =', n_comps)
        if zero_center:
            if issparse(X):
                X = X.toarray()
            X_pca = PCA(n_components=n_comps, svd_solver=svd_solver).fit_transform(X)
        else:
            sett.m(0, '... without zero-centering')
            X_pca = TruncatedSVD(n_components=n_comps).fit_transform(X)
        sett.mt(0, 'finished')
        sett.m(1, '--> to speed this up, set option exact=False')
    except ImportError:
        X_pca = _pca_fallback(X, n_comps=n_comps)
        sett.mt(0, 'preprocess: computed PCA using fallback code\n',
                '--> can be sped up by installing package scikit-learn\n',
                '    or by setting the option exact=False')
    return X_pca


def smp_norm(X):
    """
    Normalize by UMI, so that every cell has the same total read count.

    Adapted from Zheng et al. (2016), see
    https://github.com/10XGenomics/single-cell-3prime-paper/.

    Similar functions are used, for example, by Haghverdi et al. (2016) and
    Weinreb et al. (2016).

    Parameters
    ----------
    X : np.ndarray
        Expression matrix. Rows correspond to cells and columns to genes.

    Returns
    -------
    X_norm : np.ndarray
        Normalized version of the original expression matrix.
    """
    counts_per_gene = np.sum(X, axis=0)
    counts_per_cell = np.sum(X, axis=1)
    if not sp.sparse.issparse(X):
        gene_filter = counts_per_gene >= 1
        X = X * np.median(counts_per_cell) / (counts_per_cell[:, np.newaxis] + 1e-6)
    else:
        gene_filter = np.flatnonzero(counts_per_gene.A1 >= 1)
        Norm = sp.sparse.diags(np.median(counts_per_cell.A1) / (counts_per_cell.A.ravel() + 1e-6))
        X = Norm.dot(X.tobsr()).tocsr()
    return X, gene_filter


def smp_norm_weinreb16(X, max_fraction=1, mult_with_mean=False):
    """
    Normalize by UMI, so that every cell has the same total read count.

    Used, for example, by Haghverdi et al. (2016), Weinreb et al. (2016) or
    Zheng et al. (2016).

    Using Euclidian distance after this normalization will yield the same result
    as using cosine distance.
        eucl_dist(Y1, Y2) = (Y1-Y2)^2
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

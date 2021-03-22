import itertools
from math import sqrt
from typing import Collection, Optional, Union

from scipy.sparse.linalg import lsqr
import statsmodels.formula.api as smf
import scipy
import pandas as pd
import numpy as np
from numpy import array
from scipy.sparse import coo_matrix
from anndata import AnnData

from .. import logging as logg

def _generate_sphere(lib_sizes):
    # Sorts cells by their library sizes, and generates an ordering vector
    # to arrange cells in a circle based on increasing/decreasing lib size.
    nlibs = len(lib_sizes)
    o = np.argsort(lib_sizes)
    even = np.arange(0, nlibs, 2)
    odd = np.arange(1, nlibs, 2)
    out = np.r_[o[even], np.flip(o[odd])]
    return np.r_[out, out]

def _clean_size_factors(
    size_factors,
    num_detected,
):

    idx_keep = size_factors > 0
    print(np.sum(idx_keep), len(size_factors))
    if all(idx_keep):
        return size_factors

    if len(size_factors) != len(num_detected):
        raise ValueError('Size factors and num_detected should have same length')

    X = size_factors[idx_keep]
    Y = num_detected[idx_keep]


    if len(X) < 3:
        raise ValueError('need at least three points for fit')
    
    #Solving for initial conditions
    lower = np.median(X) < X
    df_lm = pd.DataFrame(
        {
            'X': X[lower],
            'Y': Y[lower]
        }
    )
    lm_fit = smf.ols(formula='Y ~ 0 + X', data=df_lm).fit()
    A = lm_fit.params['X']

    below = np.argwhere(Y < A*X)
    top = below[np.argmax(Y[below])]
    B = (A / Y[top] - 1 / X[top])

    #using logs to enforce positivity
    df_nls = pd.DataFrame(
        {
            'logA': [A for i in range(len(X))],
            'logB': [B[0] for i in range(len(X))],
            'lX': np.log(X),
            'lY': np.log(Y),
        }
    )
    #Robustifying iterations with tricube weights
    #weights = np.ones(len(Y))
    #for i in np.arange(iterations):


    
    #nls_fit = smf.ols(
    #        formula='lY ~ logA + lX - lX + np.exp(logB)',
    #        data=df_nls
    #    ).fit()
    nls_fit = None
    if nls_fit is None:
        size_factors[~idx_keep] = np.min(size_factors[idx_keep])
        return size_factors

        #init = nls_fit.coef
        #resids = np.abs(nls_fit.resid)
        #bandwidth = np.median(resids) * NMADS
        #bandwidth[bandwidth < 1e-8] = 1e-8
        #weights = 
    failed = num_detected[~idx_keep]
    coefs = nls_fit.params
    new_sf = failed/(np.exp(coefs['logA']) - coefs['np.exp(logB)'] * failed)

    new_sf[new_sf < 0] = np.max(X)

    size_factors[~idx_keep] = new_sf
    return size_factors

def _pool_size_factors(
    X, 
    avg_cell, 
    sphere, 
    sizes
):
    # A function to estimate the pooled size factors and construct the linear equations.
    num_cells, num_genes = X.shape
    if num_cells == 0:
        raise ValueError("at least one cell required for normalization")
    if num_genes == 0:
        raise ValueError('insufficient features for median computation')
    # Checking the input sizes.
    num_sizes = len(sizes)
    pool_factor = np.empty([num_sizes, num_cells], np.float64)
    if num_sizes == 0: #what to do here?
        return pool_factor

    last_size = -1
    total_size = 0

    for s in sizes:
        if (s > num_cells):
            raise ValueError(
                "each element of sizes should be within [1, number of cells]"
            )
        if s < last_size:
            raise ValueError("sizes should be sorted")
        total_size += s
        last_size = s

    # Checking pseudo cell.
    if num_genes != len(avg_cell):
        raise ValueError(
            "length of pseudo-cell vector is not the same as the number of cells"
        )

    # Checking ordering.
    if len(sphere) < num_cells * 2 - 1: #do i need the -1 here?
        raise ValueError("ordering vector is too short for number of cells")
    for o in sphere:
        if o < 0 or o >= num_cells:
            raise ValueError("elements of ordering vector are out of range")

    if scipy.sparse.issparse(X):
        X = X.A

    row_num = np.zeros(total_size * num_cells)
    col_num = np.zeros(total_size * num_cells)
    pool_factor = np.zeros(num_sizes * num_cells)
    X_reordered = X[sphere]

    is_even = num_genes % 2 == 0
    halfway = num_genes // 2
    # Running through the sliding windows.
    idx = 0
    idx_sphere = 0
    for i, win in enumerate(range(num_cells)):
        X_reordered = np.vstack(
            (X_reordered[1:, :], X_reordered[:1, :])
        )
        combined = np.zeros(num_genes)
        for j, s in enumerate(sizes):
            row_num[idx:idx+s] = win
            col_num[idx:idx+s] = sphere[idx_sphere:idx_sphere+s]
            idx += s
            combined = X_reordered[:s].sum(axis=0)
            # Computing the ratio against the reference.
            ratios = np.array(combined / avg_cell).reshape(-1)
            # Computing the median (faster than partial sort).
            if is_even:
                med_tmp = ratios[halfway]
                pool_factor[i] = (med_tmp + ratios[halfway - 1]) / 2
            else:
                pool_factor[i] = ratios[halfway]
        idx_sphere += 1
    return row_num, col_num, pool_factor

def _create_lin_system(
    X,
    avg_cell,
    sphere,
    sizes
):
    """
        # Does the heavy lifting of computing pool-based size factors
        # and creating the linear system out of the equations for each pool.
    """
    row_idx = []
    col_idx = []
    output = []

    _row_idx, _col_idx, _output = _pool_size_factors(X, avg_cell, sphere, sizes)
    row_idx.append(_row_idx)
    col_idx.append(_col_idx)
    output.append(_output)
    num_cells = X.shape[0]
    row_idx.append(
        np.arange(num_cells) + num_cells * len(sizes)
    )

    col_idx.append(
        np.arange(num_cells)
    ) 

    LOWWEIGHT = 0.000001
    output.append(
        np.repeat(np.sqrt(LOWWEIGHT) / np.sum(avg_cell), num_cells)
    )

    eqn_values = np.repeat(
        [1, sqrt(LOWWEIGHT)], 
        [len(row_idx[0]), len(row_idx[1])]
    )

    row_idx = array([y for x in row_idx for y in np.array(x)])
    col_idx = array([y for x in col_idx for y in np.array(x)])
    output = array([y for x in output for y in np.array(x)])
    
    design = coo_matrix(
        (eqn_values, (row_idx, col_idx)), shape=(len(output), num_cells)
    )

    return design, output


def _rescale_clusters(
    avg_cell_clust,
    ref_clust,
    min_mean
):
    """
        Chooses a cluster as a reference and rescales all other clusters to the reference,
        based on the 'normalization factors' computed between pseudo-cells.
    """
    num_clust = len(avg_cell_clust)
    rescaling = np.empty(num_clust)
    for clust_idx in range(num_clust):
        clust = avg_cell_clust[clust_idx]
        clust_libsize = np.sum(clust)
        ref_clust_libsize = np.sum(ref_clust)
        to_use = (
            (
                (clust / clust_libsize) 
                + (clust / ref_clust_libsize)
            ) / 2
            * (
                (clust_libsize + ref_clust_libsize) 
                / 2
            )
        )
        to_use = (to_use >= min_mean)
        if not np.all(to_use):
            clust = clust[to_use]
            ref_clust_touse= ref_clust[to_use]
        rescaled_sf = np.nanmedian(np.divide(clust, ref_clust_touse))

        if np.isinf(rescaled_sf) or rescaled_sf < 0:
            rescaled_sf = np.sum(clust) / np.sum(ref_clust_touse)
        
        rescaling[clust_idx] = rescaled_sf
    return rescaling


def _limit_cluster_size(
    clusters,
    max_size
):
    """
        Limits the maximum cluster size to avoid problems with memory in Matrix::qr().
        Done by arbitrarily splitting large clusters so that they fall below max_size.
    """
    if max_size is None:
        return clusters
    new_clusters = np.zeros(len(clusters))
    counter = 0
    for i in np.unique(clusters):
        current = np.equal(i, clusters)
        num_cells_clust = np.sum(current)

        if num_cells_clust <= max_size:
            new_clusters[current] = counter
            counter = counter + 1
            continue

        mult = np.ceil(num_cells_clust / max_size)
        realloc = _np_rep(
            np.arange(1, mult + 1) - 1 + counter,
            length=num_cells_clust
        )
        new_clusters[current] = realloc
        counter = counter + mult

    return new_clusters.astype(int)


def _split(x, f):
    count = max(f) + 1
    return tuple(
        list(itertools.compress(x, (el==i for el in f)))
        for i in range(count)
    )

def _norm_counts(X, size_factors):
    if scipy.sparse.issparse(X):
        #r,c = X.nonzero()
        #X = scipy.sparse.csr_matrix((
        #    (1.0/size_factors)[r],
        #    (r,c)), 
        #    shape=(X.shape)
        #)
        X = X.A
    
    X /= size_factors[:, None]
    return X

def _np_rep(x, reps=1, each=False, length=0):
    """ implementation of functionality of rep() and rep_len() from R

    Attributes:
        x: numpy array, which will be flattened
        reps: int, number of times x should be repeated
        each: logical; should each element be repeated reps times before the next
        length: int, length desired; if >0, overrides reps argument
        
    Taken from https://stackoverflow.com/questions/46166933/python-numpy-equivalent-of-r-rep-and-rep-len-functions
    """
    if length > 0:
        reps = np.int(np.ceil(length / x.size))
    x = np.repeat(x, reps)
    if(not each):
        x = x.reshape(-1, reps).T.ravel() 
    if length > 0:
        x = x[0:length]
    return(x)

def _per_cluster_normalize(
    X, 
    sizes, 
    scaling=None,
    min_mean=1, 
    positive=False, 
):
    '''
    Computes the normalization factors _within_ each cluster,
    along with the reference pseudo-cell used for normalization.
    '''
    if scipy.sparse.issparse(X):
        X = scipy.sparse.csr_matrix(X)
    if scaling is None:
        scaling = X.sum(axis=1)
    if np.any(scaling == 0):
        raise ValueError('Cells should have non-zero lib. sizes or scale factors')
    scaling = np.squeeze(np.asarray(scaling))
    gene_expr = _norm_counts(
        X, 
        size_factors=scaling,
    )

    avg_cell = np.array(
        (np.mean(X, axis=0) * np.mean(scaling))
    ).reshape(-1)

    high_avg = (avg_cell >= min_mean)
    use_avg_cell = avg_cell
    if (not all(high_avg)):
        gene_expr = gene_expr[:, high_avg]
        use_avg_cell = avg_cell[high_avg]

    sphere = _generate_sphere(scaling)
    sizes = sizes[sizes <= gene_expr.shape[0]]

    A, b = _create_lin_system(
        gene_expr,
        use_avg_cell,
        sphere,
        sizes
    )

    final_norm_factors = lsqr(A, b)
    final_norm_factors = final_norm_factors[0]
    final_norm_factors = final_norm_factors * scaling

    if any(final_norm_factors < 0):
        if positive:
            num_detected = np.count_nonzero(X.A, axis=1)
            final_norm_factors = _clean_size_factors(
                final_norm_factors, 
                num_detected
            )
    return final_norm_factors, avg_cell

def normalize_scran(
    adata: AnnData,
    sizes: Collection[int] = np.arange(21, 102, 5),
    clusters: Union[str, Collection[Union[str, int]], None] = None,
    ref_clust: Optional[Union[str, int]] = None,
    max_clust_size: int = 3000,
    positive: bool = True,
    debug=True,
    scaling: Union[str, Collection[float], None] = None,
    min_mean: Union[int, float] = 0.1,
    subset_gene: Union[str, None] = None,
    key_added: Optional[str] = None,
    inplace: bool = True,
):
    """\
    One sentence to describe it.

    More details, ...

    Parameters
    ----------
    adata
        ...
    sizes
        A numeric vector of pool sizes, i.e., number of cells per pool
    clusters
        An optional factor specifying which cells belong to which cluster, 
        for deconvolution within clusters.
    ref_clust
        A cluster to be used as the reference cluster for inter-cluster 
        normalization.
    max_clust_size
        An integer scalar specifying the maximum number of cells in each cluster.
    positive
        A logical scalar indicating whether linear inverse models should be used
        to enforce positive estimates.
    scaling
        A numeric scalar containing scaling factors to adjust the counts prior 
        to computing size factors.
    min_mean
        A numeric scalar specifying the minimum (library size-adjusted) average
        count of genes to be used for normalization.
    subset_gene
        An integer, logical or character vector specifying the features to use.
    ...

    Returns
    -------
    """
    if isinstance(clusters, str):
        clusters = adata.obs[clusters].cat.codes.values
    if isinstance(scaling, str):
        scaling = adata.obs[scaling]
    if isinstance(subset_gene, str):
        subset_gene = adata.var[subset_gene]

    X = adata.X
    num_cells = X.shape[0]
    num_genes = X.shape[1]
    if clusters is None: 
        clusters = np.ones(num_cells)
    clusters = _limit_cluster_size(clusters, max_clust_size)

    assert num_cells == len(clusters), 'Not all cells assigned to clusters'

    indices = _split(np.arange(len(clusters)), clusters)

    assert len(indices) != 0 and all([len(x) for x in indices]) != 0, 'Cells unassigned'
    
    if scaling is not None and len(scaling) != num_cells:
        raise ValueError('Scaling should have same length as no of cells')
    if scaling is None:
        scaling = np.ones(num_cells)
    sizes = np.sort(sizes.astype(int))

    if len(sizes) != len(set(sizes)):
        raise ValueError('There are duplicated sizes')

    #guess min_mean, should implement?
    min_mean = np.maximum(min_mean, 1e-8)

    frag_X = []
    frag_scaling = []
    for i in np.arange(len(indices)):
        idx = indices[i]
        if len(indices) > 1 or idx != np.arange(idx):
            current = X[idx, :]
        else:
            current = X
        if (subset_gene is not None):
            current =  current[:, subset_gene]
        frag_X.append(current)
        frag_scaling.append(list(scaling[idx]))

    clust_norm_factors = []
    clust_avg_cell = []
    for i, _X in enumerate(frag_X):
        _clust_norm_factors, _clust_avg_cell =_per_cluster_normalize(
            X=_X,
            scaling=None,
            sizes=sizes,
            min_mean=min_mean,
            positive=positive 
        )
        clust_norm_factors.append(_clust_norm_factors)
        clust_avg_cell.append(_clust_avg_cell)
    
    non_zeroes = []
    if ref_clust is None:
        for i, _clust_avg_cell in enumerate(clust_avg_cell):
            num_counts = np.sum(_clust_avg_cell, axis=0)
            non_zeroes.append(
                np.sum(num_counts > 0)
            )
        ref_clust = clust_avg_cell[np.argmax(non_zeroes)]
    
    rescaling_factors = _rescale_clusters(
        clust_avg_cell, 
        ref_clust,
        min_mean,
        )
    clust_norm_factors_scaled = []
    for i in range(len(clust_norm_factors)):
        clust_norm_factors_scaled.append(
            clust_norm_factors[i] * rescaling_factors[i]
        )

    clust_norm_factors_scaled = [
        y for x in clust_norm_factors_scaled for y in np.array(x)
    ]

    final_size_factors = np.ones(num_cells)
    indices = [
        y for x in indices for y in np.array(x)
    ]
    final_size_factors[indices] = clust_norm_factors_scaled
    is_positive = ((final_size_factors > 0) & (~np.isnan(final_size_factors)))
    final_size_factors /= np.mean(final_size_factors[is_positive])

    size_factors = final_size_factors
    X_norm = X/size_factors[:, None]

    if debug:
        return X_norm, size_factors, clust_norm_factors_scaled, rescaling_factors

    if inplace:
        adata.X = X_norm
        if key_added is not None:
            adata.obs[key_added] = size_factors
    else:
        return dict(X=X_norm, size_factors=size_factors)

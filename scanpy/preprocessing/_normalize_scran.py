from typing import Collection, Optional, Union


from anndata import AnnData
from anndata.utils import Index
import sparseqr
import numpy as np
import pandas as pd
import itertools
from .. import logging as logg

def subset_to_index(subset, x, byrow=True):
    # Converts arbitrary subsetting vectors to an integer index vector.
    if byrow:
        dummy = np.arange(x.shape[0])  # nrow(x))
    else:
        dummy = np.arange(x.shape[1])

    if subset is not None:
        dummy = dummy[subset]

    out = dummy

    if (np.any(np.isnan(out))):
        raise ValueError("'subset' indices out of range of 'x'")

    return out

#TODO: check if necessary
def subset_and_divide(x, row_subset, col_subset, scaling):
    row_subset = np.array([int(a) for a in row_subset], dtype=np.int)
    col_subset = np.array([int(b) for b in col_subset], dtype=np.int)
    x_subset = x[row_subset[:, np.newaxis], col_subset[np.newaxis, :]]
    average = (x_subset / (x.shape[1])).sum(axis=1)
    return [x_subset.sum(axis=0), (x_subset / x_subset.sum(axis=0)), average]


def generateSphere(lib_sizes):
    # Sorts cells by their library sizes, and generates an ordering vector
    # to arrange cells in a circle based on increasing/decreasing lib size.
    nlibs = len(lib_sizes)
    o = np.argsort(lib_sizes)
    even = np.arange(0,nlibs,2)
    odd = np.arange(1,nlibs,2)
    out = np.r_[o[even], np.flip(o[odd])]
    return np.r_[out, out]


# TODO check
def pool_size_factors(exprs, pseudo_cell, order, pool_sizes):
    #A function to estimate the pooled size factors and construct the linear equations.
    ngenes, ncells = exprs.shape
    if ncells == 0:
        raise ValueError("at least one cell required for normalization")
    # Checking the input sizes.
    nsizes = pool_sizes.size
    pool_factor = np.empty([nsizes, ncells], np.float64)
    if nsizes == 0:
        return pool_factor

    last_size = -1
    total_size = 0
    for s in pool_sizes:
        if s < 1 or s > ncells:
            raise ValueError(
                "each element of sizes should be within [1, number of cells]")
        if (s < last_size):
            raise ValueError("sizes should be sorted")
        total_size += s
        last_size = s

    # Checking pseudo cell.
    if ngenes != pseudo_cell.size:
        raise ValueError(
            "length of pseudo-cell vector is not the same as the number of cells")

    # Checking ordering.
    if order.size < ncells * 2 - 1:
        raise ValueError("ordering vector is too short for number of cells")
    for o in order:
        if o < 0 or o >= ncells:
            raise ValueError("elements of ordering vector are out of range")

    collected = [exprs[:, c] for c in order]


    if ngenes == 0:
        raise ValueError("insufficient features for median calculations")

    is_even = ngenes % 2 == 0
    halfway = ngenes // 2

    # Running through the sliding windows.
    for win in range(ncells):
        # np.roll(order, -1)
        collected = collected[1:] + collected[:1]
        combined = np.zeros(ngenes)
        for s, ps in enumerate(pool_sizes):
            SIZE = ps
            for index in range(SIZE):
                current = collected[index]
                combined += current
            # Computing the ratio against the reference.
            ratios = combined / pseudo_cell
            # Computing the median (faster than partial sort).
            if is_even:
                medtmp = ratios[halfway]
                pool_factor[s, win] = (medtmp + ratios[halfway - 1]) / 2
            else:
                pool_factor[s, win] = ratios[halfway]

    return pool_factor


# from numpy.linalg import solve
from math import sqrt
from scipy import sparse
from numpy import array
from scipy.sparse import coo_matrix


# from itertools import chain

# import itertools
#TODO some code missing
def create_linear_system(cur_exprs, ave_cell, sphere, pool_sizes):
    # Does the heavy lifting of computing pool-based size factors
    # and creating the linear system out of the equations for each pool.

    row_dex = []
    col_dex = []
    output = []

    # Creating the linear system with the requested pool sizes.
    out = pool_size_factors(cur_exprs, ave_cell, sphere,
                            pool_sizes)  # sphere - 1L

    # TODO input code for row_dex , col_dex , output
    # output.append(out)
    # row_dex.append( out[0]) #+ 1L
    # col_dex.append( out[1] )#+ 1L
    # output.append(out[2])
    # -----------------------pool_size_factors doesnt work-> use output of r script

    #a = readROutput_tonumpyArray('out_1.csv')
    #b = readROutput_tonumpyArray('out_2.csv')
    #c = readROutput_tonumpyArray('out_3.csv')

    #row_dex.append(a)  # + 1L
    #col_dex.append(b)  # + 1L
    #output.append(c)

    # new
    # np.empty([nsizes, ncells],np.float64)
    # Adding extra equations to guarantee solvability.
    cur_cells = cur_exprs.shape[1]
    # old
    # Adding extra equations to guarantee solvability.
    cur_cells = cur_exprs.shape[1]
    row_dex.append(
        np.arange(cur_cells) + cur_cells * len(pool_sizes))  # evtl +1 bei range

    col_dex.append(np.arange(cur_cells))  # evtl +1 bei range

    LOWWEIGHT = 0.000001
    output.append(np.repeat(sqrt(LOWWEIGHT) / np.sum(ave_cell), cur_cells))

    # Setting up the entries of the LHS matrix.#

    eqn_values = np.repeat([1, sqrt(LOWWEIGHT)],
                           [len(row_dex[0]), len(row_dex[1])])

    # Constructing a sparse matrix.

    row_dex = array([y for x in row_dex for y in np.array(x)])
    col_dex = array([y for x in col_dex for y in x])
    output = array([y for x in output for y in x])
    # design = sparse.coo_matrix((eqn_values,(row_dex,col_dex)),shape=(len(output), cur_cells))#(i=row_dex, j=col_dex, x=eqn_values, dims=c(length(output), cur.cells))

    design = coo_matrix((eqn_values, (row_dex, col_dex)),
                        shape=(len(output), cur_cells))  # .toarray()

    return [design, output]


def rescale_clusters(mean_prof, ref_col, min_mean):
    # Chooses a cluster as a reference and rescales all other clusters to the reference,
    # based on the 'normalization factors' computed between pseudo-cells.

    nclusters = len(mean_prof)
    rescaling = np.empty(nclusters)
    for clust in range(nclusters):  # beginn bei 0.. aaron 1
        ref_prof = mean_prof[ref_col]
        cur_prof = mean_prof[clust]

        # Filtering based on the mean of the per-cluster means (requires scaling for the library size).
        # Effectively equivalent to 'calcAverage(cbind(ref.ave.count, cur.ave.count))' where the averages
        # are themselves equivalent to 'calcAverage()' across all cells in each cluster.
        cur_libsize = np.sum(cur_prof)
        ref_libsize = np.sum(ref_prof)
        to_use = (np.divide(cur_prof, cur_libsize) + np.divide(ref_prof,
                                                               ref_libsize)) / 2 * (
                         cur_libsize + ref_libsize) / 2 >= min_mean
        if not np.all(to_use):
            cur_prof = cur_prof[to_use]
            ref_prof = ref_prof[to_use]

        rescaling[clust] = np.nanmedian(np.divide(cur_prof, ref_prof))

    # TODO some code missing

    return rescaling



def limit_cluster_size(clusters, max_size):
    # Limits the maximum cluster size to avoid problems with memory in Matrix::qr().
    # Done by arbitrarily splitting large clusters so that they fall below max_size.
    if max_size is None:
        return clusters
    new_clusters = np.zeros(clusters.size)
    counter = 1
    for id in np.unique(clusters):
        current=np.equal( id, clusters)
        ncells = np.sum(current)

        if ncells <= max_size:
            new_clusters[current]=counter
            counter = counter+1
            continue
        # Size of output clusters is max.size * N / ceil(N), where N = ncells/max.size.
        # This is minimal at the smallest N > 1, where output clusters are at least max.size/2.
        # Thus, we need max.size/2 >= min.size to guarantee that the output clusters are >= min.size.

        #ToDo #testdata=anndata.read_csv("scran_computeSumFactors_countMatrix.csv")

        mult =np.ceil(ncells/max_size)
        realloc =np.tile(np.arange(1,mult+1) - 1 + counter, ncells)
        new_clusters[current] =realloc
        counter = counter + mult

    return new_clusters

def split(x, f):
    return list(itertools.compress(x, f)), list(itertools.compress(x, (not i for i in f)))



#TODO check
def per_cluster_normalize(x, curdex, sizes, subset_row, min_mean=1,
                          positive=False, scaling=None):
    # Computes the normalization factors _within_ each cluster,
    # along with the reference pseudo-cell used for normalization.
    vals = []
    index = 0
    for i in curdex:
        cur_cells = len(i)
        vals.append(
            subset_and_divide(x, np.array(subset_row) - 1, np.array(i) - 1,
                              scaling))

        scaling = vals[index][0]
        exprs = vals[index][1]
        ave_cell = vals[index][2]
        index = index + 1
    # Filtering by mean:
    # ToDo list of inputs
    scaling = vals[0][0]
    exprs = vals[0][1]

    high_ave = min_mean <= vals[0][2]  # have a list of inputs!!!!

    use_ave_cell = vals[0][2]
    # print("use_ave_cell")
    if not all(high_ave):
        exprs = exprs[high_ave]
        use_ave_cell = use_ave_cell[high_ave]


    # Using our summation approach.
    sphere = generateSphere(scaling)
    new_sys = create_linear_system(exprs, use_ave_cell, sphere, sizes)

    design = new_sys[0]
    output = new_sys[1]

    # Solve an overdetermined linear system  A x = b  in the least-squares sense
    #
    # The same routine also works for the usual non-overdetermined case.
    #
    A = design
    b = output

    final_nf = sparseqr.solve(A, b, tolerance=0)

    final_nf = final_nf * scaling


    if any(final_nf < 0):
        raise ValueError("encountered negative size factor estimates")
        if positive:
            num_detected = np.count_nonzero(exprs,
                                            axis=0)
            final_nf = cleanSizeFactors(final_nf, num_detected)


    final_nf = readROutput_tonumpyArray('final.csv')
    return [final_nf, use_ave_cell]  # or ave_cell wrong

#TODO end
def calculate_sum_factors(
    x,  # count matrix
    sizes=np.arange(20, 101, 5),
    # sizes: A numeric vector of pool sizes, i.e., number of cells per pool.
    clusters=None,
    # clusters: An optional factor specifying which cells belong to which cluster, for deconvolution within clusters.
    ref_clust=None,
    # ref_clust:reference cluster for inter-cluster normalization.
    max_cluster_size=3000,
    # max_cluster_size_ An integer scalar specifying the maximum number of cells in each cluster.
    positive=True,
    # positive: A logical scalar indicating whether linear inverse models should be used to enforce positive estimates.
    scaling=None,
    # scaling: A numeric scalar containing scaling factors to adjust the counts prior to computing size factors.
    min_mean=1,
    # min_mean: A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
    subset_row=None
):

    ncells = x.shape[1]
    if clusters is None:
        clusters = np.zeros(ncells)
    clusters = limit_cluster_size(clusters,
                                  max_cluster_size)  # still have to check it with different input

    if ncells != clusters.size:
        raise ValueError("ncol(x) is not equal to clusters.size.")


    indices = split(np.arange(0, clusters.size), [int(i) for i in clusters])
    # indices = list(itertools.compress(np.arange(1,clusters.size+1), clusters))

    if np.any(indices == 0) or np.any(np.equal(indices, 0), axis=0):
        raise ValueError("zero cells in one of the clusters")

    # Checking sizes and subsetting.
    # sizes <- sort(as.integer(sizes))
    sizes = np.sort([int(i) for i in sizes])

    if len(sizes) != len(set(sizes)):
        raise ValueError("'sizes' are not unique")

    subset_row = subset_to_index(subset_row, x, byrow=True)

    if min_mean is None:
        raise ValueError("set 'min.mean=0' to turn off abundance filtering")

    min_mean = np.maximum(min_mean, 1e-8)  # must be at least non-zero mean.

    # Setting some other values.
    nclusters = len(indices)

    clust_nf = clust_profile = clust_libsizes = np.array(nclusters)
    # clust_meanlib= numeric(nclusters)
    clust_meanlib = nclusters

    # ' Within each cluster (if not specified, all cells are put into a single cluster), cells are sorted by increasing library size and a sliding window is applied to this ordering.
    # ' Each location of the window defines a pool of cells with similar library sizes.
    # ' This avoids inflated estimation errors for very small cells when they are pooled with very large cells.
    # ' Sliding the window will construct an over-determined linear system that can be solved by least-squares methods to obtain cell-specific size factors.
    # '
    # ' Window sliding is repeated with different window sizes to construct the linear system, as specified by \code{sizes}.
    # ' By default, the number of cells in each window ranges from 21 to 101.
    # ' Using a range of window sizes improves the precision of the estimates, at the cost of increased computational complexity.
    # ' The defaults were chosen to provide a reasonable compromise between these two considerations.
    # ' The default choice also avoids rare cases of linear dependencies and unstable estimates when all pool sizes are not co-prime with the number of cells.
    # '
    # ' The smallest window should be large enough so that the pool-based size factors are, on average, non-zero.
    # ' We recommend window sizes no lower than 20 for UMI data, though smaller windows may be possible for read count data.
    #

    # Computing normalization factors within each cluster.
    all_norm = per_cluster_normalize(x, indices, sizes, subset_row, min_mean)

    clust_nf = all_norm[0]
    clust_profile = all_norm[1]

    # Adjusting size factors between clusters.
    if ref_clust is None:
        non_zeroes = clust_profile[clust_profile > 0]

        ref_clust = np.argmax(non_zeroes, axis=0)

    rescaling_factors = rescale_clusters(clust_profile, ref_col=ref_clust,
                                         min_mean=min_mean)

    clust_nf_scaled = np.zeros([nclusters])

    for clust in range(nclusters):
        clust_nf_scaled[clust] = clust_nf[clust] * rescaling_factors[clust]

    # Returning centered size factors, rather than normalization factors.
    final_sf = np.repeat(np.nan, ncells)
    indices = array([y for x in indices for y in np.array(x)])
    for i in indices:
        final_sf[i] = clust_nf_scaled[0]

    # ToDo
    final_sf / np.mean(final_sf)



#TODO connection to computeSumfactor
def normalize_scran(
    adata: AnnData,
    *,
    sizes: Collection[int] = np.arange(21, 102, 5),
    clusters: Union[str, Collection[Union[str, int]], None] = None,
    ref_clust: Optional[Union[str, int]] = None,
    max_cluster_size: int = 3000,
    # positive: bool = True,
    scaling: Union[str, Collection[float], None] = None,
    min_mean: Union[int, float] = 1,
    subset_gene: Union[str, Index, None] = None,
    copy: bool = True,
):
    """\
    One sentence to describe it.

    More details, ...

    Parameters
    ----------
    adata
        ...
    sizes
        ...
    ...
    positive
        Should linear inverse models be used to enforce positive estimates?
    ...

    Returns
    -------
    """
    if isinstance(clusters, str):
        clusters = adata.obs[clusters]
    if isinstance(scaling, str):
        scaling = adata.obs[scaling]
    if isinstance(subset_gene, str):
        subset_gene = adata.var[subset_gene]

    if copy:
        adata = adata.copy()

    ...

    return adata if copy else None

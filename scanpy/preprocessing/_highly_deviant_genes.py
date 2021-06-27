from typing import Optional
from .._compat import Literal
from anndata import AnnData
from .. import logging as logg

from statsmodels.stats.multitest import multipletests
import multiprocessing as mp
import scipy.sparse
import scipy.stats
import numpy as np

# this whole implementation is modelled after the original R code
# https://github.com/willtownes/scrna2019/blob/master/util/functions.R


def poisson_chunk(inds, Y, n):
    """
    Compute Poisson deviance for a chunk of the count matrix.

    Returns a vector with the per-gene deviance of the chunk.

    Input:
    ------
            inds - a pair of integer indices, defining the start and end gene indices for the chunk
            Y - count matrix, raw counts, cells as rows, genes as columns
            n - what appears as sz in the R implementation, for deviance computation
    """
    # what appears here as n appears as sz in the original R implementation
    # and mu got renamed to lambda
    # subset input data to specified genes for the chunk
    Yn = Y[:, inds[0] : inds[1]]
    # turn to dense if sparse. set up appropriate dimensionality for mu
    if scipy.sparse.issparse(Yn):
        Yn = np.asarray(Yn.todense())
    mu = np.sum(Yn, axis=0) / np.sum(n)
    mu = mu[None, :]
    # in R, they skip all the nan values from the thing defined here as holder
    # seeing how it gets multiplied by the values from Y, instead of trying to find some weird skip
    # just set all the zero values at this step to 1. this way they become 0 after logging
    # and successfully go away in the multiplication. gives same results as matching R command
    holder = Yn / (n * mu)
    holder[holder == 0] = 1
    deviance = 2 * np.sum(Yn * np.log(holder), axis=0) - 2 * np.sum(Yn - n * mu, axis=0)
    return deviance


def binomial_chunk(inds, Y, n):
    """
    Compute binomial deviance for a chunk of the count matrix.

    Returns a vector with the per-gene deviance of the chunk.

    Input:
    ------
            inds - a pair of integer indices, defining the start and end gene indices for the chunk
            Y - count matrix, raw counts, cells as rows, genes as columns
            n - n in the R implementation, for deviance computation
    """
    # subset input data to specified genes for the chunk
    Yn = Y[:, inds[0] : inds[1]]
    # turn to dense if sparse. set up appropriate dimensionality for p
    if scipy.sparse.issparse(Yn):
        Yn = np.asarray(Yn.todense())
    p = np.sum(Yn, axis=0) / np.sum(n)
    p = p[None, :]
    # in R, they skip all the nan values from the thing defined here as holder
    # seeing how it gets multiplied by the values from Y, instead of trying to find some weird skip
    # just set all the zero values at this step to 1. this way they become 0 after logging
    # and successfully go away in the multiplication. gives same results as matching R command
    holder = Yn / (n * p)
    holder[holder == 0] = 1
    term1 = np.sum(Yn * np.log(holder), axis=0)
    nx = n - Yn
    # same thing, second holder
    holder = nx / (n * (1 - p))
    holder[holder == 0] = 1
    term2 = np.sum(nx * np.log(holder), axis=0)
    return 2 * (term1 + term2)


# initiating global variables for the chunk processes to use
def _pool_init(_Y, _n, _chunk_func):
    global Y, n, chunk_func
    Y = _Y
    n = _n
    chunk_func = _chunk_func


# a function that calls the requisite chunk computation based on inds and global variables
def pool_func(inds):
    return chunk_func(inds, Y, n)


def compute_deviance(Y, flavor='poisson', chunksize=1e8, n_jobs=None):
    """
    A function to compute the Poisson/binomial deviance of a count matrix.

    Returns a vector of one deviance per gene

    Input:
    ------
            Y : ``scipy.sparse`` array  or ``numpy.array``
                    Observations as rows, features as columns, raw count (integer) data
            all other input as in highly_deviant_genes()
    """
    # if no info on n_jobs, just go for all the cores
    if n_jobs is None:
        n_jobs = mp.cpu_count()
    # computation prep - obtain n(bin)/sz(poiss), store in variable called n regardless
    # and point chunk_func to appropriate function that computes the deviance of chunks
    if flavor == 'poisson':
        chunk_func = poisson_chunk
        # this n was sz in the R version
        n = np.log(np.sum(Y, axis=1))
        n = np.exp(n - np.mean(n))
    elif flavor == 'binomial':
        chunk_func = binomial_chunk
        n = np.sum(Y, axis=1)
    else:
        raise ValueError('incorrect flavor value')
    # ensure that n has the correct dimensionality, needs to be treated differently if sparse counts
    if scipy.sparse.issparse(Y):
        n = np.asarray(n)
    else:
        n = n[:, None]
    # how many genes per job? we're aiming for chunksize total values in memory at a time
    # so divide that by n_jobs and the number of cells to get genes per job
    chunkcount = np.ceil(chunksize / (Y.shape[0] * n_jobs))
    # set up chunkcount-wide index intervals for count matrix subsetting within the jobs
    inds = []
    for ind in np.arange(0, Y.shape[1], chunkcount):
        inds.append([np.int(ind), np.int(ind + chunkcount)])
    # set up the parallel computation pool with the count matrix and n
    # chunk_func is set to the corresponding Poisson/binomial computation function
    p = mp.Pool(n_jobs, _pool_init, (Y, n, chunk_func))
    deviance = p.map(pool_func, inds)
    # collapse output and return
    return np.hstack(deviance)


def highly_deviant_genes(
    adata: AnnData,
    layer: Optional[str] = None,
    flavor: Literal['poisson', 'binomial'] = 'poisson',
    alpha: Optional[float] = 0.05,
    n_top_genes: Optional[int] = None,
    chunksize: Optional[int] = 1e8,
    n_jobs: Optional[int] = None,
):
    """
    Computation of highly deviant genes [Townes19]_

    Expects count data with unexpressed genes removed.

    Identifies genes statistically significant in a Chi-square test with the number of
    cells informing the number of degrees of freedom (FDR correction via Benjamini-Hochberg).
    If ``n_top_genes`` is provided, instead reports the user-specified number of genes with
    the highest deviance.

    Parameters
    ----------
    adata
            Object to compute highly deviant genes for. Requires raw count space with
            unexpressed genes removed.
    layer
            If provided, use ``adata.layers[layer]`` for expression values instead of ``adata.X``.
    flavor
            Choose between "poisson" and "binomial" deviance
    alpha
            Significance threshold for Benjamini-Hochberg corrected Chi-square testing of
            the deviances, with the number of cells informing the number of degrees of
            freedom. Ignored if ``n_top_genes`` is provided.
    n_top_genes
            If provided, skip the statistical testing and report these many highest
            deviances as the highly deviant genes instead.
    chunksize
            How many values from the count matrix to process at once (across all parallel
            processes). Higher values mean fewer total matrix operations necessary, but
            constrained by memory availability.
    n_jobs
            How many parallel jobs to use for the parallelisation. If ``None``, use all
            available cores.

    Returns
    -------
    Updates ``.var`` with the following fields.

    highly_deviant : bool
            A Boolean indicator of highly deviant genes.
    deviance : float
            The calculated deviance for each of the genes.
    deviance_pval : float
            The Chi-square p-value of the gene being highly deviant; only computed if
            ``n_top_genes`` is ``None``.
    deviance_qval : float
            The p-value after multiple testing correction via Benjamini-Hochberg, used to
            identify the highly deviant genes; only computed if ``n_top_genes`` is ``None``.
    """
    # logging
    start = logg.info('extracting highly deviant genes')
    # run the deviance computation. store in object
    if layer is not None:
        deviance = compute_deviance(
            adata.layers[layer], flavor=flavor, chunksize=chunksize, n_jobs=n_jobs
        )
    else:
        deviance = compute_deviance(
            adata.X, flavor=flavor, chunksize=chunksize, n_jobs=n_jobs
        )
    adata.var['deviance'] = deviance
    adata.uns['hdg'] = {'flavor': flavor}
    # are we just picking the top N genes?
    if n_top_genes is not None:
        # create a vector with highly variable information to insert
        highly_deviant = np.array([False] * adata.shape[1])
        # get the indices of the elements sorted in descending order, then get the top ones
        highly_deviant[np.argsort(deviance)[::-1][:n_top_genes]] = True
        adata.var['highly_deviant'] = highly_deviant
        adata.uns['hdg']['n_top_genes'] = n_top_genes
        # logging
        logg.info('    finished', time=start)
        logg.hint(
            'added\n'
            '    \'highly_deviant\', boolean vector (adata.var)\n'
            '    \'deviance\', float vector (adata.var)\n'
        )
    else:
        # statistical testing - run a chi-square test, then get a benjamini-hochberg correction
        pvals = 1 - scipy.stats.chi2.cdf(deviance, adata.shape[1] - 1)
        qvals = multipletests(pvals, method="fdr_bh")[1]
        adata.var['deviance_pval'] = pvals
        adata.var['deviance_qval'] = qvals
        adata.var['highly_deviant'] = qvals <= alpha
        adata.uns['hdg']['alpha'] = alpha
        # logging
        logg.info('    finished', time=start)
        logg.hint(
            'added\n'
            '    \'highly_deviant\', boolean vector (adata.var)\n'
            '    \'deviance\', float vector (adata.var)\n'
            '    \'deviance_pval\', float vector (adata.var)\n'
            '    \'deviance_qval\', float vector (adata.var)\n'
        )

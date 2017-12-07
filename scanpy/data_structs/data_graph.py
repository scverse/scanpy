# Author: Alex Wolf (http://falexwolf.de)
"""Data Graph

Represent a data matrix as a weighted graph of nearest neighbor relations
(edges) among data points (nodes).
"""

import numpy as np
import scipy as sp
import scipy.spatial
import scipy.sparse
from scipy.sparse import issparse
from joblib import Parallel, delayed
from .. import settings as sett
from .. import logging as logg
from .. import utils


N_DCS = 15  # default number of diffusion components
N_PCS = 50  # default number of PCs

def add_or_update_graph_in_adata(
        adata,
        n_neighbors=None,
        n_pcs=None,
        n_dcs=None,
        knn=None,
        recompute_pca=False,
        recompute_distances=False,
        recompute_graph=False,
        n_jobs=None):
    graph = DataGraph(adata,
                      k=n_neighbors,
                      n_pcs=n_pcs,
                      n_dcs=n_dcs,
                      recompute_pca=recompute_pca,
                      recompute_distances=recompute_distances,
                      recompute_graph=recompute_graph,
                      n_jobs=n_jobs)
    if graph.fresh_compute:
        graph.update_diffmap()
        adata.uns['data_graph_distance_local'] = graph.Dsq
        adata.uns['data_graph_norm_weights'] = graph.Ktilde
        adata.obsm['X_diffmap'] = graph.rbasis[:, 1:]
        adata.obs['X_diffmap0'] = graph.rbasis[:, 0]
        adata.uns['diffmap_evals'] = graph.evals[1:]
    return graph


def no_recompute_of_graph_necessary(
        adata,
        recompute_pca=False,
        recompute_distances=False,
        recompute_graph=False,
        n_neighbors=None,
        knn=None,
        n_dcs=None):
    conditions_base = [
        not recompute_pca,
        not recompute_distances,
        not recompute_graph,
        # make sure X_diffmap is there
        'X_diffmap' in adata.obsm_keys(),
        # make sure data_graph is there
        'data_graph_norm_weights' in adata.uns
        ]
    if not all(conditions_base):
        return False
    else:
        conditions = [
            # make sure enough DCs are there
            (adata.obsm['X_diffmap'].shape[1] >= n_dcs-1
                 if n_dcs is not None else True),
            # make sure that it's sparse
            (issparse(adata.uns['data_graph_norm_weights']) == knn
                 if knn is not None else True),
            # make sure n_neighbors matches
            (n_neighbors == adata.uns[
                'data_graph_distance_local'][0].nonzero()[0].size + 1
                if n_neighbors is not None else True)]
        return all(conditions)


def get_neighbors(X, Y, k):
    Dsq = utils.comp_sqeuclidean_distance_using_matrix_mult(X, Y)
    chunk_range = np.arange(Dsq.shape[0])[:, None]
    indices_chunk = np.argpartition(Dsq, k-1, axis=1)[:, :k]
    indices_chunk = indices_chunk[chunk_range,
                                  np.argsort(Dsq[chunk_range, indices_chunk])]
    indices_chunk = indices_chunk[:, 1:]  # exclude first data point (point itself)
    distances_chunk = Dsq[chunk_range, indices_chunk]
    return indices_chunk, distances_chunk


def get_distance_matrix_and_neighbors(X, k, sparse=True, n_jobs=1):
    """Compute distance matrix in squared Euclidian norm.
    """
    if not sparse:
        if False: Dsq = utils.comp_distance(X, metric='sqeuclidean')
        else: Dsq = utils.comp_sqeuclidean_distance_using_matrix_mult(X, X)
        sample_range = np.arange(Dsq.shape[0])[:, None]
        indices = np.argpartition(Dsq, k-1, axis=1)[:, :k]
        indices = indices[sample_range, np.argsort(Dsq[sample_range, indices])]
        indices = indices[:, 1:]  # exclude first data point (point itself)
        distances = Dsq[sample_range, indices]
    elif X.shape[0] > 1e5:
        # sklearn is slower, but for large sample numbers more stable
        from sklearn.neighbors import NearestNeighbors
        sklearn_neighbors = NearestNeighbors(n_neighbors=k-1, n_jobs=n_jobs)
        sklearn_neighbors.fit(X)
        distances, indices = sklearn_neighbors.kneighbors()
        distances = distances.astype('float32')**2
    else:
        # assume we can fit at max 20000 data points into memory
        len_chunk = np.ceil(min(20000, X.shape[0]) / n_jobs).astype(int)
        n_chunks = np.ceil(X.shape[0] / len_chunk).astype(int)
        chunks = [np.arange(start, min(start + len_chunk, X.shape[0]))
                 for start in range(0, n_chunks * len_chunk, len_chunk)]
        indices = np.zeros((X.shape[0], k-1), dtype=int)
        distances = np.zeros((X.shape[0], k-1), dtype=np.float32)
        if n_jobs > 1:
            # set backend threading, said to be meaningful for computations
            # with compiled code. more important: avoids hangs
            # when using Parallel below, threading is much slower than
            # multiprocessing
            result_lst = Parallel(n_jobs=n_jobs, backend='threading')(
                delayed(get_neighbors)(X[chunk], X, k) for chunk in chunks)
        else:
            logg.info('--> can be sped up by setting `n_jobs` > 1')
        for i_chunk, chunk in enumerate(chunks):
            if n_jobs > 1:
                indices_chunk, distances_chunk = result_lst[i_chunk]
            else:
                indices_chunk, distances_chunk = get_neighbors(X[chunk], X, k)
            indices[chunk] = indices_chunk
            distances[chunk] = distances_chunk
    if sparse:
        Dsq = get_sparse_distance_matrix(indices, distances, X.shape[0], k)
    return Dsq, indices, distances


def get_sparse_distance_matrix(indices, distances, n_samples, k):
    n_neighbors = k - 1
    n_nonzero = n_samples * n_neighbors
    indptr = np.arange(0, n_nonzero + 1, n_neighbors)
    Dsq = sp.sparse.csr_matrix((distances.ravel(),
                                indices.ravel(),
                                indptr),
                                shape=(n_samples, n_samples))
    return Dsq


def get_indices_distances_from_sparse_matrix(Dsq, k):
    indices = np.zeros((Dsq.shape[0], k), dtype=int)
    distances = np.zeros((Dsq.shape[0], k), dtype=Dsq.dtype)
    for i in range(indices.shape[0]):
        neighbors = Dsq[i].nonzero()
        # account for the fact that the first neighbor is the data point
        # itself...
        indices[i][0] = i
        indices[i][1:] = neighbors[1]
        distances[i][0] = 0
        distances[i][1:] = Dsq[neighbors]
    return indices, distances


class OnFlySymMatrix():
    """Emulate a matrix where elements are calculated on the fly.
    """
    def __init__(self, get_row, shape, DC_start=0, DC_end=-1, rows=None, restrict_array=None):
        self.get_row = get_row
        self.shape = shape
        self.DC_start = DC_start
        self.DC_end = DC_end
        self.rows = {} if rows is None else rows
        self.restrict_array = restrict_array  # restrict the array to a subset

    def __getitem__(self, index):
        if isinstance(index, int) or isinstance(index, np.integer):
            if self.restrict_array is None:
                glob_index = index
            else:
                # map the index back to the global index
                glob_index = self.restrict_array[index]
            if glob_index not in self.rows:
                self.rows[glob_index] = self.get_row(glob_index)
            row = self.rows[glob_index]
            if self.restrict_array is None:
                return row
            else:
                return row[self.restrict_array]
        else:
            if self.restrict_array is None:
                glob_index_0, glob_index_1 = index
            else:
                glob_index_0 = self.restrict_array[index[0]]
                glob_index_1 = self.restrict_array[index[1]]
            if glob_index_0 not in self.rows:
                self.rows[glob_index_0] = self.get_row(glob_index_0)
            return self.rows[glob_index_0][glob_index_1]

    def restrict(self, index_array):
        """Generate a view restricted to a subset of indices.
        """
        new_shape = index_array.shape[0], index_array.shape[0]
        return OnFlySymMatrix(self.get_row, new_shape, DC_start=self.DC_start,
                              DC_end=self.DC_end,
                              rows=self.rows, restrict_array=index_array)


class DataGraph():
    """Data represented as graph of nearest neighbors.

    Represent a data matrix as a weighted graph of nearest neighbor relations
    (edges) among data points (nodes).

    This is much more efficient than, for example, the scikit-learn implementation.
    """

    def __init__(self,
                 adata,
                 k=None,
                 knn=True,
                 n_jobs=None,
                 n_pcs=None,
                 n_dcs=None,
                 recompute_pca=False,
                 recompute_distances=False,
                 recompute_graph=False,
                 flavor='haghverdi16'):
        self.sym = True  # we do not allow asymetric cases
        self.flavor = flavor  # this is to experiment around
        self.n_pcs = n_pcs if n_pcs is not None else N_PCS
        self.init_iroot_and_X(adata, recompute_pca, n_pcs)
        # use the graph in adata
        if no_recompute_of_graph_necessary(
                adata,
                recompute_pca=recompute_pca,
                recompute_distances=recompute_distances,
                recompute_graph=recompute_graph,
                n_neighbors=k,
                knn=knn,
                n_dcs=n_dcs):
            self.fresh_compute = False
            self.knn = issparse(adata.uns['data_graph_norm_weights'])
            self.Ktilde = adata.uns['data_graph_norm_weights']
            self.Dsq = adata.uns['data_graph_distance_local']
            if self.knn:
                self.k = adata.uns['data_graph_distance_local'][0].nonzero()[0].size + 1
            else:
                self.k = None  # currently do not store this, is unknown
            # for output of spectrum
            if n_dcs is None: n_dcs = adata.obsm['X_diffmap'].shape[1] + 1
            self.X_diffmap = adata.obsm['X_diffmap'][:, :n_dcs-1]
            self.evals = np.r_[1, adata.uns['diffmap_evals'][:n_dcs-1]]
            self.n_dcs = len(self.evals)
            self.rbasis = np.c_[adata.obs['X_diffmap0'].values[:, None],
                                adata.obsm['X_diffmap'][:, :n_dcs-1]]
            self.lbasis = self.rbasis
            self.Dchosen = OnFlySymMatrix(self.get_Ddiff_row,
                                          shape=(self.X.shape[0], self.X.shape[0]))
            np.set_printoptions(precision=10)
            logg.info('    using stored data graph with n_neighbors = {} and '
                      'spectrum\n    {}'
                      .format(self.k,
                              str(self.evals).replace('\n', '\n    ')))
        # recompute the graph
        else:
            self.fresh_compute = True
            self.n_dcs = n_dcs if n_dcs is not None else N_DCS
            self.k = k if k is not None else 30
            if self.k > adata.n_obs:
                self.k = 1 + int(0.5*adata.n_obs)
            logg.info('    computing data graph with n_neighbors = {} '
                      .format(self.k))
            self.evals = None
            self.rbasis = None
            self.lbasis = None
            self.X_diffmap = None
            self.Dsq = None
            self.knn = knn
            self.n_jobs = sett.n_jobs if n_jobs is None else n_jobs
            self.Dchosen = None
            if False:  # TODO
                # in case we already computed distance relations
                if not recompute_distances and 'data_graph_distance_local' in adata.uns:
                    n_neighbors = adata.uns['data_graph_distance_local'][0].nonzero()[0].size + 1
                    if (knn and issparse(adata.uns['data_graph_distance_local'])
                        and n_neighbors == self.k):
                        logg.info('    using stored distances with `n_neighbors={}`'
                                  .format(self.k))
                        self.Dsq = adata.uns['data_graph_distance_local']

    def init_iroot_directly(self, adata):
        self.iroot = None
        if 'iroot' in adata.uns:
            if adata.uns['iroot'] >= adata.n_obs:
                logg.warn('Root cell index {} does not exist for {} samples. '
                          'Is ignored.'
                          .format(adata.uns['iroot'], adata.n_obs))
            else:
                self.iroot = adata.uns['iroot']

    def init_iroot_and_X(self, adata, recompute_pca, n_pcs):
        self.X = adata.X  # might be overwritten with X_pca in the next line
        # retrieve xroot
        xroot = None
        if 'xroot' in adata.uns: xroot = adata.uns['xroot']
        elif 'xroot' in adata.var: xroot = adata.var['xroot']
        # set iroot directly
        self.init_iroot_directly(adata)
        # see whether we can set self.iroot using the full data matrix
        if xroot is not None and xroot.size == self.X.shape[1]:
            self.set_root(xroot)
        # use the full data matrix X, nothing to be done
        if self.n_pcs == 0 or self.X.shape[1] <= self.n_pcs:
            logg.info('    using data matrix X directly for building graph (no PCA)')
        # use X_pca
        else:
            # use a precomputed X_pca
            if (not recompute_pca
                and 'X_pca' in adata.obsm_keys()
                and adata.obsm['X_pca'].shape[1] >= self.n_pcs):
                logg.info('    using \'X_pca\' with n_pcs = {} for building graph'
                          .format(self.n_pcs))
            # compute X_pca
            else:
                logg.info('    compute \'X_pca\' with n_pcs = {} for building graph'
                          .format(self.n_pcs))
                from ..preprocessing import pca
                pca(adata, n_comps=self.n_pcs)
            # set the data matrix
            self.X = adata.obsm['X_pca'][:, :n_pcs]
            # see whether we can find xroot using X_pca
            if xroot is not None and xroot.size == adata.obs['X_pca'].shape[1]:
                self.set_root(xroot[:n_pcs])

    def update_diffmap(self, n_comps=None):
        """Diffusion Map as of Coifman et al. (2005) and Haghverdi et al. (2016).
        """
        if n_comps is not None:
            self.n_dcs = n_comps
            logg.info('    updating number of DCs to', self.n_dcs)
        if self.evals is None or self.evals.size < self.n_dcs:
            logg.info('    computing spectral decomposition ("diffmap") with',
                      self.n_dcs, 'components', r=True)
            self.compute_transition_matrix()
            self.embed(n_evals=self.n_dcs)
            return True
        return False

    def compute_Ddiff_all(self, n_evals=10):
        raise RuntimeError('deprecated function')
        self.embed(n_evals=n_evals)
        self.compute_M_matrix()
        self.compute_Ddiff_matrix()

    def compute_C_all(self, n_evals=10):
        self.compute_L_matrix()
        self.embed(self.L, n_evals=n_evals, sort='increase')
        evalsL = self.evals
        self.compute_Lp_matrix()
        self.compute_C_matrix()

    def spec_layout(self):
        self.compute_transition_matrix()
        self.compute_L_matrix()
        self.embed(self.L, sort='increase')
        # write results to dictionary
        ddmap = {}
        # skip the first eigenvalue/eigenvector
        ddmap['Y'] = self.rbasis[:, 1:]
        ddmap['evals'] = self.evals[1:]
        return ddmap

    def compute_distance_matrix(self):
        logg.m('computing distance matrix with n_neighbors = {}'
               .format(self.k), v=4)
        Dsq, indices, distances_sq = get_distance_matrix_and_neighbors(
            X=self.X,
            k=self.k,
            sparse=self.knn,
            n_jobs=self.n_jobs)
        self.Dsq = Dsq
        return Dsq, indices, distances_sq

    def compute_transition_matrix(self, alpha=1, recompute_distance=False):
        """Compute transition matrix.

        Parameters
        ----------
        alpha : float
            The density rescaling parameter of Coifman and Lafon (2006). Should
            in all practical applications equal 1: Then only the geometry of the
            data matters, not the sampled density.
        neglect_selfloops : bool
            Discard selfloops.

        References
        ----------
        Haghverdi et al. (2016), Coifman and Lafon (2006), Coifman et al. (2005).
        """
        if self.Dsq is None or recompute_distance:
            Dsq, indices, distances_sq = self.compute_distance_matrix()
        else:
            Dsq = self.Dsq
            indices, distances_sq = get_indices_distances_from_sparse_matrix(Dsq, self.k)
        # choose sigma, the heuristic here often makes not much
        # of a difference, but is used to reproduce the figures
        # of Haghverdi et al. (2016)
        if self.knn:
            # as the distances are not sorted except for last element
            # take median
            sigmas_sq = np.median(distances_sq, axis=1)
        else:
            # the last item is already in its sorted position as
            # argpartition puts the (k-1)th element - starting to count from
            # zero - in its sorted position
            sigmas_sq = distances_sq[:, -1]/4
        sigmas = np.sqrt(sigmas_sq)
        logg.m('determined n_neighbors =',
               self.k, 'nearest neighbors of each point', t=True, v=4)

        if self.flavor == 'unweighted':
            if not self.knn:
                raise ValueError('`flavor="unweighted"` only with `knn=True`.')
            self.Ktilde = self.Dsq.sign()
            return

        # compute the symmetric weight matrix
        if not sp.sparse.issparse(self.Dsq):
            Num = 2 * np.multiply.outer(sigmas, sigmas)
            Den = np.add.outer(sigmas_sq, sigmas_sq)
            W = np.sqrt(Num/Den) * np.exp(-Dsq/Den)
            # make the weight matrix sparse
            if not self.knn:
                self.Mask = W > 1e-14
                W[self.Mask == False] = 0
            else:
                # restrict number of neighbors to ~k
                # build a symmetric mask
                Mask = np.zeros(Dsq.shape, dtype=bool)
                for i, row in enumerate(indices):
                    Mask[i, row] = True
                    for j in row:
                        if i not in set(indices[j]):
                            W[j, i] = W[i, j]
                            Mask[j, i] = True
                # set all entries that are not nearest neighbors to zero
                W[Mask == False] = 0
                self.Mask = Mask
        else:
            W = Dsq
            for i in range(len(Dsq.indptr[:-1])):
                row = Dsq.indices[Dsq.indptr[i]:Dsq.indptr[i+1]]
                num = 2 * sigmas[i] * sigmas[row]
                den = sigmas_sq[i] + sigmas_sq[row]
                W.data[Dsq.indptr[i]:Dsq.indptr[i+1]] = np.sqrt(num/den) * np.exp(-Dsq.data[Dsq.indptr[i]: Dsq.indptr[i+1]] / den)
            W = W.tolil()
            for i, row in enumerate(indices):
                for j in row:
                    if i not in set(indices[j]):
                        W[j, i] = W[i, j]
            if False:
                W.setdiag(1)  # set diagonal to one
                logg.m('    note that now, we set the diagonal of the weight matrix to one!')
            W = W.tocsr()
        logg.m('computed W (weight matrix) with "knn" =', self.knn, t=True, v=4)

        if False:
            pl.matshow(W)
            pl.title('$ W$')
            pl.colorbar()

        # density normalization
        # as discussed in Coifman et al. (2005)
        # ensure that kernel matrix is independent of sampling density
        if alpha == 0:
            # nothing happens here, simply use the isotropic similarity matrix
            self.K = W
        else:
            # q[i] is an estimate for the sampling density at point x_i
            # it's also the degree of the underlying graph
            if not sp.sparse.issparse(W):
                q = np.sum(W, axis=0)
                # raise to power alpha
                if alpha != 1:
                    q = q**alpha
                Den = np.outer(q, q)
                self.K = W / Den
            else:
                q = np.array(W.sum(axis=0)).flatten()
                self.K = W
                for i in range(len(W.indptr[:-1])):
                    row = W.indices[W.indptr[i]: W.indptr[i+1]]
                    num = q[i] * q[row]
                    W.data[W.indptr[i]: W.indptr[i+1]] = W.data[W.indptr[i]: W.indptr[i+1]] / num
        logg.m('computed K (anisotropic kernel)', t=True, v=4)

        if not sp.sparse.issparse(self.K):
            # now compute the row normalization to build the transition matrix T
            # and the adjoint Ktilde: both have the same spectrum
            self.z = np.sum(self.K, axis=0)
            # the following is the transition matrix
            self.T = self.K / self.z[:, np.newaxis]
            # now we need the square root of the density
            self.sqrtz = np.array(np.sqrt(self.z))
            # now compute the density-normalized Kernel
            # it's still symmetric
            szszT = np.outer(self.sqrtz, self.sqrtz)
            self.Ktilde = self.K / szszT
        else:
            self.z = np.array(self.K.sum(axis=0)).flatten()
            # now we need the square root of the density
            self.sqrtz = np.array(np.sqrt(self.z))
            # now compute the density-normalized Kernel
            # it's still symmetric
            self.Ktilde = self.K
            for i in range(len(self.K.indptr[:-1])):
                row = self.K.indices[self.K.indptr[i]: self.K.indptr[i+1]]
                num = self.sqrtz[i] * self.sqrtz[row]
                self.Ktilde.data[self.K.indptr[i]: self.K.indptr[i+1]] = self.K.data[self.K.indptr[i]: self.K.indptr[i+1]] / num
        logg.m('computed Ktilde (normalized anistropic kernel)', v=4)

    def compute_L_matrix(self):
        """Graph Laplacian for K.
        """
        self.L = np.diag(self.z) - self.K
        logg.info('compute graph Laplacian')

    def embed(self, matrix=None, n_evals=15, sym=None, sort='decrease'):
        """Compute eigen decomposition of matrix.

        Parameters
        ----------
        matrix : np.ndarray
            Matrix to diagonalize.
        n_evals : int
            Number of eigenvalues/vectors to be computed, set n_evals = 0 if
            you need all eigenvectors.
        sym : bool
            Instead of computing the eigendecomposition of the assymetric
            transition matrix, computed the eigendecomposition of the symmetric
            Ktilde matrix.

        Writes attributes
        -----------------
        evals : np.ndarray
            Eigenvalues of transition matrix
        lbasis : np.ndarray
            Matrix of left eigenvectors (stored in columns).
        rbasis : np.ndarray
             Matrix of right eigenvectors (stored in columns).
             self.rbasis is projection of data matrix on right eigenvectors,
             that is, the projection on the diffusion components.
             these are simply the components of the right eigenvectors
             and can directly be used for plotting.
        """
        np.set_printoptions(precision=10)
        if sym is None: sym = self.sym
        self.rbasisBool = True
        if matrix is None: matrix = self.Ktilde
        # compute the spectrum
        if n_evals == 0:
            evals, evecs = sp.linalg.eigh(matrix)
        else:
            n_evals = min(matrix.shape[0]-1, n_evals)
            # ncv = max(2 * n_evals + 1, int(np.sqrt(matrix.shape[0])))
            ncv = None
            which = 'LM' if sort == 'decrease' else 'SM'
            # it pays off to increase the stability with a bit more precision
            matrix = matrix.astype(np.float64)
            evals, evecs = sp.sparse.linalg.eigsh(matrix, k=n_evals,
                                                  which=which, ncv=ncv)
            evals, evecs = evals.astype(np.float32), evecs.astype(np.float32)
        if sort == 'decrease':
            evals = evals[::-1]
            evecs = evecs[:, ::-1]
        if logg.verbosity_greater_or_equal_than(4):
            logg.m('computed eigenvalues', t=True, v=4)
        else:
            logg.info('    eigenvalues of transition matrix')
        logg.info('   ', str(evals).replace('\n', '\n    '))
        # assign attributes
        self.evals = evals
        count_ones = sum([1 for v in self.evals if v == 1])
        if count_ones > len(self.evals)/2:
            logg.warn('Transition matrix has many irreducible blocks!')
        if sym:
            self.rbasis = self.lbasis = evecs
        else:
            # The eigenvectors of T are stored in self.rbasis and self.lbasis
            # and are simple trafos of the eigenvectors of Ktilde.
            # rbasis and lbasis are right and left eigenvectors, respectively
            self.rbasis = np.array(evecs / self.sqrtz[:, np.newaxis])
            self.lbasis = np.array(evecs * self.sqrtz[:, np.newaxis])
            # normalize in L2 norm
            # note that, in contrast to that, a probability distribution
            # on the graph is normalized in L1 norm
            # therefore, the eigenbasis in this normalization does not correspond
            # to a probability distribution on the graph
            if False:
                self.rbasis /= np.linalg.norm(self.rbasis, axis=0, ord=2)
                self.lbasis /= np.linalg.norm(self.lbasis, axis=0, ord=2)
        # init on-the-fly computed distance "matrix"
        self.Dchosen = OnFlySymMatrix(self.get_Ddiff_row,
                                      shape=self.Dsq.shape)

    def _get_M_row_chunk(self, i_range):
        from ..cython import utils_cy
        M_chunk = np.zeros((len(i_range), self.X.shape[0]), dtype=np.float32)
        for i_cnt, i in enumerate(i_range):
            if False:  # not much slower, but slower
                M_chunk[i_cnt] = self.get_M_row(j)
            else:
                M_chunk[i_cnt] = utils_cy.get_M_row(i, self.evals, self.rbasis, self.lbasis)
        return M_chunk

    def compute_M_matrix(self):
        """The M matrix is the matrix that results from summing over all powers of
        T in the subspace without the first eigenspace.

        See Haghverdi et al. (2016).
        """
        if self.n_jobs >= 4:  # if we have enough cores, skip this step
            return            # TODO: make sure that this is really the best strategy
        logg.m('    try computing "M" matrix using up to 90% of `sett.max_memory`')
        if True:  # Python version
            self.M = sum([self.evals[l]/(1-self.evals[l])
                          * np.outer(self.rbasis[:, l], self.lbasis[:, l])
                          for l in range(1, self.evals.size)])
            self.M += np.outer(self.rbasis[:, 0], self.lbasis[:, 0])
        else:  # Cython version
            used_memory, _ = logg.get_memory_usage()
            memory_for_M = self.X.shape[0]**2 * 23 / 8 / 1e9  # in GB
            logg.m('    max memory =', sett.max_memory,
                   ' / used memory = {:.1f}'.format(used_memory),
                   ' / memory_for_M = {:.1f}'.format(memory_for_M))
            if used_memory + memory_for_M < 0.9 * sett.max_memory:
                logg.m(0, '    allocate memory and compute M matrix')
                len_chunk = np.ceil(self.X.shape[0] / self.n_jobs).astype(int)
                n_chunks = np.ceil(self.X.shape[0] / len_chunk).astype(int)
                chunks = [np.arange(start, min(start + len_chunk, self.X.shape[0]))
                         for start in range(0, n_chunks * len_chunk, len_chunk)]
                # parallel computing does not seem to help
                if False:  # self.n_jobs > 1:
                    # here backend threading is not necessary, and seems to slow
                    # down everything considerably
                    result_lst = Parallel(n_jobs=self.n_jobs, backend='threading')(
                        delayed(self._get_M_row_chunk)(chunk)
                        for chunk in chunks)
                self.M = np.zeros((self.X.shape[0], self.X.shape[0]),
                                  dtype=np.float32)
                for i_chunk, chunk in enumerate(chunks):
                    if False:  # self.n_jobs > 1:
                        M_chunk = result_lst[i_chunk]
                    else:
                        M_chunk = self._get_M_row_chunk(chunk)
                    self.M[chunk] = M_chunk
                # the following did not work
                # filename = sett.writedir + 'tmp.npy'
                # np.save(filename, self.M)
                # self.M = filename
                sett.mt(0, 'finished computation of M')
            else:
                logg.m('not enough memory to compute M, using "on-the-fly" computation')

    def compute_Ddiff_matrix(self):
        """Returns the distance matrix in the Diffusion Pseudotime metric.

        See Haghverdi et al. (2016).

        Notes
        -----
        - Is based on M matrix.
        - self.Ddiff[self.iroot,:] stores diffusion pseudotime as a vector.
        """
        if self.M.shape[0] > 1000:
            logg.m('--> high number of dimensions for computing DPT distance matrix\n'
                   '    computing PCA with 50 components')
            from ..preprocessing import pca
            self.M = pca(self.M, n_comps=50, mute=True)
        self.Ddiff = sp.spatial.distance.squareform(sp.spatial.distance.pdist(self.M))
        logg.m('computed Ddiff distance matrix', t=True)
        self.Dchosen = self.Ddiff

    def _get_Ddiff_row_chunk(self, m_i, j_range):
        from ..cython import utils_cy
        M = self.M  # caching with a file on disk did not work
        d_i = np.zeros(len(j_range))
        for j_cnt, j in enumerate(j_range):
            if False:  # not much slower, but slower
                m_j = self.get_M_row(j)
                d_i[j_cnt] = sp.spatial.distance.cdist(m_i[None, :], m_j[None, :])
            else:
                if M is None:
                    m_j = utils_cy.get_M_row(j, self.evals, self.rbasis, self.lbasis)
                else:
                    m_j = M[j]
                d_i[j_cnt] = utils_cy.c_dist(m_i, m_j)
        return d_i

    def get_Ddiff_row(self, i):
        if not self.sym:
            raise ValueError('Not bug-free implemented! '
                             'Computation needs to be adjusted if sym=False.')
        row = sum([(self.evals[l]/(1-self.evals[l])
                     * (self.rbasis[i, l] - self.lbasis[:, l]))**2
                   # account for float32 precision
                    for l in range(0, self.evals.size) if self.evals[l] < 0.999999])
        row += sum([(self.rbasis[i, l] - self.lbasis[:, l])**2
                    for l in range(0, self.evals.size) if self.evals[l] >= 0.999999])
        return np.sqrt(row)

    def get_Ddiff_row_deprecated(self, i):
        from ..cython import utils_cy
        if self.M is None:
            m_i = utils_cy.get_M_row(i, self.evals, self.rbasis, self.lbasis)
        else:
            m_i = self.M[i]
        len_chunk = np.ceil(self.X.shape[0] / self.n_jobs).astype(int)
        n_chunks = np.ceil(self.X.shape[0] / len_chunk).astype(int)
        chunks = [np.arange(start, min(start + len_chunk, self.X.shape[0]))
                  for start in range(0, n_chunks * len_chunk, len_chunk)]
        if self.n_jobs >= 4:  # problems with high memory calculations, we skip computing M above
            # here backend threading is not necessary, and seems to slow
            # down everything considerably
            result_lst = Parallel(n_jobs=self.n_jobs)(
                delayed(self._get_Ddiff_row_chunk)(m_i, chunk)
                for chunk in chunks)
        d_i = np.zeros(self.X.shape[0])
        for i_chunk, chunk in enumerate(chunks):
            if self.n_jobs >= 4: d_i_chunk = result_lst[i_chunk]
            else: d_i_chunk = self._get_Ddiff_row_chunk(m_i, chunk)
            d_i[chunk] = d_i_chunk
        return d_i

    def compute_Lp_matrix(self):
        """See Fouss et al. (2006) and von Luxburg et al. (2007).

        See Proposition 6 in von Luxburg (2007) and the inline equations
        right in the text above.
        """
        self.Lp = sum([1/self.evals[i]
                      * np.outer(self.rbasis[:, i], self.lbasis[:, i])
                      for i in range(1, self.evals.size)])
        sett.mt(0, 'computed pseudoinverse of Laplacian')

    def compute_C_matrix(self):
        """See Fouss et al. (2006) and von Luxburg et al. (2007).

        This is the commute-time matrix. It's a squared-euclidian distance
        matrix in \mathbb{R}^n.
        """
        self.C = np.repeat(np.diag(self.Lp)[:, np.newaxis],
                           self.Lp.shape[0], axis=1)
        self.C += np.repeat(np.diag(self.Lp)[np.newaxis, :],
                            self.Lp.shape[0], axis=0)
        self.C -= 2*self.Lp
        # the following is much slower
        # self.C = np.zeros(self.Lp.shape)
        # for i in range(self.Lp.shape[0]):
        #     for j in range(self.Lp.shape[1]):
        #         self.C[i, j] = self.Lp[i, i] + self.Lp[j, j] - 2*self.Lp[i, j]
        volG = np.sum(self.z)
        self.C *= volG
        sett.mt(0, 'computed commute distance matrix')
        self.Dchosen = self.C

    def compute_MFP_matrix(self):
        """See Fouss et al. (2006).

        This is the mean-first passage time matrix. It's not a distance.

        Mfp[i, k] := m(k|i) in the notation of Fouss et al. (2006). This
        corresponds to the standard notation for transition matrices (left index
        initial state, right index final state, i.e. a right-stochastic
        matrix, with each row summing to one).
        """
        self.MFP = np.zeros(self.Lp.shape)
        for i in range(self.Lp.shape[0]):
            for k in range(self.Lp.shape[1]):
                for j in range(self.Lp.shape[1]):
                    self.MFP[i, k] += (self.Lp[i, j] - self.Lp[i, k]
                                       - self.Lp[k, j] + self.Lp[k, k]) * self.z[j]
        sett.mt(0, 'computed mean first passage time matrix')
        self.Dchosen = self.MFP

    def set_pseudotime(self):
        """Return pseudotime with respect to root point.
        """
        self.pseudotime = self.Dchosen[self.iroot].copy()
        self.pseudotime /= np.max(self.pseudotime)

    def set_root(self, xroot):
        """Determine the index of the root cell.

        Given an expression vector, find the observation index that is closest
        to this vector.

        Parameters
        ----------
        xroot : np.ndarray
            Vector that marks the root cell, the vector storing the initial
            condition, only relevant for computing pseudotime.
        """
        if self.X.shape[1] != xroot.size:
            raise ValueError('The root vector you provided does not have the '
                             'correct dimension. Make sure you provide the dimension-'
                             'reduced version, if you provided X_pca.')
        # this is the squared distance
        dsqroot = 1e10
        iroot = 0
        for i in range(self.X.shape[0]):
            diff = self.X[i, :] - xroot
            dsq = diff.dot(diff)
            if dsq < dsqroot:
                dsqroot = dsq
                iroot = i
                if np.sqrt(dsqroot) < 1e-10: break
        logg.m('setting root index to', iroot, v=4)
        if self.iroot is not None and iroot != self.iroot:
            logg.warn('Changing index of iroot from {} to {}.'.format(self.iroot, iroot))
        self.iroot = iroot
        return self.iroot

    def _test_embed(self):
        """
        Checks and tests for embed.
        """
        # pl.semilogy(w,'x',label=r'$ \widetilde K$')
        # pl.show()
        if sett.verbosity > 2:
            # output of spectrum of K for comparison
            w, v = np.linalg.eigh(self.K)
            sett.mi('spectrum of K (kernel)')
        if sett.verbosity > 3:
            # direct computation of spectrum of T
            w, vl, vr = sp.linalg.eig(self.T, left=True)
            sett.mi('spectrum of transition matrix (should be same as of Ktilde)')

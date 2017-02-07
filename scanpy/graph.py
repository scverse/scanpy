# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Data Graph
"""   

# standard modules
from collections import OrderedDict as odict
# scientific modules
import numpy as np
import scipy as sp
import matplotlib
from .compat.matplotlib import pyplot as pl
# scanpy modules
from . import settings as sett
from . import plotting as plott
from . import utils

class DataGraph(object):
    """ 
    Class for treating graphs.
    """

    def __init__(self, ddata_or_X_or_Dsq, params):
        """ 
        """
        isddata = isinstance(ddata_or_X_or_Dsq, dict)
        if isddata:
            ddata = ddata_or_X_or_Dsq
        else:
            X_or_Dsq = ddata_or_X_or_Dsq
        if isddata:
            self.Dsq = None
            if 'Xpca' in ddata:
                self.X = ddata['Xpca']
                sett.m(0, '--> using Xpca for building graph')
            else:
                self.X = ddata['X']
                sett.m(0, '--> using X for building graph')
        else:
            if X_or_Dsq.shape[0] == X_or_Dsq.shape[1]:
                sett.m(0,'--> computing data graph from distance matrix')
                self.Dsq = X_or_Dsq
                self.X = None
            else:
                sett.m(0,'--> computing data graph from data matrix')
                self.X = X_or_Dsq
                self.Dsq = None
        self.N = self.X.shape[0]
        self.params = params
        if self.params['sigma'] > 0:
            self.params['method'] = 'global'
        else:
            self.params['method'] = 'local'            

    def diffmap(self):
        """
        Diffusion Map as of Coifman et al. (2005) incorparting
        suggestions of Haghverdi et al. (2016).
        """
        sett.mt(0, 'start computing Diffusion Map with method',
                 '"'+self.params['method']+'"')
        self.compute_transition_matrix()
        self.embed()
        # write results to dictionary
        ddmap = {}
        # skip the first eigenvalue/eigenvector
        ddmap['Y'] = self.rbasis[:, 1:]
        ddmap['evals'] = self.evals[1:]
        return ddmap

    def spec_layout(self):
        """
        """
        self.compute_transition_matrix()
        self.compute_L_matrix()
        self.embed(self.L, sort='increase')
        # write results to dictionary
        ddmap = {}
        # skip the first eigenvalue/eigenvector
        ddmap['Y'] = self.rbasis[:, 1:]
        ddmap['evals'] = self.evals[1:]
        return ddmap

    def compute_transition_matrix(self, weighted=True, 
                                  neglect_selfloops=False, alpha=1):
        """ 
        Compute similarity matrix and transition matrix.

        Notes
        -----
        In the code, the following two parameters are set.

        alpha : float
            The density rescaling parameter of Coifman and Lafon (2006). Should
            in all practical applications equal 1: Then only the geometry of the
            data matters, not the sampling density.
        neglect_selfloops : bool 
            Discard selfloops.

        See also
        --------
        Also Haghverdi et al. (2016, 2015) and Coifman and Lafon (2006) and
        Coifman et al. (2005).
        """
        if self.Dsq is None:
            if self.X.shape[1] > 1000:
                sett.m(0, '--> high number of dimensions for computing distance matrix\n'
                       '    consider preprocessing using PCA')
            # compute distance matrix in squared Euclidian norm
            Dsq = utils.comp_distance(self.X, metric='sqeuclidean')
        else:
            Dsq = self.Dsq
        if self.params['method'] == 'local':
            # choose sigma (width of a Gaussian kernel) according to the
            # distance of the kth nearest neighbor of each point, including the
            # point itself in the count
            k = self.params['k']
            # deterimine the distance of the k nearest neighbors
            indices = np.zeros((Dsq.shape[0], k), dtype=np.int_)
            distances_sq = np.zeros((Dsq.shape[0], k), dtype=np.float_)
            for irow, row in enumerate(Dsq):
                # the last item is already in its sorted position as
                # argpartition puts the (k-1)th element - starting to count from
                # zero - in its sorted position
                idcs = np.argpartition(row, k-1)[:k]
                indices[irow] = idcs
                distances_sq[irow] = np.sort(row[idcs])
            # choose sigma, the heuristic here often makes not much 
            # of a difference, but is used to reproduce the figures
            # of Haghverdi et al. (2016)
            if self.params['knn']:
                # as the distances are not sorted except for last element
                # take median
                sigmas_sq = np.median(distances_sq, axis=1)
            else:
                # the last item is already in its sorted position as
                # argpartition puts the (k-1)th element - starting to count from
                # zero - in its sorted position
                sigmas_sq = distances_sq[:,-1]/4
            sigmas = np.sqrt(sigmas_sq)
            sett.mt(0, 'determined k =', k, 'nearest neighbors of each point')
        elif self.params['method'] == 'standard':
            sigmas = self.params['sigma'] * np.ones(self.N)
            sigmas_sq = sigmas**2

        # compute the symmetric weight matrix
        Num = 2 * np.multiply.outer(sigmas,sigmas)
        Den = np.add.outer(sigmas_sq,sigmas_sq)
        W = np.sqrt(Num/Den) * np.exp(-Dsq/Den)
        # make the weight matrix sparse
        if not self.params['knn']:
            self.Mask = W > 1e-14
            W[self.Mask == False] = 0
        else:
            # restrict number of neighbors to k
            Mask = np.zeros(Dsq.shape, dtype=bool)
            for irow,row in enumerate(indices):
                Mask[irow, row] = True
                for j in row:
                    if irow not in indices[j]:
                        Mask[j, irow] = True
            # set all entries that are not nearest neighbors to zero
            W[Mask == False] = 0
            self.Mask = Mask

        if not weighted:
            W = Mask.astype(float)
        sett.mt(0,'computed W (weight matrix) with "knn" =', self.params['knn'])

        # neglect self-loops
        # also proposed by Haghverdi et al. (2015)
        # notice that then, the kernel does not encode a notion of similarity
        # then anymore and is not positive semidefinite anymore in practice, it
        # doesn't matter too much
        if neglect_selfloops:
            np.fill_diagonal(W, 0)
            
        if False:
            pl.matshow(W)
            pl.title('$ W$')
            pl.colorbar()
            
        # density normalization 
        # as discussed in Coifman et al. (2005)
        # ensure that kernel matrix is independent of sampling density
        if alpha == 0:
            # nothing happens here, simply use the isotropic similarity matrix
            self.K = np.array(W)
        else:
            # q[i] is an estimate for the sampling density at point x_i
            # it's also the degree of the underlying graph
            q = np.sum(W, axis=0)
            # raise to power alpha
            if alpha != 1: 
                q = q**alpha
            Den = np.outer(q, q)
            self.K = W / Den
        sett.mt(0,'computed K (anisotropic kernel)')
        if False:
            pl.matshow(self.K)
            pl.title('$ K$')
            pl.colorbar()

        # now compute the row normalization to build the transition matrix T
        # and the adjoint Ktilde: both have the same spectrum
        self.z = np.sum(self.K, axis=0)
        # the following is the transition matrix
        self.T = self.K/self.z[:, np.newaxis]
        # now we need the square root of the density
        self.sqrtz = np.array(np.sqrt(self.z))
        # now compute the density-normalized Kernel
        # it's still symmetric
        szszT = np.outer(self.sqrtz, self.sqrtz)
        self.Ktilde = self.K / szszT
        sett.mt(0,'computed Ktilde (normalized anistropic kernel)')

        if False:
            pl.matshow(self.Ktilde)
            pl.title('$ \widetilde K$')
            pl.colorbar()
            pl.show()

    def compute_L_matrix(self):
        """
        Graph Laplacian for K.
        """
        self.L = np.diag(self.z) - self.K

    def embed(self, matrix=None, number=10, sym=True, sort='decrease'):
        """ 
        Compute eigen decomposition of matrix.

        Parameters
        ----------
        matrix : np.ndarray
            Matrix to diagonalize.
        number : int  
            Number of eigenvalues/vectors to be computed, set number = 0 if
            you need all eigenvectors.
        sym : bool
            Instead of computing the eigendecomposition of the assymetric 
            transition matrix, computed the eigendecomposition of the symmetric
            Ktilde matrix.

        Writes class members
        --------------------
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
        np.set_printoptions(precision=3)
        self.rbasisBool = True
        if matrix is None:
            matrix = self.Ktilde
        # compute the spectrum
        if number == 0:
            evals, evecs = sp.linalg.eigh(matrix)
        else:
            number = min(matrix.shape[0]-1, number)
            # ncv = max(2 * number + 1, int(np.sqrt(matrix.shape[0])))
            ncv = None
            which = 'LM' if sort == 'decrease' else 'SM'
            evals, evecs = sp.sparse.linalg.eigsh(matrix, k=number, which=which,
                                                  ncv=ncv)
        if sort == 'decrease':
            evals = evals[::-1]
            evecs = evecs[:, ::-1]
        sett.mt(0, 'computed eigenvalues:')
        sett.m(0, evals)
        sett.m(1, 'computed', number, 'eigenvalues. if you want more increase the'
               'parameter "number" or set it to zero, to compute all eigenvalues')
        # assign attributes
        self.evals = evals
        if sym:
            self.rbasis = self.lbasis = evecs
        else:
            # The eigenvectors of T are stored in self.rbasis and self.lbasis 
            # and are simple trafos of the eigenvectors of Ktilde.
            # rbasis and lbasis are right and left eigenvectors, respectively
            self.rbasis = np.array(evecs/self.sqrtz[:,np.newaxis])
            self.lbasis = np.array(evecs*self.sqrtz[:,np.newaxis])
            # normalize in L2 norm
            # note that, in contrast to that, a probability distribution
            # on the graph is normalized in L1 norm
            # therefore, the eigenbasis in this normalization does not correspond
            # to a probability distribution on the graph
            self.rbasis /= np.linalg.norm(self.rbasis,axis=0,ord=2)
            self.lbasis /= np.linalg.norm(self.lbasis,axis=0,ord=2)

    def compute_M_matrix(self):
        """ 
        The M matrix is the matrix that results from summing over all powers of
        T in the subspace without the first eigenspace.

        See Haghverdi et al. (2016).
        """
        # the projected inverse therefore is
        self.M = sum([self.evals[i]/(1-self.evals[i])
                      * np.outer(self.rbasis[:, i], self.lbasis[:, i])
                      for i in range(1, self.evals.size)])
        sett.mt(0,'computed M matrix')
        if False:
            pl.matshow(self.Ktilde)
            pl.title('Ktilde')
            pl.colorbar()
            pl.matshow(self.M)
            pl.title('M')
            pl.colorbar()

    def compute_Lp_matrix(self):
        """
        See Fouss et al. (2006) and von Luxburg et al. (2007).

        See Proposition 6 in von Luxburg (2007) and the inline equations
        right in the text above.
        """
        self.Lp = sum([1/self.evals[i]
                      * np.outer(self.rbasis[:,i], self.lbasis[:,i])
                      for i in range(1, self.evals.size)])
        sett.mt(0,'computed pseudoinverse of Laplacian')

    def compute_Ddiff_matrix(self):
        """ 
        Returns the distance matrix in the Diffusion Pseudotime metric.

        See Haghverdi et al. (2016).
        
        Notes
        -----
        - Is based on the M matrix. 
        - self.Ddiff[self.iroot,:] stores diffusion pseudotime as a vector.
        """
        if self.M.shape[0] > 1000 and self.params['nr_pcs'] == 0:
            sett.m(0, '--> high number of dimensions for computing DPT distance matrix\n'
                   '    by setting nr_pcs > 0 you can speed up the computation')
        if self.params['nr_pcs'] > 0 and self.M.shape[0] > self.params['nr_pcs']:
            import scanpy.preprocess as pp
            self.M = pp(self.M, 'pca', nr_comps=self.params['nr_pcs'])
        self.Ddiff = sp.spatial.distance.pdist(self.M)
        self.Ddiff = sp.spatial.distance.squareform(self.Ddiff)
        sett.mt(0, 'computed Ddiff distance matrix')
        self.Dchosen = self.Ddiff

    def compute_C_matrix(self):
        """
        See Fouss et al. (2006) and von Luxburg et al. (2007).

        This is the commute-time matrix. It's a squared-euclidian distance
        matrix in \mathbb{R}^n.
        """
        self.C = np.zeros(self.Lp.shape)
        for i in range(self.Lp.shape[0]):
            for j in range(self.Lp.shape[1]):
                self.C[i, j] = self.Lp[i, i] + self.Lp[j, j] - 2*self.Lp[i, j]
        volG = np.sum(self.z)
        self.C *= volG
        sett.mt(0,'computed commute distance matrix')
        self.Dchosen =  self.C

    def compute_MFP_matrix(self):
        """
        See Fouss et al. (2006).

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
        sett.mt(0,'computed mean first passage time matrix')
        self.Dchosen = self.MFP

    def set_pseudotimes(self):
        """
        Return pseudotimes with respect to root point.
        """
        self.pseudotimes = self.Dchosen[self.iroot]/np.max(self.Dchosen[self.iroot])

    def find_root(self,xroot):
        """ 
        Determine the index of the root cell.

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
                             'reduced version, if you provided Xpca.')
        # this is the squared distance
        dsqroot = 1e10
        self.iroot = 0
        for i in range(self.N):
            diff = self.X[i, :]-xroot
            dsq = diff.dot(diff)
            if  dsq < dsqroot:
                dsqroot = dsq
                self.iroot = i
                if np.sqrt(dsqroot) < 1e-10:
                    sett.m(2,'root found at machine prec')
                    break
        sett.m(1,'sample',self.iroot,'has distance',np.sqrt(dsqroot),'from root')
        return self.iroot

    def _test_embed(self):
        """
        Checks and tests for embed.
        """
        # pl.semilogy(w,'x',label=r'$ \widetilde K$')
        # pl.show()
        if sett.verbosity > 2:
            # output of spectrum of K for comparison
            w,v = np.linalg.eigh(self.K)
            sett.mi('spectrum of K (kernel)')
        if sett.verbosity > 3:
            # direct computation of spectrum of T
            w,vl,vr = sp.linalg.eig(self.T,left=True)
            sett.mi('spectrum of transition matrix (should be same as of Ktilde)')

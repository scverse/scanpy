# Author: F. Alex Wolf (http://falexwolf.de)
"""Diffusion Pseudotime Analysis
"""

import sys
import numpy as np
import scipy as sp
import scipy.sparse
from .. import logging as logg
from ..data_structs import data_graph


def dpt(adata, n_branchings=0, k=30, knn=True, n_pcs=50, n_pcs_post=30,
        allow_branching_at_root=False, n_jobs=None, recompute_diffmap=False,
        flavor='haghverdi16', copy=False):
    """Diffusion Pseudotime

    Infer progression of cells, identify branching subgroups.

    Reference
    ---------
    - Diffusion Pseudotime: Haghverdi et al., Nature Methods 13, 3971 (2016).

    See also
    --------
    - Diffusion Maps: Coifman et al., PNAS 102, 7426 (2005).
    - Diffusion Maps applied to single-cell data: Haghverdi et al., Bioinformatics
      31, 2989 (2015).
    - Diffusion Maps as a flavour of spectral clustering: von Luxburg,
      arXiv:0711.0189 (2007).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata.add['xroot'] : np.ndarray
            Root of stochastic process on data points (root cell), specified
            as expression vector of shape X.shape[1].
        adata.smp['X_pca']: np.ndarray
            PCA representation of the data matrix (result of preprocessing with
            PCA). If it exists in adata, dpt will use this instead of adata.X.
        adata.smp['X_diffmap']: np.ndarray
            Diffmap representation of the data matrix (result of running
            `diffmap`).  Will be used if option `recompute_diffmap` is False.
    n_branchings : int, optional (default: 1)
        Number of branchings to detect.
    k : int, optional (default: 30)
        Number of nearest neighbors on the knn graph. If knn == False, set the
        Gaussian kernel width to the distance of the kth neighbor.
    knn : bool, optional (default: True)
        If True, use a hard threshold to restrict the number of neighbors to
        k, that is, consider a knn graph. Otherwise, use a Gaussian Kernel
        to assign low weights to neighbors more distant than the kth nearest
        neighbor.
    n_pcs: int, optional (default: 50)
        Use n_pcs PCs to compute the Euclidian distance matrix, which is the
        basis for generating the graph. Set to 0 if you don't want preprocessing
        with PCA.
    n_pcs_post: int, optional (default: 30)
        Use n_pcs_post PCs to compute the DPT distance matrix. This speeds up
        the computation at almost no loss of accuracy. Set to 0 if you don't
        want postprocessing with PCA.
    allow_branching_at_root : bool, optional (default: False)
        Allow to have branching directly at root point.
    n_jobs : int or None (default: None)
        Number of cpus to use for parallel processing (default: sett.n_jobs).
    recompute_diffmap : bool, (default: False)
        Recompute diffusion maps.
    copy : bool, optional (default: False)
        Copy instance before computation and return a copy. Otherwise, perform
        computation inplace and return None.

    Notes
    -----
    Writes the following arrays as sample annotation to adata.smp.
        dpt_pseudotime : np.ndarray
            Array of dim (number of samples) that stores the pseudotime of each
            cell, that is, the DPT distance with respect to the root cell.
        dpt_groups : np.ndarray of dtype string
            Array of dim (number of samples) that stores the subgroup id ('0',
            '1', ...) for each cell. The groups  typically correspond to
            'progenitor cells', 'undecided cells' or 'branches' of a process.
    Writes the following additional arrays as unstructured annotation to adata.
        X_diffmap : np.ndarray
            Array of shape (number of samples) x (number of eigen
            vectors). DiffMap representation of data, which is the right eigen
            basis of the transition matrix with eigenvectors as columns.
        dpt_evals : np.ndarray
            Array of size (number of eigen vectors). Eigenvalues of transition matrix.
    """
    adata = adata.copy() if copy else adata
    if 'xroot' not in adata.add:
        msg = \
   '''DPT requires specifying the expression "xroot" of a root cell.

   Either
       adata.add['xroot'] = adata.X[root_cell_index, :]
   where "root_cell_index" is the integer index of the root cell, or
       adata.add['xroot'] = adata[root_cell_name, :].X
   where "root_cell_name" is the name (a string) of the root cell.'''
        sys.exit(msg)
    if n_branchings == 0:
        logg.m('set parameter `n_branchings` > 0 to detect branchings', v='hint')
    dpt = DPT(adata, k=k, knn=knn, n_pcs=n_pcs, n_pcs_post=n_pcs_post,
              n_jobs=n_jobs, recompute_diffmap=recompute_diffmap,
              n_branchings=n_branchings, allow_branching_at_root=allow_branching_at_root,
              flavor=flavor)
    # diffusion map
    ddmap = dpt.diffmap()
    adata.smp['X_diffmap'] = ddmap['X_diffmap']
    # also store the 0th comp, which is skipped for plotting
    adata.smp['X_diffmap0'] = dpt.rbasis[:, 0]
    adata.add['diffmap_evals'] = ddmap['evals']
    if knn: adata.add['distance'] = dpt.Dsq
    logg.m('perform Diffusion Pseudotime analysis', r=True)
    if False:
        # compute M matrix of cumulative transition probabilities,
        # see Haghverdi et al. (2016)
        dpt.compute_M_matrix()
    # compute DPT distance matrix, which we refer to as 'Ddiff'
    if False:  # we do not compute the full Ddiff matrix, only the elements we need
        dpt.compute_Ddiff_matrix()
    dpt.set_pseudotime()  # pseudotimes are distances from root point
    adata.add['iroot'] = dpt.iroot  # update iroot, might have changed when subsampling, for example
    adata.smp['dpt_pseudotime'] = dpt.pseudotime
    # detect branchings and partition the data into segments
    dpt.branchings_segments()
    # vector of length n_groups
    adata.add['dpt_groups_names'] = np.array(dpt.segs_names_unique, dtype='U')
    for itips, tips in enumerate(dpt.segs_tips):
        if tips[0] == -1: adata.add['dpt_groups_names'][itips] = '?'
        if dpt.segs_undecided[itips]: adata.add['dpt_groups_names'][itips] = '?'
    # vector of length n_samples of groupnames
    adata.smp['dpt_groups'] = dpt.segs_names.astype('U')
    # the ordering according to segments and pseudotime
    adata.smp['dpt_order'] = dpt.indices
    # the changepoints - marking different segments - in the ordering above
    adata.add['dpt_changepoints'] = dpt.changepoints
    # the tip points of segments
    adata.add['dpt_grouptips'] = dpt.segs_tips
    # the tree/graph adjacency matrix
    adata.add['dpt_groups_adjacency'] = dpt.segs_adjacency
    logg.m('finished', t=True, end=' ')
    logg.m('and added\n'
           '    "dpt_pseudotime", stores pseudotime (adata.smp),\n'
           '    "dpt_groups", the segments of the tree-like trajectory (adata.smp),\n'
           '    "dpt_order", is an index array for sorting the cells (adata.smp),\n'
           '    "dpt_grouptips", stores the indices of tip cells (adata.add)')
    return adata if copy else None


class DPT(data_graph.DataGraph):
    """Diffusion Pseudotime.
    """

    def __init__(self, adata_or_X, k=30, knn=True,
                 n_jobs=1, n_pcs=50, n_pcs_post=30,
                 recompute_diffmap=None, n_branchings=0,
                 allow_branching_at_root=False, flavor='haghverdi16'):
        super(DPT, self).__init__(adata_or_X, k=k, knn=knn, n_pcs=n_pcs,
                                  n_pcs_post=n_pcs_post, n_jobs=n_jobs,
                                  recompute_diffmap=recompute_diffmap,
                                  flavor=flavor)
        self.n_branchings = n_branchings
        self.allow_branching_at_root = allow_branching_at_root

    def branchings_segments(self):
        """Detect branchings and partition the data into corresponding segments.

        Detect all branchings up to `n_branchings`.

        Writes
        ------
        segs : np.ndarray
            Array of dimension (number of segments) x (number of data
            points). Each row stores a mask array that defines a segment.
        segs_tips : np.ndarray
            Array of dimension (number of segments) x 2. Each row stores the
            indices of the two tip points of each segment.
        segs_names : np.ndarray
            Array of dimension (number of data points). Stores an integer label
            for each segment.
        """
        # logg.m('weights', self.evals[1:]/(1-self.evals[1:]))
        # logg.m('evals', self.evals[1:])
        self.detect_branchings()
        self.check_segments()
        self.postprocess_segments()
        # self.order_segments()
        self.set_segs_names()
        self.order_pseudotime()

    def detect_branchings(self):
        """Detect all branchings up to `n_branchings`.

        Writes Attributes
        -----------------
        segs : np.ndarray
            List of integer index arrays.
        segs_tips : np.ndarray
            List of indices of the tips of segments.
        """
        logg.m('... detect', self.n_branchings, 'branchings')
        # a segment is a subset of points of the data set (defined by the
        # indices of the points in the segment)
        # initialize the search for branchings with a single segment,
        # that is, get the indices of the whole data set
        indices_all = np.arange(self.X.shape[0], dtype=int)
        # let's keep a list of segments, the first segment to add is the
        # whole data set
        segs = [indices_all]
        # a segment can as well be defined by the two points that have maximal
        # distance in the segment, the "tips" of the segment
        #
        # the rest of the points in the segment is then defined by demanding
        # them to "be close to the line segment that connects the tips", that
        # is, for such a point, the normalized added distance to both tips is
        # smaller than one:
        #     (D[tips[0],i] + D[tips[1],i])/D[tips[0],tips[1] < 1
        # of course, this condition is fulfilled by the full cylindrical
        # subspace surrounding that line segment, where the radius of the
        # cylinder can be infinite
        #
        # if D denotes a euclidian distance matrix, a line segment is a linear
        # object, and the name "line" is justified. if we take the
        # diffusion-based distance matrix Dchosen, which approximates geodesic
        # distance, with "line", we mean the shortest path between two points,
        # which can be highly non-linear in the original space
        #
        # let us define the tips of the whole data set
        if False:  # this is safe, but not compatible with on-the-fly computation
            tips_all = np.array(np.unravel_index(np.argmax(self.Dchosen), self.Dchosen.shape))
        else:
            tip_0 = np.argmax(self.Dchosen[self.iroot])
            tips_all = np.array([tip_0, np.argmax(self.Dchosen[tip_0])])
        # we keep a list of the tips of each segment
        segs_tips = [tips_all]
        segs_undecided = [True]
        segs_adjacency = [[]]
        for ibranch in range(self.n_branchings):
            # this is dangerous!
            # make sure that in each iteration the correct Dchosen
            # is used, also, Dchosen saves elements...,
            # calling OnFlySymMatrix destroys that, can we save something?
            # if ibranch > 0:
            #     DC_start = ibranch + 1
            #     logg.m('... constraining to DC >=', DC_start)
            #     self.Dchosen = data_graph.OnFlySymMatrix(self.get_Ddiff_row,
            #                                              shape=self.Dchosen.shape,
            #                                              DC_start=DC_start)
            # out of the list of segments, determine the segment
            # that most strongly deviates from a straight line
            # and provide the three tip points that span the triangle
            # of maximally distant points
            iseg, tips3 = self.select_segment(segs, segs_tips, segs_undecided)
            logg.m('... detect branching {}:'.format(ibranch + 1),
                   'split group', iseg,
                   'with tip cells [third start end] =', tips3, t=True)
            # detect branching and update segs and segs_tips
            self.detect_branching(segs, segs_tips, segs_undecided,
                                  segs_adjacency, iseg, tips3)
        # get back to what we had in the beginning
        # self.Dchosen = data_graph.OnFlySymMatrix(self.get_Ddiff_row,
        #                                          shape=self.Dchosen.shape)
        # store as class members
        self.segs = segs
        self.segs_tips = segs_tips
        self.segs_undecided = segs_undecided
        # the following is a bit too much, but this allows easy storage
        self.segs_adjacency = np.zeros((len(segs), len(segs)), dtype=bool)
        for i, seg_adjacency in enumerate(segs_adjacency):
            self.segs_adjacency[i, seg_adjacency] = 1

    def select_segment(self, segs, segs_tips, segs_undecided):
        """Out of a list of line segments, choose segment that has the most
        distant second data point.

        Assume the distance matrix Ddiff is sorted according to seg_idcs.
        Compute all the distances.

        Returns
        -------
        iseg : int
            Index identifying the position within the list of line segments.
        tips3 : int
            Positions of tips within chosen segment.
        """
        scores_tips = np.zeros((len(segs), 4))
        allindices = np.arange(self.X.shape[0], dtype=int)
        for iseg, seg in enumerate(segs):
            # do not consider too small segments
            if segs_tips[iseg][0] == -1: continue
            # restrict distance matrix to points in segment
            if not isinstance(self.Dchosen, data_graph.OnFlySymMatrix):
                Dseg = self.Dchosen[np.ix_(seg, seg)]
            else:
                Dseg = self.Dchosen.restrict(seg)
            third_maximizer = None
            if segs_undecided[iseg]:
                # check that none of our tips "connects" with a tip of the
                # other segments
                for jseg in range(len(segs)):
                    if jseg != iseg:
                        # take the inner tip, the "second tip" of the segment
                        for itip in range(2):
                            if (self.Dchosen[segs_tips[jseg][1], segs_tips[iseg][itip]]
                                < 0.5 * self.Dchosen[segs_tips[iseg][~itip], segs_tips[iseg][itip]]):
                                logg.m('... group', iseg, 'with tip', segs_tips[iseg][itip],
                                       'connects with', jseg, 'with tip', segs_tips[jseg][1])
                                logg.m('    do not use the tip for "triangulation"')
                                third_maximizer = itip
            # map the global position to the position within the segment
            tips = [np.where(allindices[seg] == tip)[0][0]
                    for tip in segs_tips[iseg]]
            # find the third point on the segment that has maximal
            # added distance from the two tip points
            dseg = Dseg[tips[0]] + Dseg[tips[1]]
            # add this point to tips, it's a third tip, we store it at the first
            # position in an array called tips3
            third_tip = np.argmax(dseg)
            if third_maximizer is not None:
                # find a fourth point that has maximal distance to all three
                dseg += Dseg[third_tip]
                fourth_tip = np.argmax(dseg)
                tips[itip] = fourth_tip
                dseg -= Dseg[tips[itip]]
            tips3 = np.insert(tips, 0, third_tip)
            # compute the score as ratio of the added distance to the third tip,
            # to what it would be if it were on the straight line between the
            # two first tips, given by Dseg[tips[:2]]
            # if we did not normalize, there would be a danger of simply
            # assigning the highest score to the longest segment
            score = dseg[tips3[0]] / Dseg[tips3[1], tips3[2]]
            print('... group', iseg, 'score', score, 'n_points', len(seg))
            if len(seg) < 50: score = 0
            # logg.m('... do not consider groups with less than 50 points for splitting')
            # write result
            scores_tips[iseg, 0] = score
            scores_tips[iseg, 1:] = tips3
        iseg = np.argmax(scores_tips[:, 0])
        tips3 = scores_tips[iseg, 1:].astype(int)
        return iseg, tips3

    def postprocess_segments(self):
        """Convert the format of the segment class members."""
        # make segs a list of mask arrays, it's easier to store
        # as there is a hdf5 equivalent
        for iseg, seg in enumerate(self.segs):
            mask = np.zeros(self.X.shape[0], dtype=bool)
            mask[seg] = True
            self.segs[iseg] = mask
        # convert to arrays
        self.segs = np.array(self.segs)
        self.segs_tips = np.array(self.segs_tips)

    def check_segments(self):
        """Perform checks on segments and sort them according to pseudotime."""
        # find the segment that contains the root cell
        for iseg, seg in enumerate(self.segs):
            if self.iroot in seg:
                isegroot = iseg
                break
        # check whether the root cell is one of the tip cells of the
        # segment, if not we need to introduce a new branching, directly
        # at the root cell
        if self.iroot not in self.segs_tips[isegroot]:
            # if it's not exactly a tip, but very close to it,
            # just keep it as it is
            dist_to_root = self.Dchosen[self.iroot, self.segs_tips[iseg]]
            # otherwise, allow branching at root
            if (False and np.min(dist_to_root) > 0.01*self.Dchosen[tuple(self.segs_tips[iseg])]
                and self.allow_branching_at_root):
                logg.m('... adding branching directly at root')
                allindices = np.arange(self.X.shape[0], dtype=int)
                tips3_global = np.insert(self.segs_tips[iseg], 0, self.iroot)
                # map the global position to the position within the segment
                tips3 = np.array([np.where(allindices[self.segs[iseg]] == tip)[0][0]
                                  for tip in tips3_global])
                # detect branching and update self.segs and self.segs_tips
                self.segs, self.segs_tips = self.detect_branching(self.segs,
                                                                 self.segs_tips,
                                                                 iseg, tips3)

    def order_segments(self):
        """Order segments according to average pseudotime."""
        # there are different options for computing the score
        if False:
            # minimum of pseudotime in the segment
            score = np.min
        if True:
            # average pseudotime
            score = np.average
        # score segments by minimal pseudotime
        seg_scores = []
        for seg in self.segs:
            seg_scores.append(score(self.pseudotime[seg]))
        indices = np.argsort(seg_scores)
        # order segments by minimal pseudotime
        self.segs = self.segs[indices]
        self.segs_tips = self.segs_tips[indices]

    def set_segs_names(self):
        """Return a single array that stores integer segment labels."""
        segs_names = np.zeros(self.X.shape[0], dtype=np.int8)
        self.segs_names_unique = []
        for iseg, seg in enumerate(self.segs):
            segs_names[seg] = iseg
            self.segs_names_unique.append(iseg)
        self.segs_names = segs_names

    def order_pseudotime(self):
        """Define indices that reflect segment and pseudotime order.

        Writes
        ------
        indices : np.ndarray
            Index array of shape n, which stores an ordering of the data points
            with respect to increasing segment index and increasing pseudotime.
        changepoints : np.ndarray
            Index array of shape len(ssegs)-1, which stores the indices of
            points where the segment index changes, with respect to the ordering
            of indices.
        """
        # within segs_tips, order tips according to pseudotime
        for itips, tips in enumerate(self.segs_tips):
            if tips[0] != -1:
                indices = np.argsort(self.pseudotime[tips])
                self.segs_tips[itips] = self.segs_tips[itips][indices]
            else:
                logg.m('... group', itips, 'is very small')
        # sort indices according to segments
        indices = np.argsort(self.segs_names)
        segs_names = self.segs_names[indices]
        # find changepoints of segments
        changepoints = np.arange(indices.size-1)[np.diff(segs_names) == 1] + 1
        pseudotime = self.pseudotime[indices]
        for iseg, seg in enumerate(self.segs):
            # only consider one segment, it's already ordered by segment
            seg_sorted = seg[indices]
            # consider the pseudotime on this segment and sort them
            seg_indices = np.argsort(pseudotime[seg_sorted])
            # within the segment, order indices according to increasing pseudotime
            indices[seg_sorted] = indices[seg_sorted][seg_indices]
        # define class members
        self.indices = indices
        self.changepoints = changepoints

    def detect_branching(self, segs, segs_tips, segs_undecided, segs_adjacency, iseg, tips3):
        """Detect branching on given segment.

        Call function _detect_branching and perform bookkeeping on segs and
        segs_tips.

        Parameters
        ----------
        segs : list of np.ndarray
            Dchosen distance matrix restricted to segment.
        segs_tips : list of np.ndarray
            Stores all tip points for the segments in segs.
        iseg : int
            Position of segment under study in segs.
        tips3 : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.

        Returns
        -------
        segs : list of np.ndarray
            Updated list of segments.
        segs_tips : list of np.ndarray
            Updated list of segs_tips.
        """
        seg = segs[iseg]
        # restrict distance matrix to points in segment
        if not isinstance(self.Dchosen, data_graph.OnFlySymMatrix):
            Dseg = self.Dchosen[np.ix_(seg, seg)]
        else:
            Dseg = self.Dchosen.restrict(seg)
        # given the three tip points and the distance matrix detect the
        # branching on the segment, return the list ssegs of segments that
        # are defined by splitting this segment
        ssegs, ssegs_tips = self._detect_branching(Dseg, tips3)
        # map back to global indices
        for iseg_new, seg_new in enumerate(ssegs):
            if len(np.flatnonzero(seg_new)) > 1:  # terrible hack
                ssegs[iseg_new] = seg[seg_new]
                if ssegs_tips[iseg_new][0] != -1:
                    ssegs_tips[iseg_new] = seg[ssegs_tips[iseg_new]]
        # remove previous segment
        segs.pop(iseg)
        segs_tips.pop(iseg)
        segs_undecided.pop(iseg)
        # insert undecided cells at same position
        segs.insert(iseg, ssegs[-1])
        segs_tips.insert(iseg, ssegs_tips[-1])
        segs_undecided.insert(iseg, True)
        # append new segments
        segs += ssegs[:-1]
        segs_tips += ssegs_tips[:-1]
        segs_undecided += [False, False, False]
        # establish edges
        # each of the new segments connects with the undecided cells
        segs_adjacency[iseg] += list(range(len(segs_adjacency), len(segs_adjacency) + 3))
        segs_adjacency += [[iseg], [iseg], [iseg]]
        # check whether one of the tips of the undecided cells
        # connects to the previous segments, or whether one
        # of the new 'proper segments' connects with the previous cells
        for jseg in range(len(segs) - 3, len(segs)):
            for kseg in range(len(segs) - 3):
                if kseg != iseg:
                    for jtip in range(2):
                        for ktip in range(2):
                            j_index = segs_tips[jseg][jtip]
                            k_index = segs_tips[kseg][ktip]
                            if ((sp.sparse.issparse(self.Dsq) and self.Dsq[j_index, k_index] > 0)
                                or (not sp.sparse.issparse(self.Dsq)
                                    and self.Dsq[j_index, k_index]
                                    <= 5 * np.partition(self.Dsq[j_index], self.k-1)[self.k-1])):
                                segs_adjacency[jseg].append(kseg)
                                segs_adjacency[kseg].append(jseg)
                                logg.m('... adding edge between groups', jseg, kseg)
                                try:
                                    segs_adjacency[iseg].remove(kseg)
                                    segs_adjacency[kseg].remove(iseg)
                                except ValueError:
                                    pass
        return segs, segs_tips, segs_undecided, segs_adjacency

    def _detect_branching(self, Dseg, tips):
        """Detect branching on given segment.

        Call function __detect_branching three times for all three orderings of
        tips. Points that do not belong to the same segment in all three
        orderings are assigned to a fourth segment. The latter is, by Haghverdi
        et al. (2016) referred to as 'undecided cells'.

        Parameters
        ----------
        Dseg : np.ndarray
            Dchosen distance matrix restricted to segment.
        tips : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.

        Returns
        -------
        ssegs : list of np.ndarray
            List of segments obtained from splitting the single segment defined
            via the first two tip cells.
        ssegs_tips : list of np.ndarray
            List of tips of segments in ssegs.
        """
        if True:
            ssegs = self._detect_branching_single(Dseg, tips)
        else:  # not needed, just for conserving the idea
            ssegs = self._detect_branching_versions(Dseg, tips)
        # make sure that each data point has a unique association with a segment
        masks = np.zeros((3, Dseg.shape[0]), dtype=bool)
        for iseg, seg in enumerate(ssegs):
            masks[iseg][seg] = True
        nonunique = np.sum(masks, axis=0) > 1
        #
        ssegs = []
        for iseg, mask in enumerate(masks):
            mask[nonunique] = False
            ssegs.append(np.arange(Dseg.shape[0], dtype=int)[mask])
        # compute new tips within new segments
        ssegs_tips = []
        for inewseg, newseg in enumerate(ssegs):
            if len(np.flatnonzero(newseg)) > 1:
                secondtip = newseg[np.argmax(Dseg[tips[inewseg]][newseg])]
                ssegs_tips.append([tips[inewseg], secondtip])
            else:
                ssegs_tips.append(np.array([-1, -1]))
                # logg.m('... group', inewseg, 'contains less than 4 data points')
        # add the points not associated with a clear seg to ssegs
        mask = np.zeros(Dseg.shape[0], dtype=bool)
        # all points assigned to segments (flatten ssegs)
        mask[[i for l in ssegs for i in l]] = True
        # append all the points that have not been assigned, in Haghverdi et
        # al. (2016), we call them 'undecided cells'
        undecided_cells = np.arange(Dseg.shape[0], dtype=int)[mask == False]
        ssegs.append(undecided_cells)
        tip_0 = undecided_cells[np.argmax(Dseg[undecided_cells[0]][undecided_cells])]
        tip_1 = undecided_cells[np.argmax(Dseg[tip_0][undecided_cells])]
        ssegs_tips.append([tip_0, tip_1])
        # ssegs_tips.append([-1, -1])
        return ssegs, ssegs_tips

    def _detect_branching_single(self, Dseg, tips):
        """Detect branching on given segment.
        """
        # compute branchings using different starting points the first index of
        # tips is the starting point for the other two, the order does not
        # matter
        ssegs = []
        # permutations of tip cells
        ps = [[0, 1, 2],  # start by computing distances from the first tip
              [1, 2, 0],  #             -"-                       second tip
              [2, 0, 1]]  #             -"-                       third tip
        for i, p in enumerate(ps):
            ssegs.append(self.__detect_branching(Dseg, tips[p]))
        return ssegs

    def _detect_branching_versions(self, Dseg, tips):
        """Detect branching on given segment using three different versions.

        Did not prove useful, currently.
        """
        # compute branchings using different starting points the first index of
        # tips is the starting point for the other two, the order does not
        # matter
        ssegs_versions = []
        # permutations of tip cells
        ps = [[0, 1, 2],  # start by computing distances from the first tip
              [1, 2, 0],  #             -"-                       second tip
              [2, 0, 1]]  #             -"-                       third tip
        # invert permutations
        inv_ps = [[0, 1, 2],
                  [2, 0, 1],
                  [1, 2, 0]]
        for i, p in enumerate(ps):
            ssegs = self.__detect_branching(Dseg,
                                            tips[p])
            ssegs_versions.append(np.array(ssegs)[inv_ps[i]])
        ssegs = []
        # run through all three assignments of segments, and keep
        # only those assignments that were found in all three runs
        for inewseg, newseg_versions in enumerate(np.array(ssegs_versions).T):
            if len(newseg_versions) == 3:
                newseg = np.intersect1d(np.intersect1d(newseg_versions[0],
                                                       newseg_versions[1]),
                                        newseg_versions[2])
            else:
                newseg = newseg_versions[0]
            ssegs.append(newseg)
        return ssegs

    def __detect_branching(self, Dseg, tips):
        """Detect branching on given segment.

        Compute point that maximizes kendall tau correlation of the sequences of
        distances to the second and the third tip, respectively, when 'moving
        away' from the first tip: tips[0]. 'Moving away' means moving in the
        direction of increasing distance from the first tip.

        Parameters
        ----------
        Dseg : np.ndarray
            Dchosen distance matrix restricted to segment.
        tips : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.

        Returns
        -------
        ssegs : list of np.ndarray
            List of segments obtained from "splitting away the first tip cell".
        """
        # sort distance from first tip point
        # then the sequence of distances Dseg[tips[0]][idcs] increases
        idcs = np.argsort(Dseg[tips[0]])
        # consider now the sequence of distances from the other
        # two tip points, which only increase when being close to `tips[0]`
        # where they become correlated
        # at the point where this happens, we define a branching point
        if True:
            imax = self.kendall_tau_split(Dseg[tips[1]][idcs],
                                          Dseg[tips[2]][idcs])
        if False:
            # if we were in euclidian space, the following should work
            # as well, but here, it doesn't because the scales in Dseg are
            # highly different, one would need to write the following equation
            # in terms of an ordering, such as exploited by the kendall
            # correlation method above
            imax = np.argmin(Dseg[tips[0]][idcs]
                             + Dseg[tips[1]][idcs]
                             + Dseg[tips[2]][idcs])
        # init list to store new segments
        ssegs = []
        # first new segment: all points until, but excluding the branching point
        # increasing the following slightly from imax is a more conservative choice
        # as the criterion based on normalized distances, which follows below,
        # is less stable
        ibranch = imax + 2  # this used to be imax + 1!
        # ibranch = int(0.95 * imax)
        return idcs[:ibranch]
        # ssegs.append(idcs[:ibranch])
        # TODO get rid of the following heuristics
        # define nomalized distances to tip points for the rest of the data
        # dist1 = Dseg[tips[1], idcs[ibranch:]] / Dseg[tips[1], idcs[ibranch-1]]
        # dist2 = Dseg[tips[2], idcs[ibranch:]] / Dseg[tips[2], idcs[ibranch-1]]
        # assign points according to whether being closer to tip cell 1 or 2
        # ssegs.append(idcs[ibranch:][dist1 <= dist2])
        # ssegs.append(idcs[ibranch:][dist1 > dist2])
        # return ssegs

    def kendall_tau_split(self, a, b):
        """Return splitting index that maximizes correlation in the sequences.

        Compute difference in Kendall tau for all splitted sequences.

        For each splitting index i, compute the difference of the two
        correlation measures kendalltau(a[:i], b[:i]) and
        kendalltau(a[i:], b[i:]).

        Returns the splitting index that maximizes
            kendalltau(a[:i], b[:i]) - kendalltau(a[i:], b[i:])

        Parameters
        ----------
        a, b : np.ndarray
            One dimensional sequences.

        Returns
        -------
        i : int
            Splitting index according to above description.
        """
        if a.size != b.size:
            raise ValueError('a and b need to have the same size')
        if a.ndim != b.ndim != 1:
            raise ValueError('a and b need to be one-dimensional arrays')
        import scipy as sp
        min_length = 5
        n = a.size
        idx_range = np.arange(min_length, a.size-min_length-1, dtype=int)
        corr_coeff = np.zeros(idx_range.size)
        pos_old = sp.stats.kendalltau(a[:min_length], b[:min_length])[0]
        neg_old = sp.stats.kendalltau(a[min_length:], b[min_length:])[0]
        for ii, i in enumerate(idx_range):
            if True:
                # compute differences in concordance when adding a[i] and b[i]
                # to the first subsequence, and removing these elements from
                # the second subsequence
                diff_pos, diff_neg = self._kendall_tau_diff(a, b, i)
                pos = pos_old + self._kendall_tau_add(i, diff_pos, pos_old)
                neg = neg_old + self._kendall_tau_subtract(n-i, diff_neg, neg_old)
                pos_old = pos
                neg_old = neg
            if False:
                # computation using sp.stats.kendalltau, takes much longer!
                # just for debugging purposes
                pos = sp.stats.kendalltau(a[:i+1], b[:i+1])[0]
                neg = sp.stats.kendalltau(a[i+1:], b[i+1:])[0]
            if False:
                # the following is much slower than using sp.stats.kendalltau,
                # it is only good for debugging because it allows to compute the
                # tau-a version, which does not account for ties, whereas
                # sp.stats.kendalltau computes tau-b version, which accounts for
                # ties
                pos = sp.stats.mstats.kendalltau(a[:i], b[:i], use_ties=False)[0]
                neg = sp.stats.mstats.kendalltau(a[i:], b[i:], use_ties=False)[0]
            corr_coeff[ii] = pos - neg
        iimax = np.argmax(corr_coeff)
        imax = min_length + iimax
        corr_coeff_max = corr_coeff[iimax]
        if corr_coeff_max < 0.3:
            logg.m('... is root itself, never obtain significant correlation')
        return imax

    def _kendall_tau_add(self, len_old, diff_pos, tau_old):
        """Compute Kendall tau delta.

        The new sequence has length len_old + 1.

        Parameters
        ----------
        len_old : int
            The length of the old sequence, used to compute tau_old.
        diff_pos : int
            Difference between concordant and non-concordant pairs.
        tau_old : float
            Kendall rank correlation of the old sequence.
        """
        return 2./(len_old+1)*(float(diff_pos)/len_old-tau_old)

    def _kendall_tau_subtract(self, len_old, diff_neg, tau_old):
        """Compute Kendall tau delta.

        The new sequence has length len_old - 1.

        Parameters
        ----------
        len_old : int
            The length of the old sequence, used to compute tau_old.
        diff_neg : int
            Difference between concordant and non-concordant pairs.
        tau_old : float
            Kendall rank correlation of the old sequence.
        """
        return 2./(len_old-2)*(-float(diff_neg)/(len_old-1)+tau_old)

    def _kendall_tau_diff(self, a, b, i):
        """Compute difference in concordance of pairs in split sequences.

        Consider splitting a and b at index i.

        Parameters
        ----------
        a, b : np.ndarray

        Returns
        -------
        diff_pos, diff_neg : int, int
            Difference between concordant and non-concordant pairs for both
            subsequences.
        """
        # compute ordering relation of the single points a[i] and b[i]
        # with all previous points of the sequences a and b, respectively
        a_pos = np.zeros(a[:i].size, dtype=int)
        a_pos[a[:i] > a[i]] = 1
        a_pos[a[:i] < a[i]] = -1
        b_pos = np.zeros(b[:i].size, dtype=int)
        b_pos[b[:i] > b[i]] = 1
        b_pos[b[:i] < b[i]] = -1
        diff_pos = np.dot(a_pos, b_pos).astype(float)

        # compute ordering relation of the single points a[i] and b[i]
        # with all later points of the sequences
        a_neg = np.zeros(a[i:].size, dtype=int)
        a_neg[a[i:] > a[i]] = 1
        a_neg[a[i:] < a[i]] = -1
        b_neg = np.zeros(b[i:].size, dtype=int)
        b_neg[b[i:] > b[i]] = 1
        b_neg[b[i:] < b[i]] = -1
        diff_neg = np.dot(a_neg, b_neg)

        return diff_pos, diff_neg

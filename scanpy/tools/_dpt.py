from typing import Tuple

import numpy as np
import pandas as pd
import scipy as sp
import networkx as nx
from natsort import natsorted

from .. import logging as logg
from ..logging import _settings_verbosity_greater_or_equal_than
from ..neighbors import Neighbors, OnFlySymMatrix


def _diffmap(adata, n_comps=15):
    logg.info('computing Diffusion Maps using n_comps={}(=n_dcs)'.format(n_comps), r=True)
    dpt = DPT(adata)
    dpt.compute_transitions()
    dpt.compute_eigen(n_comps=n_comps)
    adata.obsm['X_diffmap'] = dpt.eigen_basis
    adata.uns['diffmap_evals'] = dpt.eigen_values
    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    logg.hint('added\n'
              '    \'X_diffmap\', diffmap coordinates (adata.obsm)\n'
              '    \'diffmap_evals\', eigenvalues of transition matrix (adata.uns)')


def dpt(adata, n_dcs=10, n_branchings=0, min_group_size=0.01,
        allow_kendall_tau_shift=True, copy=False):
    """Infer progression of cells through geodesic distance along the graph [Haghverdi16]_ [Wolf19]_.

    Reconstruct the progression of a biological process from snapshot
    data. `Diffusion Pseudotime` has been introduced by [Haghverdi16]_ and
    implemented within Scanpy [Wolf18]_. Here, we use a further developed
    version, which is able to deal with disconnected graphs [Wolf19]_ and can
    be run in a `hierarchical` mode by setting the parameter
    `n_branchings>1`. We recommend, however, to only use
    :func:`~scanpy.api.tl.dpt` for computing pseudotime (`n_branchings=0`) and
    to detect branchings via :func:`~scanpy.api.paga`. For pseudotime, you need
    to annotate your data with a root cell. For instance::

        adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_types'] == 'Stem')[0]

    This requires to run :func:`~scanpy.api.pp.neighbors`, first. In order to
    reproduce the original implementation of DPT, use `method=='gauss'` in
    this. Using the default `method=='umap'` only leads to minor quantitative
    differences, though.

    .. versionadded:: 1.1

    :func:`~scanpy.api.tl.dpt` also requires to run
    :func:`~scanpy.api.tl.diffmap` first. As previously,
    :func:`~scanpy.api.tl.dpt` came with a default parameter of ``n_dcs=10`` but
    :func:`~scanpy.api.tl.diffmap` has a default parameter of ``n_comps=15``,
    you need to pass ``n_comps=10`` in :func:`~scanpy.api.tl.diffmap` in order
    to exactly reproduce previous :func:`~scanpy.api.tl.dpt` results.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_dcs : `int`, optional (default: 10)
        The number of diffusion components to use.
    n_branchings : `int`, optional (default: 0)
        Number of branchings to detect.
    min_group_size : [0, 1] or `float`, optional (default: 0.01)
        During recursive splitting of branches ('dpt groups') for `n_branchings`
        > 1, do not consider groups that contain less than `min_group_size` data
        points. If a float, `min_group_size` refers to a fraction of the total
        number of data points.
    allow_kendall_tau_shift : `bool`, optional (default: `True`)
        If a very small branch is detected upon splitting, shift away from
        maximum correlation in Kendall tau criterion of [Haghverdi16]_ to
        stabilize the splitting.
    copy : `bool`, optional (default: `False`)
        Copy instance before computation and return a copy. Otherwise, perform
        computation inplace and return None.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    If `n_branchings==0`, no field `dpt_groups` will be written.

    **dpt_pseudotime** : :class:`pandas.Series` (`adata.obs`, dtype `float`)
        Array of dim (number of samples) that stores the pseudotime of each
        cell, that is, the DPT distance with respect to the root cell.
    **dpt_groups** : :class:`pandas.Series` (`adata.obs`, dtype `category`)
        Array of dim (number of samples) that stores the subgroup id ('0',
        '1', ...) for each cell. The groups  typically correspond to
        'progenitor cells', 'undecided cells' or 'branches' of a process.

    Notes
    -----
    The tool is similar to the R package `destiny` of [Angerer16]_.
    """
    # standard errors, warnings etc.
    adata = adata.copy() if copy else adata
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` and `tl.diffmap` first.')
    if 'iroot' not in adata.uns and 'xroot' not in adata.var:
        logg.warn(
            'No root cell found. To compute pseudotime, pass the index or '
            'expression vector of a root cell, one of:\n'
            '    adata.uns[\'iroot\'] = root_cell_index\n'
            '    adata.var[\'xroot\'] = adata[root_cell_name, :].X')
    if 'X_diffmap' not in adata.obsm.keys():
        logg.warn('Trying to run `tl.dpt` without prior call of `tl.diffmap`. '
                  'Falling back to `tl.diffmap` with default parameters.')
        _diffmap(adata)
    # start with the actual computation
    dpt = DPT(adata, n_dcs=n_dcs, min_group_size=min_group_size,
              n_branchings=n_branchings,
              allow_kendall_tau_shift=allow_kendall_tau_shift)
    logg.info('computing Diffusion Pseudotime using n_dcs={}'.format(n_dcs), r=True)
    if n_branchings > 1: logg.info('    this uses a hierarchical implementation')
    if dpt.iroot is not None:
        dpt._set_pseudotime()  # pseudotimes are distances from root point
        adata.uns['iroot'] = dpt.iroot  # update iroot, might have changed when subsampling, for example
        adata.obs['dpt_pseudotime'] = dpt.pseudotime
    # detect branchings and partition the data into segments
    if n_branchings > 0:
        dpt.branchings_segments()
        adata.obs['dpt_groups'] = pd.Categorical(
            values=dpt.segs_names.astype('U'),
            categories=natsorted(np.array(dpt.segs_names_unique).astype('U')))
        # the "change points" separate segments in the ordering above
        adata.uns['dpt_changepoints'] = dpt.changepoints
        # the tip points of segments
        adata.uns['dpt_grouptips'] = dpt.segs_tips
        # the ordering according to segments and pseudotime
        ordering_id = np.zeros(adata.n_obs, dtype=int)
        for count, idx in enumerate(dpt.indices): ordering_id[idx] = count
        adata.obs['dpt_order'] = ordering_id
        adata.obs['dpt_order_indices'] = dpt.indices
    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    logg.hint('added\n'
           + ('    \'dpt_pseudotime\', the pseudotime (adata.obs)'
              if dpt.iroot is not None else '')
           + ('\n    \'dpt_groups\', the branching subgroups of dpt (adata.obs)\n'
              + '    \'dpt_order\', cell order (adata.obs)'
              if n_branchings > 0 else ''))
    return adata if copy else None


class DPT(Neighbors):
    """Hierarchical Diffusion Pseudotime.
    """

    def __init__(self, adata, n_dcs=None, min_group_size=0.01,
                 n_branchings=0, allow_kendall_tau_shift=False):
        super(DPT, self).__init__(adata, n_dcs=n_dcs)
        self.flavor = 'haghverdi16'
        self.n_branchings = n_branchings
        self.min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * self._adata.shape[0])
        self.passed_adata = adata  # just for debugging purposes
        self.choose_largest_segment = False
        self.allow_kendall_tau_shift = allow_kendall_tau_shift

    def branchings_segments(self):
        """Detect branchings and partition the data into corresponding segments.

        Detect all branchings up to `n_branchings`.

        Writes
        ------
        segs : np.ndarray
            Array of dimension (number of segments) × (number of data
            points). Each row stores a mask array that defines a segment.
        segs_tips : np.ndarray
            Array of dimension (number of segments) × 2. Each row stores the
            indices of the two tip points of each segment.
        segs_names : np.ndarray
            Array of dimension (number of data points). Stores an integer label
            for each segment.
        """
        self.detect_branchings()
        self.postprocess_segments()
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
        logg.m('    detect', self.n_branchings,
               'branching' + ('' if self.n_branchings == 1 else 's'))
        # a segment is a subset of points of the data set (defined by the
        # indices of the points in the segment)
        # initialize the search for branchings with a single segment,
        # that is, get the indices of the whole data set
        indices_all = np.arange(self._adata.shape[0], dtype=int)
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
            tips_all = np.array(np.unravel_index(np.argmax(self.distances_dpt), self.distances_dpt.shape))
        else:
            if self.iroot is not None:
                tip_0 = np.argmax(self.distances_dpt[self.iroot])
            else:
                tip_0 = np.argmax(self.distances_dpt[0])
            tips_all = np.array([tip_0, np.argmax(self.distances_dpt[tip_0])])
        # we keep a list of the tips of each segment
        segs_tips = [tips_all]
        segs_connects = [[]]
        segs_undecided = [True]
        segs_adjacency = [[]]
        logg.m('    do not consider groups with less than {} points for splitting'
               .format(self.min_group_size))
        for ibranch in range(self.n_branchings):
            iseg, tips3 = self.select_segment(segs, segs_tips, segs_undecided)
            if iseg == -1:
                logg.m('    partitioning converged')
                break
            logg.m('    branching {}:'.format(ibranch + 1),
                   'split group', iseg)  # [third start end]
            # detect branching and update segs and segs_tips
            self.detect_branching(segs, segs_tips,
                                  segs_connects,
                                  segs_undecided,
                                  segs_adjacency, iseg, tips3)
        # store as class members
        self.segs = segs
        self.segs_tips = segs_tips
        self.segs_undecided = segs_undecided
        # the following is a bit too much, but this allows easy storage
        self.segs_adjacency = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
        self.segs_connects = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=int)
        for i, seg_adjacency in enumerate(segs_adjacency):
            self.segs_connects[i, seg_adjacency] = segs_connects[i]
        for i in range(len(segs)):
            for j in range(len(segs)):
                self.segs_adjacency[i, j] = self.distances_dpt[self.segs_connects[i, j],
                                                         self.segs_connects[j, i]]
        self.segs_adjacency = self.segs_adjacency.tocsr()
        self.segs_connects = self.segs_connects.tocsr()

    def check_adjacency(self):
        n_edges_per_seg = np.sum(self.segs_adjacency > 0, axis=1).A1
        for n_edges in range(1, np.max(n_edges_per_seg) + 1):
            for iseg in range(self.segs_adjacency.shape[0]):
                if n_edges_per_seg[iseg] == n_edges:
                    neighbor_segs = self.segs_adjacency[iseg].todense().A1
                    closest_points_other_segs = [seg[np.argmin(self.distances_dpt[self.segs_tips[iseg][0], seg])]
                                                 for seg in self.segs]
                    seg = self.segs[iseg]
                    closest_points_in_segs = [seg[np.argmin(self.distances_dpt[tips[0], seg])]
                                              for tips in self.segs_tips]
                    distance_segs = [self.distances_dpt[closest_points_other_segs[ipoint], point]
                                     for ipoint, point in enumerate(closest_points_in_segs)]
                    # exclude the first point, the segment itself
                    closest_segs = np.argsort(distance_segs)[1:n_edges+1]
                    # update adjacency matrix within the loop!
                    # self.segs_adjacency[iseg, neighbor_segs > 0] = 0
                    # self.segs_adjacency[iseg, closest_segs] = np.array(distance_segs)[closest_segs]
                    # self.segs_adjacency[neighbor_segs > 0, iseg] = 0
                    # self.segs_adjacency[closest_segs, iseg] = np.array(distance_segs)[closest_segs].reshape(len(closest_segs), 1)
                    # n_edges_per_seg = np.sum(self.segs_adjacency > 0, axis=1).A1
                    print(iseg, distance_segs, closest_segs)
                    # print(self.segs_adjacency)
        # self.segs_adjacency.eliminate_zeros()

    def select_segment(self, segs, segs_tips, segs_undecided) -> Tuple[int, int]:
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
        allindices = np.arange(self._adata.shape[0], dtype=int)
        for iseg, seg in enumerate(segs):
            # do not consider too small segments
            if segs_tips[iseg][0] == -1: continue
            # restrict distance matrix to points in segment
            if not isinstance(self.distances_dpt, OnFlySymMatrix):
                Dseg = self.distances_dpt[np.ix_(seg, seg)]
            else:
                Dseg = self.distances_dpt.restrict(seg)
            third_maximizer = None
            if segs_undecided[iseg]:
                # check that none of our tips "connects" with a tip of the
                # other segments
                for jseg in range(len(segs)):
                    if jseg != iseg:
                        # take the inner tip, the "second tip" of the segment
                        for itip in range(2):
                            if (self.distances_dpt[segs_tips[jseg][1], segs_tips[iseg][itip]]
                                < 0.5 * self.distances_dpt[segs_tips[iseg][~itip], segs_tips[iseg][itip]]):
                                # logg.m('    group', iseg, 'with tip', segs_tips[iseg][itip],
                                #        'connects with', jseg, 'with tip', segs_tips[jseg][1], v=4)
                                # logg.m('    do not use the tip for "triangulation"', v=4)
                                third_maximizer = itip
            # map the global position to the position within the segment
            tips = [np.where(allindices[seg] == tip)[0][0]
                    for tip in segs_tips[iseg]]
            # find the third point on the segment that has maximal
            # added distance from the two tip points
            dseg = Dseg[tips[0]] + Dseg[tips[1]]
            if not np.isfinite(dseg).any():
                continue
            # add this point to tips, it's a third tip, we store it at the first
            # position in an array called tips3
            third_tip = np.argmax(dseg)
            if third_maximizer is not None:
                # find a fourth point that has maximal distance to all three
                dseg += Dseg[third_tip]
                fourth_tip = np.argmax(dseg)
                if fourth_tip != tips[0] and fourth_tip != third_tip:
                    tips[1] = fourth_tip
                    dseg -= Dseg[tips[1]]
                else:
                    dseg -= Dseg[third_tip]
            tips3 = np.append(tips, third_tip)
            # compute the score as ratio of the added distance to the third tip,
            # to what it would be if it were on the straight line between the
            # two first tips, given by Dseg[tips[:2]]
            # if we did not normalize, there would be a danger of simply
            # assigning the highest score to the longest segment
            score = dseg[tips3[2]] / Dseg[tips3[0], tips3[1]]
            score = len(seg) if self.choose_largest_segment else score  # simply the number of points
            logg.m('    group', iseg, 'score', score, 'n_points', len(seg),
                   '(too small)' if len(seg) < self.min_group_size else '', v=4)
            if len(seg) <= self.min_group_size: score = 0
            # write result
            scores_tips[iseg, 0] = score
            scores_tips[iseg, 1:] = tips3
        iseg = np.argmax(scores_tips[:, 0])
        if scores_tips[iseg, 0] == 0: return -1, None
        tips3 = scores_tips[iseg, 1:].astype(int)
        return iseg, tips3

    def postprocess_segments(self):
        """Convert the format of the segment class members."""
        # make segs a list of mask arrays, it's easier to store
        # as there is a hdf5 equivalent
        for iseg, seg in enumerate(self.segs):
            mask = np.zeros(self._adata.shape[0], dtype=bool)
            mask[seg] = True
            self.segs[iseg] = mask
        # convert to arrays
        self.segs = np.array(self.segs)
        self.segs_tips = np.array(self.segs_tips)

    def set_segs_names(self):
        """Return a single array that stores integer segment labels."""
        segs_names = np.zeros(self._adata.shape[0], dtype=np.int8)
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
        if self.iroot is not None:
            for itips, tips in enumerate(self.segs_tips):
                if tips[0] != -1:
                    indices = np.argsort(self.pseudotime[tips])
                    self.segs_tips[itips] = self.segs_tips[itips][indices]
                else:
                    logg.m('    group', itips, 'is very small', v=4)
        # sort indices according to segments
        indices = np.argsort(self.segs_names)
        segs_names = self.segs_names[indices]
        # find changepoints of segments
        changepoints = np.arange(indices.size-1)[np.diff(segs_names) == 1] + 1
        if self.iroot is not None:
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

    def detect_branching(self, segs, segs_tips, segs_connects, segs_undecided, segs_adjacency,
                         iseg, tips3):
        """Detect branching on given segment.

        Updates all list parameters inplace.

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
        """
        seg = segs[iseg]
        # restrict distance matrix to points in segment
        if not isinstance(self.distances_dpt, OnFlySymMatrix):
            Dseg = self.distances_dpt[np.ix_(seg, seg)]
        else:
            Dseg = self.distances_dpt.restrict(seg)
        # given the three tip points and the distance matrix detect the
        # branching on the segment, return the list ssegs of segments that
        # are defined by splitting this segment
        result = self._detect_branching(Dseg, tips3, seg)
        ssegs, ssegs_tips, ssegs_adjacency, ssegs_connects, trunk = result
        # map back to global indices
        for iseg_new, seg_new in enumerate(ssegs):
            ssegs[iseg_new] = seg[seg_new]
            ssegs_tips[iseg_new] = seg[ssegs_tips[iseg_new]]
            ssegs_connects[iseg_new] = list(seg[ssegs_connects[iseg_new]])
        # remove previous segment
        segs.pop(iseg)
        segs_tips.pop(iseg)
        # insert trunk/undecided_cells at same position
        segs.insert(iseg, ssegs[trunk])
        segs_tips.insert(iseg, ssegs_tips[trunk])
        # append other segments
        segs += [seg for iseg, seg in enumerate(ssegs) if iseg != trunk]
        segs_tips += [seg_tips for iseg, seg_tips in enumerate(ssegs_tips) if iseg != trunk]
        if len(ssegs) == 4:
            # insert undecided cells at same position
            segs_undecided.pop(iseg)
            segs_undecided.insert(iseg, True)
        # correct edges in adjacency matrix
        n_add = len(ssegs) - 1
        prev_connecting_segments = segs_adjacency[iseg].copy()
        if self.flavor == 'haghverdi16':
            segs_adjacency += [[iseg] for i in range(n_add)]
            segs_connects += [seg_connects for iseg, seg_connects in enumerate(ssegs_connects) if iseg != trunk]
            prev_connecting_points = segs_connects[iseg]
            for jseg_cnt, jseg in enumerate(prev_connecting_segments):
                iseg_cnt = 0
                for iseg_new, seg_new in enumerate(ssegs):
                    if iseg_new != trunk:
                        pos = segs_adjacency[jseg].index(iseg)
                        connection_to_iseg = segs_connects[jseg][pos]
                        if connection_to_iseg in seg_new:
                            kseg = len(segs) - n_add + iseg_cnt
                            segs_adjacency[jseg][pos] = kseg
                            pos_2 = segs_adjacency[iseg].index(jseg)
                            segs_adjacency[iseg].pop(pos_2)
                            idx = segs_connects[iseg].pop(pos_2)
                            segs_adjacency[kseg].append(jseg)
                            segs_connects[kseg].append(idx)
                            break
                        iseg_cnt += 1
            segs_adjacency[iseg] += list(range(len(segs_adjacency) - n_add, len(segs_adjacency)))
            segs_connects[iseg] += ssegs_connects[trunk]
        else:
            segs_adjacency += [[] for i in range(n_add)]
            segs_connects += [[] for i in range(n_add)]
            kseg_list = [iseg] + list(range(len(segs) - n_add, len(segs)))
            for jseg in prev_connecting_segments:
                pos = segs_adjacency[jseg].index(iseg)
                distances = []
                closest_points_in_jseg = []
                closest_points_in_kseg = []
                for kseg in kseg_list:
                    reference_point_in_k = segs_tips[kseg][0]
                    closest_points_in_jseg.append(segs[jseg][np.argmin(self.distances_dpt[reference_point_in_k, segs[jseg]])])
                    # do not use the tip in the large segment j, instead, use the closest point
                    reference_point_in_j = closest_points_in_jseg[-1]  # segs_tips[jseg][0]
                    closest_points_in_kseg.append(segs[kseg][np.argmin(self.distances_dpt[reference_point_in_j, segs[kseg]])])
                    distances.append(self.distances_dpt[closest_points_in_jseg[-1], closest_points_in_kseg[-1]])
                    # print(jseg, '(', segs_tips[jseg][0], closest_points_in_jseg[-1], ')',
                    #       kseg, '(', segs_tips[kseg][0], closest_points_in_kseg[-1], ') :', distances[-1])
                idx = np.argmin(distances)
                kseg_min = kseg_list[idx]
                segs_adjacency[jseg][pos] = kseg_min
                segs_connects[jseg][pos] = closest_points_in_kseg[idx]
                pos_2 = segs_adjacency[iseg].index(jseg)
                segs_adjacency[iseg].pop(pos_2)
                segs_connects[iseg].pop(pos_2)
                segs_adjacency[kseg_min].append(jseg)
                segs_connects[kseg_min].append(closest_points_in_jseg[idx])
            # if we split two clusters, we need to check whether the new segments connect to any of the other
            # old segments
            # if not, we add a link between the new segments, if yes, we add two links to connect them at the
            # correct old segments
            do_not_attach_kseg = False
            for kseg in kseg_list:
                distances = []
                closest_points_in_jseg = []
                closest_points_in_kseg = []
                jseg_list = [jseg for jseg in range(len(segs))
                             if jseg != kseg and jseg not in prev_connecting_segments]
                for jseg in jseg_list:
                    reference_point_in_k = segs_tips[kseg][0]
                    closest_points_in_jseg.append(segs[jseg][np.argmin(self.distances_dpt[reference_point_in_k, segs[jseg]])])
                    # do not use the tip in the large segment j, instead, use the closest point
                    reference_point_in_j = closest_points_in_jseg[-1]  # segs_tips[jseg][0]
                    closest_points_in_kseg.append(segs[kseg][np.argmin(self.distances_dpt[reference_point_in_j, segs[kseg]])])
                    distances.append(self.distances_dpt[closest_points_in_jseg[-1], closest_points_in_kseg[-1]])
                idx = np.argmin(distances)
                jseg_min = jseg_list[idx]
                if jseg_min not in kseg_list:
                    segs_adjacency_sparse = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
                    for i, seg_adjacency in enumerate(segs_adjacency):
                        segs_adjacency_sparse[i, seg_adjacency] = 1
                    G = nx.Graph(segs_adjacency_sparse)
                    paths_all = nx.single_source_dijkstra_path(G, source=kseg)
                    if jseg_min not in paths_all:
                        segs_adjacency[jseg_min].append(kseg)
                        segs_connects[jseg_min].append(closest_points_in_kseg[idx])
                        segs_adjacency[kseg].append(jseg_min)
                        segs_connects[kseg].append(closest_points_in_jseg[idx])
                        logg.m('    attaching new segment', kseg, 'at', jseg_min)
                        # if we split the cluster, we should not attach kseg
                        do_not_attach_kseg = True
                    else:
                        logg.m('    cannot attach new segment', kseg, 'at', jseg_min,
                               '(would produce cycle)')
                        if kseg != kseg_list[-1]:
                            logg.m('        continue')
                            continue
                        else:
                            logg.m('        do not add another link')
                            break
                if jseg_min in kseg_list and not do_not_attach_kseg:
                    segs_adjacency[jseg_min].append(kseg)
                    segs_connects[jseg_min].append(closest_points_in_kseg[idx])
                    segs_adjacency[kseg].append(jseg_min)
                    segs_connects[kseg].append(closest_points_in_jseg[idx])
                    break
        segs_undecided += [False for i in range(n_add)]

    def _detect_branching(self, Dseg: np.ndarray, tips: np.ndarray, seg_reference=None):
        """Detect branching on given segment.

        Call function __detect_branching three times for all three orderings of
        tips. Points that do not belong to the same segment in all three
        orderings are assigned to a fourth segment. The latter is, by Haghverdi
        et al. (2016) referred to as 'undecided cells'.

        Parameters
        ----------
        Dseg
            Dchosen distance matrix restricted to segment.
        tips
            The three tip points. They form a 'triangle' that contains the data.

        Returns
        -------
        ssegs : list of np.ndarray
            List of segments obtained from splitting the single segment defined
            via the first two tip cells.
        ssegs_tips : list of np.ndarray
            List of tips of segments in ssegs.
        """
        if self.flavor == 'haghverdi16':
            ssegs = self._detect_branching_single_haghverdi16(Dseg, tips)
        elif self.flavor == 'wolf17_tri':
            ssegs = self._detect_branching_single_wolf17_tri(Dseg, tips)
        elif self.flavor == 'wolf17_bi' or self.flavor == 'wolf17_bi_un':
            ssegs = self._detect_branching_single_wolf17_bi(Dseg, tips)
        else:
            raise ValueError('`flavor` needs to be in {"haghverdi16", "wolf17_tri", "wolf17_bi"}.')
        # make sure that each data point has a unique association with a segment
        masks = np.zeros((len(ssegs), Dseg.shape[0]), dtype=bool)
        for iseg, seg in enumerate(ssegs):
            masks[iseg][seg] = True
        nonunique = np.sum(masks, axis=0) > 1
        ssegs = []
        for iseg, mask in enumerate(masks):
            mask[nonunique] = False
            ssegs.append(np.arange(Dseg.shape[0], dtype=int)[mask])
        # compute new tips within new segments
        ssegs_tips = []
        for inewseg, newseg in enumerate(ssegs):
            if len(np.flatnonzero(newseg)) <= 1:
                logg.warn('detected group with only {} cells'.format(np.flatnonzero(newseg)))
            secondtip = newseg[np.argmax(Dseg[tips[inewseg]][newseg])]
            ssegs_tips.append([tips[inewseg], secondtip])
        undecided_cells = np.arange(Dseg.shape[0], dtype=int)[nonunique]
        if len(undecided_cells) > 0:
            ssegs.append(undecided_cells)
            # establish the connecting points with the other segments
            ssegs_connects = [[], [], [], []]
            for inewseg, newseg_tips in enumerate(ssegs_tips):
                reference_point = newseg_tips[0]
                # closest cell to the new segment within undecided cells
                closest_cell = undecided_cells[np.argmin(Dseg[reference_point][undecided_cells])]
                ssegs_connects[inewseg].append(closest_cell)
                # closest cell to the undecided cells within new segment
                closest_cell = ssegs[inewseg][np.argmin(Dseg[closest_cell][ssegs[inewseg]])]
                ssegs_connects[-1].append(closest_cell)
            # also compute tips for the undecided cells
            tip_0 = undecided_cells[np.argmax(Dseg[undecided_cells[0]][undecided_cells])]
            tip_1 = undecided_cells[np.argmax(Dseg[tip_0][undecided_cells])]
            ssegs_tips.append([tip_0, tip_1])
            ssegs_adjacency = [[3], [3], [3], [0, 1, 2]]
            trunk = 3
        elif len(ssegs) == 3:
            reference_point = np.zeros(3, dtype=int)
            reference_point[0] = ssegs_tips[0][0]
            reference_point[1] = ssegs_tips[1][0]
            reference_point[2] = ssegs_tips[2][0]
            closest_points = np.zeros((3, 3), dtype=int)
            # this is another strategy than for the undecided_cells
            # here it's possible to use the more symmetric procedure
            # shouldn't make much of a difference
            closest_points[0, 1] = ssegs[1][np.argmin(Dseg[reference_point[0]][ssegs[1]])]
            closest_points[1, 0] = ssegs[0][np.argmin(Dseg[reference_point[1]][ssegs[0]])]
            closest_points[0, 2] = ssegs[2][np.argmin(Dseg[reference_point[0]][ssegs[2]])]
            closest_points[2, 0] = ssegs[0][np.argmin(Dseg[reference_point[2]][ssegs[0]])]
            closest_points[1, 2] = ssegs[2][np.argmin(Dseg[reference_point[1]][ssegs[2]])]
            closest_points[2, 1] = ssegs[1][np.argmin(Dseg[reference_point[2]][ssegs[1]])]
            added_dist = np.zeros(3)
            added_dist[0] = Dseg[closest_points[1, 0], closest_points[0, 1]] + Dseg[closest_points[2, 0], closest_points[0, 2]]
            added_dist[1] = Dseg[closest_points[0, 1], closest_points[1, 0]] + Dseg[closest_points[2, 1], closest_points[1, 2]]
            added_dist[2] = Dseg[closest_points[1, 2], closest_points[2, 1]] + Dseg[closest_points[0, 2], closest_points[2, 0]]
            trunk = np.argmin(added_dist)
            ssegs_adjacency = [[trunk] if i != trunk else
                               [j for j in range(3) if j != trunk]
                               for i in range(3)]
            ssegs_connects = [[closest_points[i, trunk]] if i != trunk else
                              [closest_points[trunk, j] for j in range(3) if j != trunk]
                              for i in range(3)]
        else:
            trunk = 0
            ssegs_adjacency = [[1], [0]]
            reference_point_in_0 = ssegs_tips[0][0]
            closest_point_in_1 = ssegs[1][np.argmin(Dseg[reference_point_in_0][ssegs[1]])]
            reference_point_in_1 = closest_point_in_1  # ssegs_tips[1][0]
            closest_point_in_0 = ssegs[0][np.argmin(Dseg[reference_point_in_1][ssegs[0]])]
            ssegs_connects = [[closest_point_in_1], [closest_point_in_0]]
        return ssegs, ssegs_tips, ssegs_adjacency, ssegs_connects, trunk

    def _detect_branching_single_haghverdi16(self, Dseg, tips):
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
            ssegs.append(self.__detect_branching_haghverdi16(Dseg, tips[p]))
        return ssegs

    def _detect_branching_single_wolf17_tri(self, Dseg, tips):
        # all pairwise distances
        dist_from_0 = Dseg[tips[0]]
        dist_from_1 = Dseg[tips[1]]
        dist_from_2 = Dseg[tips[2]]
        closer_to_0_than_to_1 = dist_from_0 < dist_from_1
        closer_to_0_than_to_2 = dist_from_0 < dist_from_2
        closer_to_1_than_to_2 = dist_from_1 < dist_from_2
        masks = np.zeros((2, Dseg.shape[0]), dtype=bool)
        masks[0] = closer_to_0_than_to_1
        masks[1] = closer_to_0_than_to_2
        segment_0 = np.sum(masks, axis=0) == 2
        masks = np.zeros((2, Dseg.shape[0]), dtype=bool)
        masks[0] = ~closer_to_0_than_to_1
        masks[1] = closer_to_1_than_to_2
        segment_1 = np.sum(masks, axis=0) == 2
        masks = np.zeros((2, Dseg.shape[0]), dtype=bool)
        masks[0] = ~closer_to_0_than_to_2
        masks[1] = ~closer_to_1_than_to_2
        segment_2 = np.sum(masks, axis=0) == 2
        ssegs = [segment_0, segment_1, segment_2]
        return ssegs

    def _detect_branching_single_wolf17_bi(self, Dseg, tips):
        dist_from_0 = Dseg[tips[0]]
        dist_from_1 = Dseg[tips[1]]
        closer_to_0_than_to_1 = dist_from_0 < dist_from_1
        ssegs = [closer_to_0_than_to_1, ~closer_to_0_than_to_1]
        return ssegs

    def __detect_branching_haghverdi16(self, Dseg, tips):
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
        if imax > 0.95 * len(idcs) and self.allow_kendall_tau_shift:
            # if "everything" is correlated (very large value of imax), a more
            # conservative choice amounts to reducing this
            logg.warn('shifting branching point away from maximal kendall-tau correlation (suppress this with `allow_kendall_tau_shift=False`)')
            ibranch = int(0.95 * imax)
        else:
            # otherwise, a more conservative choice is the following
            ibranch = imax + 1
        return idcs[:ibranch]

    def kendall_tau_split(self, a, b) -> int:
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
            logg.m('    is root itself, never obtain significant correlation', v=4)
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

    def _kendall_tau_diff(self, a: np.ndarray, b: np.ndarray, i) -> Tuple[int, int]:
        """Compute difference in concordance of pairs in split sequences.

        Consider splitting a and b at index i.

        Parameters
        ----------
        a
            ?
        b
            ?

        Returns
        -------
        diff_pos
            Difference between concordant pairs for both subsequences.
        diff_neg
            Difference between non-concordant pairs for both subsequences.
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

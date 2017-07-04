# Author: F. Alex Wolf (http://falexwolf.de)
"""Extremal Graph Abstraction
"""

import sys
import numpy as np
import scipy as sp
import networkx as nx
import scipy.sparse
from .. import logging as logg
from ..data_structs import data_graph


def ega(adata, n_splits=0, k=30, knn=True, n_pcs=50, n_dcs=10,
        min_group_size=0.01, n_jobs=None, recompute_diffmap=False,
        recompute_pca=False, flavor='bi_constrained',
        attachedness_measure='n_connecting_edges',
        copy=False):
    """Extremal Graph Abstraction

    Infer an abstraction of the relations of subgroups in the data through
    extremal graph abstraction. The result is a hierarchy of cell subgroups
    reprsented as a tree on which each node corresponds to a cell subgroup. The
    tree induces an ordering between nodes and the cells are ordered within each
    node.

    Reference
    ---------
    Wolf et al., bioRxiv (2017)

    See also
    --------
    - Diffusion Pseudotime: Haghverdi et al., Nature Methods 13, 3971 (2016).

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
    n_splits : int, optional (default: 1)
        Number of splits to detect.
    k : int, optional (default: 30)
        Number of nearest neighbors on the knn graph.
    n_pcs: int, optional (default: 50)
        Use n_pcs PCs to compute the Euclidian distance matrix, which is the
        basis for generating the graph. Set to 0 if you don't want preprocessing
        with PCA.
    n_jobs : int or None (default: None)
        Number of cpus to use for parallel processing (default: sett.n_jobs).
    recompute_diffmap : bool, (default: False)
        Recompute diffusion maps.
    recompute_pca : bool, (default: False)
        Recompute PCA.
    attachedness_measure : {'n_connecting_edges', 'd_full_pairwise'} (default: 'n_connecting_edges')
        How to measure attachedness. Parameter only for developing.
    copy : bool, optional (default: False)
        Copy instance before computation and return a copy. Otherwise, perform
        computation inplace and return None.

    Notes
    -----
    Writes the following arrays as sample annotation to adata.smp.
        ega_adjacency : sparse csr matrix
            Array of dim (number of samples) that stores the pseudotime of each
            cell, that is, the DPT distance with respect to the root cell.
        ega_groups : np.ndarray of dtype string
            Array of dim (number of samples) that stores the subgroup id ('0',
            '1', ...) for each cell.
    """
    adata = adata.copy() if copy else adata
    root_cell_was_passed = True
    if 'xroot' not in adata.add and 'xroot' not in adata.var:
        root_cell_was_passed = False
        logg.m('... no root cell found, no computation of pseudotime')
        msg = \
    '''To enable computation of pseudotime, pass the expression "xroot" of a root cell.
    Either add
        adata.var['xroot'] = adata.X[root_cell_index, :]
    where `root_cell_index` is the integer index of the root cell, or
        adata.var['xroot'] = adata[root_cell_name, :].X
    where `root_cell_name` is the name (a string) of the root cell.'''
        logg.hint(msg)
    if n_splits == 0:
        logg.hint('set parameter `n_splits` > 0 to detect splits')
    if n_splits > 1:
        logg.info('... running extremal graph abstraction (EGA)')
    ega = EGA(adata, k=k, n_pcs=n_pcs,
              min_group_size=min_group_size,
              n_jobs=n_jobs,
              recompute_diffmap=recompute_diffmap,
              recompute_pca=recompute_pca,
              n_splits=n_splits,
              attachedness_measure=attachedness_measure,
              flavor=flavor)
    ddmap = ega.diffmap(n_comps=n_dcs)
    adata.smp['X_diffmap'] = ddmap['X_diffmap']
    # also store the 0th comp, which is skipped for plotting
    adata.smp['X_diffmap0'] = ega.rbasis[:, 0]
    adata.add['diffmap_evals'] = ddmap['evals']
    adata.add['distance'] = ega.Dsq
    logg.info('perform Diffusion Pseudotime analysis', r=True)
    if False:
        # compute M matrix of cumulative transition probabilities,
        # see Haghverdi et al. (2016)
        ega.compute_M_matrix()
    # compute EGA distance matrix, which we refer to as 'Ddiff'
    if False:  # we do not compute the full Ddiff matrix, only the elements we need
        ega.compute_Ddiff_matrix()
    if root_cell_was_passed:
        ega.set_pseudotime()  # pseudotimes are distances from root point
        adata.add['iroot'] = ega.iroot  # update iroot, might have changed when subsampling, for example
        adata.smp['ega_pseudotime'] = ega.pseudotime
    # detect splits and partition the data into segments
    ega.splits_segments()
    # vector of length n_groups
    adata.add['ega_groups_names'] = [str(n) for n in ega.segs_names_unique]
    # for itips, tips in enumerate(ega.segs_tips):
    #     # if tips[0] == -1: adata.add['ega_groups_names'][itips] = '?'
    #     if ega.segs_undecided[itips]: adata.add['ega_groups_names'][itips] += '?'
    # vector of length n_samples of groupnames
    adata.smp['ega_groups'] = ega.segs_names.astype('U')
    # the ordering according to segments and pseudotime
    adata.smp['ega_order'] = ega.indices
    # the changepoints - marking different segments - in the ordering above
    adata.add['ega_changepoints'] = ega.changepoints
    # the tip points of segments
    adata.add['ega_grouptips'] = ega.segs_tips
    # the tree/graph adjacency matrix
    adata.add['ega_groups_adjacency'] = ega.segs_adjacency
    logg.info('finished', t=True, end=' ')
    logg.info('and added\n'
           + ('    "ega_pseudotime", stores pseudotime (adata.smp),\n' if root_cell_was_passed else '')
           + '    "ega_groups", the segments of trajectories a long a tree (adata.smp),\n'
           '    "ega_groups_adjacency", the adjacency matrix defining the tree (adata.add),\n'
           '    "ega_order", an index array for sorting the cells (adata.smp),\n'
           '    "ega_grouptips", the indices of tip cells (adata.add)')
    return adata if copy else None


class EGA(data_graph.DataGraph):
    """Extremal Graph Abstraction
    """

    def __init__(self, adata_or_X, k=30,
                 n_jobs=1, n_pcs=50,
                 min_group_size=20,
                 recompute_pca=None,
                 recompute_diffmap=None, n_splits=0,
                 attachedness_measure='n_connecting_edges',
                 flavor='haghverdi16'):
        super(EGA, self).__init__(adata_or_X, k=k, n_pcs=n_pcs,
                                  n_jobs=n_jobs,
                                  recompute_pca=recompute_pca,
                                  recompute_diffmap=recompute_diffmap,
                                  flavor=flavor)
        self.n_splits = n_splits
        self.min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * self.X.shape[0])
        self.passed_adata = adata_or_X  # just for debugging purposes
        self.choose_largest_segment = True
        self.attachedness_measure = attachedness_measure

    def splits_segments(self):
        """Detect splits and partition the data into corresponding segments.

        Detect all splits up to `n_splits`.

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
        self.detect_splits()
        self.postprocess_segments()
        self.set_segs_names()
        self.order_pseudotime()

    def detect_splits(self):
        """Detect all splits up to `n_splits`.

        Writes Attributes
        -----------------
        segs : np.ndarray
            List of integer index arrays.
        segs_tips : np.ndarray
            List of indices of the tips of segments.
        """
        logg.info('... detect', self.n_splits,
                  'branching' + ('' if self.n_splits == 1 else 's'))
        indices_all = np.arange(self.X.shape[0], dtype=int)
        segs = [indices_all]
        if False:  # this is safe, but not compatible with on-the-fly computation
            tips_all = np.array(np.unravel_index(np.argmax(self.Dchosen), self.Dchosen.shape))
        else:
            if self.iroot is not None:
                tip_0 = np.argmax(self.Dchosen[self.iroot])
            else:
                tip_0 = np.argmax(self.Dchosen[0])
            tips_all = np.array([tip_0, np.argmax(self.Dchosen[tip_0])])
        # we keep a list of the tips of each segment
        segs_tips = [tips_all]
        # print('init tips', segs_tips)
        segs_undecided = [True]
        segs_adjacency = [[]]
        segs_distances = np.zeros((1, 1))
        segs_adjacency_nodes = [{}]
        logg.info('... do not consider groups with less than {} points for splitting'
                  .format(self.min_group_size))
        for ibranch in range(self.n_splits):
            if 'unconstrained' in self.flavor:
                iseg, new_tips = self.select_segment(segs, segs_tips, segs_undecided)
                if iseg == -1:
                    logg.info('... partitioning converged')
                    break
                logg.info('... branching {}:'.format(ibranch + 1),
                          'split group', iseg)
                segs_distances = self.do_split(segs, segs_tips,
                                               segs_undecided,
                                               segs_adjacency,
                                               segs_distances,
                                               iseg, new_tips)
            else:
                logg.info('... split', ibranch + 1)
                stop, segs_distances = self.do_split_constrained(segs, segs_tips,
                                                                 segs_adjacency,
                                                                 segs_adjacency_nodes,
                                                                 segs_distances)
                if stop: break
        # store as class members
        self.segs = segs
        self.segs_tips = segs_tips
        self.segs_undecided = segs_undecided
        # the following is a bit too much, but this allows easy storage
        self.segs_adjacency = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
        for i, neighbors in enumerate(segs_adjacency):
            self.segs_adjacency[i, neighbors] = segs_distances[i][neighbors]
        self.segs_adjacency = self.segs_adjacency.tocsr()

        # print(self.segs_adjacency)
        # # print(segs_distances)
        # from .. import plotting as pl
        # pl.matrix(np.log1p(segs_distances), show=False)
        # import matplotlib.pyplot as pl
        # for i, neighbors in enumerate(segs_adjacency):
        #     pl.scatter([i for j in neighbors], neighbors, color='green')
        # pl.show()
        # pl.figure()
        # for i, ds in enumerate(segs_distances):
        #     ds = np.log1p(ds)
        #     x = [i for j, d in enumerate(ds) if i != j]
        #     y = [d for j, d in enumerate(ds) if i != j]
        #     pl.scatter(x, y, color='gray')
        #     neighbors = segs_adjacency[i]
        #     pl.scatter([i for j in neighbors],
        #                ds[neighbors], color='green')
        # pl.show()
        # print(self.segs_adjacency)
        # print([len(s) for s in self.segs])
        # if 'uncontract' not in self.flavor:
        #     self.contract_segments()
        # print(self.segs_adjacency)
        # print([len(s) for s in self.segs])
        # self.check_adjacency()
        # print('new')
        # print(self.segs_adjacency)

    def contract_segments(self):
        for i in range(1000):
            i, n = self.propose_segments_to_contract()
            if i != 0 or n != 0:
                G = nx.Graph(self.segs_adjacency)
                G_contracted = nx.contracted_nodes(G, n, i, self_loops=False)
                self.segs_adjacency = nx.to_scipy_sparse_matrix(G_contracted)
                self.segs[n] = np.append(self.segs[n], self.segs[i])
                self.segs.pop(i)
            else:
                break

    def propose_segments_to_contract(self):
        # nodes with two edges
        n_edges_per_seg = np.sum(self.segs_adjacency > 0, axis=1).A1
        for i in range(len(self.segs)):
            if n_edges_per_seg[i] == 2:
                neighbors = self.segs_adjacency[i].nonzero()[1]
                for neighbors_edges in range(1, 20):
                    for n_cnt, n in enumerate(neighbors):
                        if n_edges_per_seg[n] == neighbors_edges:
                            logg.info('merging segment', i, 'into', n, '(two edges)')
                            return i, n
        # nodes with a very small cell number into
        for i in range(len(self.segs)):
            if len(self.segs[i]) < self.min_group_size/2:
                neighbors = self.segs_adjacency[i].nonzero()[1]
                neighbor_sizes = [len(self.segs[n]) for n in neighbors]
                n = neighbors[np.argmax(neighbor_sizes)]
                logg.info('merging segment', i, 'into', n,
                          '(smaller than `min_group_size/2` = {})'.format(self.min_group_size/2))
                return i, n
        return 0, 0

    def check_adjacency(self):
        n_edges_per_seg = np.sum(self.segs_adjacency > 0, axis=1).A1
        for n_edges in range(1, np.max(n_edges_per_seg) + 1):
            for iseg in range(self.segs_adjacency.shape[0]):
                if n_edges_per_seg[iseg] == n_edges:
                    neighbor_segs = self.segs_adjacency[iseg].todense().A1
                    closest_points_other_segs = [seg[np.argmin(self.Dchosen[self.segs_tips[iseg][0], seg])]
                                                 for seg in self.segs]
                    seg = self.segs[iseg]
                    closest_points_in_segs = [seg[np.argmin(self.Dchosen[tips[0], seg])]
                                              for tips in self.segs_tips]
                    distance_segs = [self.Dchosen[closest_points_other_segs[ipoint], point]
                                     for ipoint, point in enumerate(closest_points_in_segs)]
                    # exclude the first point, the segment itself
                    closest_segs = np.argsort(distance_segs)[1:n_edges+1]
                    # update adjacency matrix within the loop!
                    # self.segs_adjacency[iseg, neighbor_segs > 0] = 0
                    # self.segs_adjacency[iseg, closest_segs] = np.array(distance_segs)[closest_segs]
                    # self.segs_adjacency[neighbor_segs > 0, iseg] = 0
                    # self.segs_adjacency[closest_segs, iseg] = np.array(distance_segs)[closest_segs].reshape(len(closest_segs), 1)
                    # n_edges_per_seg = np.sum(self.segs_adjacency > 0, axis=1).A1
                    # print(iseg, distance_segs, closest_segs)
                    # print(self.segs_adjacency)
        # self.segs_adjacency.eliminate_zeros()

    def do_split_constrained(self, segs, segs_tips, segs_adjacency, segs_adjacency_nodes,
                             segs_distances):
        isegs = np.argsort([len(seg) for seg in segs])[::-1]
        if len(segs[isegs[0]]) < self.min_group_size:
            return True, segs_distances
        for iseg in isegs:
            seg = segs[iseg]
            logg.info('... splitting group {} with size {}'.format(iseg, len(seg)))
            jsegs = [jseg for jseg in range(len(segs)) if jseg != iseg]
            dtip = np.zeros(len(seg))
            for jseg in jsegs:
                if len(segs_tips[jseg]) > 0:
                    jtip = segs_tips[jseg][0]
                    dtip += self.Dchosen[jtip, seg]
            if len(jsegs) > 0: dtip /= len(jsegs)
            itip = segs_tips[iseg][0]
            dtip += self.Dchosen[itip, seg]
            new_itip = np.argmax(dtip)
            if False:
                new_seg = np.ones(len(seg), dtype=bool)
                for jseg in range(len(segs)):
                    if len(segs_tips[jseg]) > 0:
                        jtip = segs_tips[jseg][0]
                        closer_to_jtip_than_to_new_itip = self.Dchosen[jtip, seg] < self.Dchosen[seg[new_itip], seg]
                        new_seg[closer_to_jtip_than_to_new_itip] = False
            else:
                new_seg = self.Dchosen[seg[new_itip], seg] < self.Dchosen[itip, seg]
            ssegs = [new_seg, ~new_seg]
            ssegs_tips = [[new_itip], []]
            l = len(np.flatnonzero(ssegs[0]))
            if l != 0 and len(seg) - l != 0:
                break
        logg.info('        obtained two groups with sizes {} and {}'
                  .format(l, len(seg)-l))
        trunk = 1
        # map back to global indices
        # print('        tips of split seg before map', ssegs_tips)
        for iseg_new, seg_new in enumerate(ssegs):
            ssegs[iseg_new] = seg[seg_new]
            ssegs_tips[iseg_new] = seg[ssegs_tips[iseg_new]]
        # add iseg tip
        if len(segs_tips[iseg]) > 0: ssegs_tips[trunk] = [segs_tips[iseg][0]]
        # print('        tips of split seg', ssegs_tips)
        # remove previous segment
        segs.pop(iseg)
        segs_tips.pop(iseg)
        # insert trunk at same position
        segs.insert(iseg, ssegs[trunk])
        segs_tips.insert(iseg, ssegs_tips[trunk])
        # append other segments
        segs += [seg for iseg, seg in enumerate(ssegs) if iseg != trunk]
        segs_tips += [seg_tips for iseg, seg_tips in enumerate(ssegs_tips) if iseg != trunk]
        # correct edges in adjacency matrix
        n_add = len(ssegs) - 1
        new_shape = (segs_distances.shape[0] + n_add, segs_distances.shape[1] + n_add)
        # segs_distances.resize() throws an error!
        segs_distances_help = segs_distances.copy()
        segs_distances = np.zeros((new_shape))
        segs_distances[np.ix_(range(segs_distances_help.shape[0]),
                              range(segs_distances_help.shape[1]))] = segs_distances_help
        segs_distances = self.adjust_adjacency(iseg, n_add,
                                               segs,
                                               segs_tips,
                                               segs_adjacency,
                                               segs_adjacency_nodes,
                                               segs_distances, iseg)
        return False, segs_distances

    def select_segment(self, segs, segs_tips, segs_undecided):
        """Out of a list of line segments, choose segment that has the most
        distant second data point.

        Assume the distance matrix Ddiff is sorted according to seg_idcs.
        Compute all the distances.

        Returns
        -------
        iseg : int
            Index identifying the position within the list of line segments.
        new_tips : int
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
            # map the global position to the position within the segment
            tips = [np.where(allindices[seg] == tip)[0][0]
                    for tip in segs_tips[iseg]]
            # find the third point on the segment that has maximal
            # added distance from the two tip points
            dseg = Dseg[tips[0]] + Dseg[tips[1]]
            third_tip = np.argmax(dseg)
            new_tips = np.append(tips, third_tip)
            # compute the score as ratio of the added distance to the third tip,
            # to what it would be if it were on the straight line between the
            # two first tips, given by Dseg[tips[:2]]
            # if we did not normalize, there would be a danger of simply
            # assigning the highest score to the longest segment
            if 'bi' in self.flavor:
                score = Dseg[new_tips[0], new_tips[1]]
            else:
                score = dseg[new_tips[2]] / Dseg[new_tips[0], new_tips[1]] * len(seg)
            score = len(seg) if self.choose_largest_segment else score  # simply the number of points
            # self.choose_largest_segment = False
            logg.m('... group', iseg, 'score', score, 'n_points', len(seg),
                   '(too small)' if len(seg) < self.min_group_size else '', v=4)
            if len(seg) <= self.min_group_size: score = 0
            # write result
            scores_tips[iseg, 0] = score
            scores_tips[iseg, 1:] = new_tips
        iseg = np.argmax(scores_tips[:, 0])
        if scores_tips[iseg, 0] == 0: return -1, None
        new_tips = scores_tips[iseg, 1:].astype(int)
        return iseg, new_tips

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

    def do_split(self, segs, segs_tips, segs_undecided, segs_adjacency,
                 segs_distances, iseg, new_tips):
        """Detect branching on given segment.

        Updates all list parameters inplace.

        Call function _do_split and perform bookkeeping on segs and
        segs_tips.

        Parameters
        ----------
        segs : list of np.ndarray
            Dchosen distance matrix restricted to segment.
        segs_tips : list of np.ndarray
            Stores all tip points for the segments in segs.
        iseg : int
            Position of segment under study in segs.
        new_tips : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.
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
        result = self._do_split(Dseg, new_tips, seg, segs_tips)
        ssegs, ssegs_tips, ssegs_adjacency, trunk = result
        # map back to global indices
        for iseg_new, seg_new in enumerate(ssegs):
            ssegs[iseg_new] = seg[seg_new]
            ssegs_tips[iseg_new] = seg[ssegs_tips[iseg_new]]
        # remove previous segment
        segs.pop(iseg)
        segs_tips.pop(iseg)
        # insert trunk at same position
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
        new_shape = (segs_distances.shape[0] + n_add, segs_distances.shape[1] + n_add)
        # segs_distances.resize() throws an error!
        segs_distances_help = segs_distances.copy()
        segs_distances = np.zeros((new_shape))
        segs_distances[np.ix_(range(segs_distances_help.shape[0]),
                              range(segs_distances_help.shape[1]))] = segs_distances_help
        segs_distances = self.adjust_adjacency(iseg, n_add,
                                               segs,
                                               segs_tips,
                                               segs_adjacency,
                                               segs_distances)
        segs_undecided += [False for i in range(n_add)]
        # need to return segs_distances as inplace formulation doesn't work
        return segs_distances

    def compute_attachedness(self, jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes):
        distances = []
        closest_points_in_jseg = []
        closest_points_in_kseg = []
        if self.attachedness_measure == 'd_single_pairwise':
            for kseg in kseg_list:
                reference_point_in_kseg = segs_tips[kseg][0]
                closest_points_in_jseg.append(segs[jseg][np.argmin(self.Dchosen[reference_point_in_kseg, segs[jseg]])])
                reference_point_in_jseg = closest_points_in_jseg[-1]
                closest_points_in_kseg.append(segs[kseg][np.argmin(self.Dchosen[reference_point_in_jseg, segs[kseg]])])
                distances.append(self.Dchosen[closest_points_in_jseg[-1], closest_points_in_kseg[-1]])
                logg.m('   ',
                       jseg, '(tip: {}, clos: {})'.format(segs_tips[jseg][0], closest_points_in_jseg[-1]),
                       kseg, '(tip: {}, clos: {})'.format(segs_tips[kseg][0], closest_points_in_kseg[-1]),
                       '->', distances[-1], v=4)
        elif self.attachedness_measure == 'd_full_pairwise':
            for kseg in kseg_list:
                closest_distance = 1e12
                closest_point_in_jseg = 0
                closest_point_in_kseg = 0
                for reference_point_in_kseg in segs[kseg]:
                    closest_point_in_jseg_test = segs[jseg][np.argmin(self.Dchosen[reference_point_in_kseg, segs[jseg]])]
                    if self.Dchosen[reference_point_in_kseg, closest_point_in_jseg_test] < closest_distance:
                        closest_point_in_jseg = closest_point_in_jseg_test
                        closest_point_in_kseg = reference_point_in_kseg
                        closest_distance = self.Dchosen[reference_point_in_kseg, closest_point_in_jseg_test]
                closest_points_in_kseg.append(closest_point_in_kseg)
                closest_points_in_jseg.append(closest_point_in_jseg)
                distances.append(closest_distance)
                logg.m('   ',
                       jseg, '(tip: {}, clos: {})'.format(segs_tips[jseg][0], closest_points_in_jseg[-1]),
                       kseg, '(tip: {}, clos: {})'.format(segs_tips[kseg][0], closest_points_in_kseg[-1]),
                       '->', distances[-1], v=4)
        elif self.attachedness_measure == 'ed_full_pairwise':
            for kseg in kseg_list:
                closest_similarity = 1e12
                closest_point_in_jseg = 0
                closest_point_in_kseg = 0
                for reference_point_in_kseg in segs[kseg]:
                    closest_point_in_jseg_test = segs[jseg][np.argmax(self.Ktilde[reference_point_in_kseg, segs[jseg]])]
                    if self.Ktilde[reference_point_in_kseg, closest_point_in_jseg_test] > closest_similarity:
                        closest_point_in_jseg = closest_point_in_jseg_test
                        closest_point_in_kseg = reference_point_in_kseg
                        closest_similarity = self.Ktilde[reference_point_in_kseg, closest_point_in_jseg_test]
                closest_points_in_kseg.append(closest_point_in_kseg)
                closest_points_in_jseg.append(closest_point_in_jseg)
                closest_distance = 1/closest_similarity
                distances.append(closest_distance)
                logg.m('   ',
                       jseg, '(tip: {}, clos: {})'.format(segs_tips[jseg][0], closest_points_in_jseg[-1]),
                       kseg, '(tip: {}, clos: {})'.format(segs_tips[kseg][0], closest_points_in_kseg[-1]),
                       '->', distances[-1], v=4)
        elif self.attachedness_measure == 'n_connecting_edges_brute_force':
            # this is a very slow implementation!!
            segs_jseg = set(segs[jseg])
            for kseg in kseg_list:
                connectedness = 0
                for reference_point_in_kseg in segs[kseg]:
                    for j in self.Ktilde[reference_point_in_kseg].nonzero()[1]:
                        if j in segs_jseg:
                            # print(reference_point_in_kseg, j, end=' | ')
                            connectedness += 1
                distances.append(1./(connectedness+1))
            logg.m(' ', jseg, '-', kseg_list, '->', distances, v=4)
        else:
            raise ValueError('unknown attachedness measure')
        return distances, closest_points_in_jseg, closest_points_in_kseg

    def trace_existing_connections(self, jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes, trunk):
        j_connects = segs_adjacency_nodes[jseg].copy()
        connectedness = [0, 0]
        not_trunk = 1 if trunk == 0 else 1
        kseg_trunk = set(segs[kseg_list[trunk]])
        kseg_not_trunk = set(segs[kseg_list[not_trunk]])
        for j_connect, connects in j_connects.items():
            for point_connect, seg_connect in connects:
                if seg_connect == kseg_list[trunk]:
                    if point_connect in kseg_trunk:
                        connectedness[0] += 1
                    elif point_connect in kseg_not_trunk:
                        if j_connect not in segs_adjacency_nodes[jseg]:
                            segs_adjacency_nodes[jseg][j_connect] = []
                        idx = segs_adjacency_nodes[jseg][j_connect].index((point_connect, kseg_list[trunk]))
                        segs_adjacency_nodes[jseg][j_connect][idx] = (point_connect, kseg_list[not_trunk])
                        if point_connect not in segs_adjacency_nodes[kseg_list[not_trunk]]:
                            segs_adjacency_nodes[kseg_list[not_trunk]][point_connect] = []
                        segs_adjacency_nodes[kseg_list[not_trunk]][point_connect].append((j_connect, jseg))
                        # clean up the dictionary for trunk
                        idx = segs_adjacency_nodes[kseg_list[trunk]][point_connect].index((j_connect, jseg))
                        segs_adjacency_nodes[kseg_list[trunk]][point_connect].pop(idx)
                        if len(segs_adjacency_nodes[kseg_list[trunk]][point_connect]) == 0:
                            del segs_adjacency_nodes[kseg_list[trunk]][point_connect]
                        connectedness[1] += 1
                    else:
                        print('should not occur!')
                        print('iseg is', kseg_list[trunk], 'and jseg is', jseg)
                        print('jseg connect is', j_connects[j_connect], 'from', j_connect)
                        print('point_connect is', point_connect, 'should be in', kseg_list[trunk])
                        for i, seg in enumerate(segs):
                            if point_connect in seg:
                                print('but found point_connect in', i)
                        raise ValueError('invalid state of algorithm')
        distances = [1/(1+c) for c in connectedness]
        print(' ', jseg, '-', kseg_list, '->', distances)
        return distances

    def establish_new_connections(self, kseg_list, segs, segs_adjacency_nodes):
        kseg_loop_idx = 0 if len(segs[kseg_list[0]]) < len(segs[kseg_list[1]]) else 1
        kseg_loop = kseg_list[kseg_loop_idx]
        kseg_test = kseg_list[0 if kseg_loop_idx == 1 else 1]
        seg_loop = segs[kseg_loop]
        seg_test = set(segs[kseg_test])
        connections = 0
        for p in seg_loop:
            p_neighbors = set(self.Ktilde[p].nonzero()[1])
            for q in p_neighbors:
                if q in seg_test:
                    if p not in segs_adjacency_nodes[kseg_loop]:
                        segs_adjacency_nodes[kseg_loop][p] = []
                    segs_adjacency_nodes[kseg_loop][p].append((q, kseg_test))
                    if q not in segs_adjacency_nodes[kseg_test]:
                        segs_adjacency_nodes[kseg_test][q] = []
                    segs_adjacency_nodes[kseg_test][q].append((p, kseg_loop))
                    connections += 1
                    if p in {22, 135} and q in {22, 135}:
                        print('!!!!', kseg_loop, kseg_test)
        distance = 1/(1+connections)
        print(' ', kseg_loop, '-', kseg_test, '->', distance)
        return distance

    def adjust_adjacency(self, iseg, n_add, segs, segs_tips, segs_adjacency, segs_adjacency_nodes,
                         segs_distances, trunk):
        prev_connecting_segments = segs_adjacency[iseg].copy()
        segs_adjacency += [[] for i in range(n_add)]
        segs_adjacency_nodes += [{} for i in range(n_add)]
        kseg_list = [iseg] + list(range(len(segs) - n_add, len(segs)))
        if self.attachedness_measure == 'n_connecting_edges':
            jseg_list = [jseg for jseg in range(len(segs)) if jseg not in kseg_list]
            for jseg in jseg_list:
                distances = self.trace_existing_connections(jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes, 0)
                segs_distances[jseg, kseg_list] = distances
                segs_distances[kseg_list, jseg] = distances
            distance = self.establish_new_connections(kseg_list, segs, segs_adjacency_nodes)
            segs_distances[kseg_list[0], kseg_list[1]] = distance
            segs_distances[kseg_list[1], kseg_list[0]] = distance
        logg.info('... treat existing connections')
        for jseg in prev_connecting_segments:
            if self.attachedness_measure != 'n_connecting_edges':
                result = self.compute_attachedness(jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes)
                distances, closest_points_in_jseg, closest_points_in_kseg = result
                segs_distances[jseg, kseg_list] = distances
                segs_distances[kseg_list, jseg] = distances
            idx = np.argmin(segs_distances[jseg, kseg_list])
            kseg_min = kseg_list[idx]
            pos = segs_adjacency[jseg].index(iseg)
            segs_adjacency[jseg][pos] = kseg_min
            pos_2 = segs_adjacency[iseg].index(jseg)
            segs_adjacency[iseg].pop(pos_2)
            segs_adjacency[kseg_min].append(jseg)
        # in case the segment we split should correspond to two "clusters", we
        # need to check whether the new segments connect to any of the other old
        # segments
        # if not, we add a link between the new segments, if yes, we add two
        # links to connect them at the correct old segments
        logg.info('... treat new connections')
        do_not_attach_kseg = False
        continue_after_distance_compute = False
        for kseg in kseg_list:
            jseg_list = [jseg for jseg in range(len(segs))
                         if jseg != kseg and jseg not in prev_connecting_segments]
            if self.attachedness_measure != 'n_connecting_edges':
                result = self.compute_attachedness(kseg, jseg_list, segs, segs_tips, segs_adjacency_nodes)
                distances, closest_points_in_kseg, closest_points_in_jseg = result
                segs_distances[kseg, jseg_list] = distances
                segs_distances[jseg_list, kseg] = distances
            if continue_after_distance_compute: continue
            idx = np.argmin(segs_distances[kseg, jseg_list])
            jseg_min = jseg_list[idx]
            if jseg_min not in kseg_list:
                segs_adjacency_sparse = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
                for i, neighbors in enumerate(segs_adjacency):
                    segs_adjacency_sparse[i, neighbors] = 1
                G = nx.Graph(segs_adjacency_sparse)
                paths_all = nx.single_source_dijkstra_path(G, source=kseg)
                if jseg_min not in paths_all:
                    segs_adjacency[jseg_min].append(kseg)
                    segs_adjacency[kseg].append(jseg_min)
                    logg.info('    attaching new segment', kseg, 'at', jseg_min)
                    # if we split the cluster, we should not attach kseg
                    do_not_attach_kseg = True
                else:
                    logg.info('    cannot attach new segment', kseg, 'at', jseg_min,
                              '(would produce cycle)')
                    if kseg != kseg_list[-1]:
                        logg.info('        continue')
                        continue
                    else:
                        logg.info('        do not add another link')
                        conintue_after_distance_compute = True
            if jseg_min in kseg_list and not do_not_attach_kseg:
                segs_adjacency[jseg_min].append(kseg)
                segs_adjacency[kseg].append(jseg_min)
                continue_after_distance_compute = True
                logg.info('    attaching new segment', kseg, 'with new segment', jseg_min)
        return segs_distances

    def _do_split(self, Dseg, tips, seg_reference, old_tips):
        """Detect branching on given segment.

        Call function __do_split three times for all three orderings of
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
        if 'tri' in self.flavor:
            ssegs = self._do_split_single_wolf17_tri(Dseg, tips)
        elif 'bi' in self.flavor:
            ssegs = self._do_split_single_wolf17_bi(Dseg, tips)
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
            secondtip = newseg[np.argmax(Dseg[tips[inewseg]][newseg])]
            ssegs_tips.append([tips[inewseg], secondtip])
        undecided_cells = np.arange(Dseg.shape[0], dtype=int)[nonunique]
        if len(undecided_cells) > 0:
            ssegs.append(undecided_cells)
            # establish the connecting points with the other segments
            for inewseg, newseg_tips in enumerate(ssegs_tips):
                reference_point = newseg_tips[0]
                # closest cell to the new segment within undecided cells
                closest_cell = undecided_cells[np.argmin(Dseg[reference_point][undecided_cells])]
                # closest cell to the undecided cells within new segment
                closest_cell = ssegs[inewseg][np.argmin(Dseg[closest_cell][ssegs[inewseg]])]
            # also compute tips for the undecided cells
            tip_0 = undecided_cells[np.argmax(Dseg[undecided_cells[0]][undecided_cells])]
            tip_1 = undecided_cells[np.argmax(Dseg[tip_0][undecided_cells])]
            ssegs_tips.append([tip_0, tip_1])
            ssegs_adjacency = [[3], [3], [3], [0, 1, 2]]
            trunk = 3
            # import matplotlib.pyplot as pl
            # for iseg_new, seg_new in enumerate(ssegs):
            #     pl.figure()
            #     pl.scatter(self.passed_adata.smp['X_diffmap'][:, 0], self.passed_adata.smp['X_diffmap'][:, 1], s=1, c='grey')
            #     pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][seg_new, 0], self.passed_adata.smp['X_diffmap'][seg_reference][seg_new, 1], marker='x', s=2, c='blue')
            #     # pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][tips[iseg_new], 0], self.passed_adata.smp['X_diffmap'][seg_reference][tips[iseg_new], 1], marker='x', c='black')
            #     # pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][second_tip[iseg_new], 0], self.passed_adata.smp['X_diffmap'][seg_reference][second_tip[iseg_new], 1], marker='o', c='black')
            #     pl.xticks([])
            #     pl.yticks([])
            #     # pl.savefig('./figs/cutting_off_tip={}.png'.format(iseg_new))
            # pl.show()
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
            # print(ssegs_adjacency)
            # import matplotlib.pyplot as pl
            # for iseg_new, seg_new in enumerate(ssegs):
            #     pl.figure()
            #     pl.scatter(self.passed_adata.smp['X_diffmap'][:, 0], self.passed_adata.smp['X_diffmap'][:, 1], s=1, c='grey')
            #     pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][seg_new, 0], self.passed_adata.smp['X_diffmap'][seg_reference][seg_new, 1], marker='x', s=2, c='blue')
            #     # pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][tips[iseg_new], 0], self.passed_adata.smp['X_diffmap'][seg_reference][tips[iseg_new], 1], marker='x', c='black')
            #     # pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][second_tip[iseg_new], 0], self.passed_adata.smp['X_diffmap'][seg_reference][second_tip[iseg_new], 1], marker='o', c='black')
            #     for i in range(3):
            #         if i != iseg_new:
            #             pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][closest_points[iseg_new, i], 0],
            #                        self.passed_adata.smp['X_diffmap'][seg_reference][closest_points[iseg_new, i], 1], marker='o', c='black')
            #             pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][closest_points[i, iseg_new], 0],
            #                        self.passed_adata.smp['X_diffmap'][seg_reference][closest_points[i, iseg_new], 1], marker='x', c='black')
            #     pl.xticks([])
            #     pl.yticks([])
            #     # pl.savefig('./figs/cutting_off_tip={}.png'.format(iseg_new))
            # pl.show()
            # print('trunk', trunk)
        else:
            trunk = 0
            ssegs_adjacency = [[1], [0]]
            reference_point_in_0 = ssegs_tips[0][0]
            closest_point_in_1 = ssegs[1][np.argmin(Dseg[reference_point_in_0][ssegs[1]])]
            reference_point_in_1 = closest_point_in_1  # ssegs_tips[1][0]
            closest_point_in_0 = ssegs[0][np.argmin(Dseg[reference_point_in_1][ssegs[0]])]
        return ssegs, ssegs_tips, ssegs_adjacency, trunk

    def _do_split_single_haghverdi16(self, Dseg, tips):
        """Detect branching on given segment.
        """
        # compute splits using different starting points the first index of
        # tips is the starting point for the other two, the order does not
        # matter
        ssegs = []
        # permutations of tip cells
        ps = [[0, 1, 2],  # start by computing distances from the first tip
              [1, 2, 0],  #             -"-                       second tip
              [2, 0, 1]]  #             -"-                       third tip
        # import matplotlib.pyplot as pl
        for i, p in enumerate(ps):
            ssegs.append(self.__do_split_haghverdi16(Dseg, tips[p]))
        return ssegs

    def _do_split_single_wolf17_tri(self, Dseg, tips):
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

    def _do_split_single_wolf17_bi(self, Dseg, tips):
        dist_from_0 = Dseg[tips[0]]
        dist_from_1 = Dseg[tips[1]]
        if True:
            closer_to_0_than_to_1 = dist_from_0 < dist_from_1
            ssegs = [closer_to_0_than_to_1, ~closer_to_0_than_to_1]
        else:
            time = dist_from_0 - dist_from_1
            idcs = np.argsort(time)
            i = np.argmax(np.diff(time[idcs]))
            ssegs = [idcs[:i+1], idcs[i+1:]]
        return ssegs

    def __do_split_haghverdi16(self, Dseg, tips):
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
            logg.m('... is root itself, never obtain significant correlation', v=4)
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

# Author: F. Alex Wolf (http://falexwolf.de)
"""
"""

import sys, os
import numpy as np
import scipy as sp
import networkx as nx
import scipy.sparse
from .. import logging as logg
from ..data_structs import ann_data
from ..data_structs import data_graph
from .louvain import louvain
from ..plotting import utils as pl_utils


MINIMAL_REALIZED_ATTACHEDNESS = 0.05


def aga(adata,
        node_groups='louvain',
        n_nodes=None,
        n_neighbors=30,
        n_pcs=50,
        n_dcs=10,
        resolution=1,
        recompute_pca=False,
        recompute_distances=False,
        recompute_graph=False,
        recompute_louvain=False,
        attachedness_measure='random_walk',
        n_jobs=None,
        copy=False):
    """Approximate Graph Abstraction

    Infer the relations of subgroups in the data through approximate graph
    abstraction. The result is a much simpler graph where each node corresponds
    to a cell subgroup. The tree induces an ordering between nodes and the cells
    are ordered within each node.

    Reference
    ---------
    Wolf et al., bioRxiv (2017)

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata.add['xroot'] : np.ndarray
            Root of stochastic process on data points (root cell), specified
            as expression vector of shape adata.n_smps.
        adata.smp['X_pca']: np.ndarray
            PCA representation of the data matrix (result of preprocessing with
            PCA). Will be used if option `recompute_pca` is False.
        adata.smp['X_diffmap']: np.ndarray
            Diffmap representation of the data matrix (result of running
            `diffmap`). Will be used if option `recompute_graph` is False.
    node_groups : any categorical sample annotation or {'louvain', 'segments'}, optional (default: 'louvain')
        Criterion to determine the resoluting partitions of the
        graph/data. 'louvain' uses the louvain algorithm and optimizes
        modularity of the graph, 'segments' uses a bipartioning
        criterium that is loosely inspired by hierarchical clustering. You can
        also pass your predefined groups by choosing any sample annotation.
    n_nodes : int or None, optional (default: None)
        Number of nodes in the abstracted graph. Except when choosing
        'segments' for `node_groups`, for which `n_nodes` defaults to
        `n_nodes=1`, `n_nodes` defaults to the number of groups implied by the
        choice of `node_groups`.
    n_neighbors : int or None, optional (default: None)
        Number of nearest neighbors on the knn graph. See the default
        of data_graph (usually 30).
    n_pcs : int, optional (default: 50)
        Use n_pcs PCs to compute the Euclidian distance matrix, which is the
        basis for generating the graph. Set to 0 if you don't want preprocessing
        with PCA.
    n_dcs : int, optional (default: 10)
        Number of diffusion components (very similar to eigen vectors of
        adjacency matrix) to use for distance computations.
    resolution : float, optional (default: 1.0)
        For Louvain algorithm. Note that you should set `recompute_louvain` to
        True if changing this to recompute.
    recompute_graph : bool, optional (default: False)
        Recompute single-cell graph. Only then `n_neighbors` has an effect if
        there is already a cached `distance` or `X_diffmap` in adata.
    recompute_pca : bool, optional (default: False)
        Recompute PCA.
    recompute_louvain : bool, optional (default: False)
        When changing the `resolution` parameter, you should set this to True.
    attachedness_measure : {'connectedness', 'random_walk'}, optional (default: 'random_walk')
        How to measure attachedness.
    n_jobs : int or None (default: None)
        Number of cpus to use for parallel processing (default: sett.n_jobs).
    copy : bool, optional (default: False)
        Copy instance before computation and return a copy. Otherwise, perform
        computation inplace and return None.

    Notes
    -----
    Writes the following.
        aga_adjacency : sparse csr matrix
            Array of dim (number of samples) that stores the pseudotime of each
            cell, that is, the DPT distance with respect to the root cell.
        aga_groups : np.ndarray of dtype string
            Array of dim (number of samples) that stores the subgroup id ('0',
            '1', ...) for each cell.
        aga_pseudotime : np.ndarray of dtype float
            Array of dim (number of samples) that stores a pseudotime from a
            root node, if the latter was passed.
    """
    adata = adata.copy() if copy else adata
    fresh_compute_louvain = False
    if (node_groups == 'louvain'
        and ('louvain_groups' not in adata.smp_keys()
             or recompute_louvain
             or not data_graph.no_recompute_of_graph_necessary(
            adata,
            recompute_pca=recompute_pca,
            recompute_distances=recompute_distances,
            recompute_graph=recompute_graph,
            n_neighbors=n_neighbors,
            n_dcs=n_dcs))):
        louvain(adata,
                resolution=resolution,
                n_neighbors=n_neighbors,
                recompute_pca=recompute_pca,
                recompute_graph=recompute_graph,
                n_pcs=n_pcs,
                n_dcs=n_dcs)
        fresh_compute_louvain = True
    clusters = node_groups
    if node_groups == 'louvain': clusters = 'louvain_groups'
    logg.info('running Approximate Graph Abstraction (AGA)', r=True)
    if ('iroot' not in adata.add
        and 'xroot' not in adata.add
        and 'xroot' not in adata.var):
        logg.info('    no root cell found, no computation of pseudotime')
        msg = \
    '''To enable computation of pseudotime, pass the index or expression vector
    of a root cell. Either add
        adata.add['iroot'] = root_cell_index
    or (robust to subsampling)
        adata.var['xroot'] = adata.X[root_cell_index, :]
    where "root_cell_index" is the integer index of the root cell, or
        adata.var['xroot'] = adata[root_cell_name, :].X
    where "root_cell_name" is the name (a string) of the root cell.'''
        logg.hint(msg)
    aga = AGA(adata,
              clusters=clusters,
              n_neighbors=n_neighbors,
              n_pcs=n_pcs,
              n_dcs=n_dcs,
              min_group_size=20/resolution,
              n_jobs=n_jobs,
              # we do not need to recompute things both in the louvain
              # call above and here
              recompute_graph=recompute_graph and not fresh_compute_louvain,
              recompute_distances=recompute_distances and not fresh_compute_louvain,
              recompute_pca=recompute_pca and not fresh_compute_louvain,
              n_nodes=n_nodes,
              attachedness_measure=attachedness_measure)
    updated_diffmap = aga.update_diffmap()
    adata.smp['X_diffmap'] = aga.rbasis[:, 1:]
    adata.smp['X_diffmap0'] = aga.rbasis[:, 0]
    adata.add['diffmap_evals'] = aga.evals[1:]
    adata.add['data_graph_distance_local'] = aga.Dsq
    adata.add['data_graph_norm_weights'] = aga.Ktilde
    if aga.iroot is not None:
        aga.set_pseudotime()  # pseudotimes are random walk distances from root point
        adata.add['iroot'] = aga.iroot  # update iroot, might have changed when subsampling, for example
        adata.smp['aga_pseudotime'] = aga.pseudotime
    # detect splits and partition the data into segments
    aga.splits_segments()
    # vector of length n_samples of group names
    adata.smp['aga_groups'] = aga.segs_names.astype('U')
    # vectors of length n_groups
    adata.add['aga_groups_names'] = np.array([str(n) for n in aga.segs_names_unique])
    adata.add['aga_groups_sizes'] = aga.segs_sizes
    # the ordering according to groups and pseudotime
    adata.smp['aga_indices'] = aga.indices
    # the changepoints - marking different segments - in the ordering above
    adata.add['aga_changepoints'] = aga.changepoints
    # the tip points of segments
    # adata.add['aga_grouptips'] = aga.segs_tips
    # the tree/graph adjacency matrix
    adata.add['aga_adjacency'] = aga.segs_adjacency
    if fresh_compute_louvain:
        adata.smp['louvain_groups'] = adata.smp['aga_groups']
        adata.add['louvain_groups_names'] = adata.add['aga_groups_names']
    if (clusters not in {'segments', 'unconstrained_segments'}
        and not fresh_compute_louvain):
        adata.add['aga_groups_original'] = clusters
        adata.add['aga_groups_names_original'] = np.array(aga.segs_names_original)
        if clusters + '_colors' not in adata.add:
            pl_utils.add_colors_for_categorical_sample_annotation(adata, clusters)
        colors_original = []
        if clusters + '_names' not in adata.add:
            from natsort import natsorted
            adata.add[clusters + '_names'] = natsorted(np.unique(adata.smp[clusters]))
        name_list = list(adata.add[clusters + '_names'])
        for name in aga.segs_names_original:
            idx = name_list.index(name)
            colors_original.append(adata.add[clusters + '_colors'][idx])
        adata.add['aga_groups_colors_original'] = np.array(colors_original)
    adata.add['aga_distances'] = aga.segs_distances
    adata.add['aga_attachedness'] = aga.segs_attachedness
    adata.add['aga_attachedness_absolute'] = aga.segs_attachedness_absolute
    adata.add['aga_adjacency_absolute'] = aga.segs_adjacency_absolute
    logg.info('    finished', t=True, end=' ')
    logg.info('and added\n'
              '    "aga_adjacency", adjacency matrix defining the abstracted graph (adata.add),\n'
              '    "aga_groups", groups corresponding to nodes of abstracted graph (adata.smp)'
              + (',\n    "aga_pseudotime", pseudotime with respect to root cell (adata.smp)' if aga.iroot is not None else ''))
    return adata if copy else None


def aga_contract_graph(adata, min_group_size=0.01, max_n_contractions=1000, copy=False):
    """Contract the abstracted graph.
    """
    adata = adata.copy() if copy else adata
    if 'aga_adjacency' not in adata.add: raise ValueError('run tool aga first!')
    min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * adata.n_smps)
    logg.info('contract graph using `min_group_size={}`'.format(min_group_size))

    def propose_nodes_to_contract(adjacency, node_groups):
        # nodes with two edges
        n_edges_per_seg = np.sum(adjacency > 0, axis=1).A1
        for i in range(adjacency.shape[0]):
            if n_edges_per_seg[i] == 2:
                neighbors = adjacency[i].nonzero()[1]
                for neighbors_edges in range(1, 20):
                    for n_cnt, n in enumerate(neighbors):
                        if n_edges_per_seg[n] == neighbors_edges:
                            logg.m('merging node {} into {} (two edges)'
                                   .format(i, n), v=4)
                            return i, n
        # node groups with a very small cell number
        for i in range(adjacency.shape[0]):
            if node_groups[str(i) == node_groups].size < min_group_size:
                neighbors = adjacency[i].nonzero()[1]
                neighbor_sizes = [node_groups[str(n) == node_groups].size for n in neighbors]
                n = neighbors[np.argmax(neighbor_sizes)]
                logg.m('merging node {} into {} '
                       '(smaller than `min_group_size` = {})'
                       .format(i, n, min_group_size), v=4)
                return i, n
        return 0, 0

    def contract_nodes(adjacency, node_groups):
        for count in range(max_n_contractions):
            i, n = propose_nodes_to_contract(adjacency, node_groups)
            if i != 0 or n != 0:
                G = nx.Graph(adjacency)
                G_contracted = nx.contracted_nodes(G, n, i, self_loops=False)
                adjacency = nx.to_scipy_sparse_matrix(G_contracted)
                node_groups[str(i) == node_groups] = str(n)
                for j in range(i+1, G.size()+1):
                    node_groups[str(j) == node_groups] = str(j-1)
            else:
                break
        return adjacency, node_groups

    size_before = adata.add['aga_adjacency'].shape[0]
    adata.add['aga_adjacency'], adata.smp['aga_groups'] = contract_nodes(
        adata.add['aga_adjacency'], adata.smp['aga_groups'])
    adata.add['aga_groups_names'] = np.unique(adata.smp['aga_groups'])
    for key in ['aga_attachedness', 'aga_groups_original',
                'aga_groups_names_original', 'aga_groups_colors_original']:
        if key in adata.add: del adata.add[key]
    logg.info('    contracted graph from {} to {} nodes'
              .format(size_before, adata.add['aga_adjacency'].shape[0]))
    logg.m('removed adata.add["aga_attachedness"]', v=4)
    return adata if copy else None


class AGA(data_graph.DataGraph):
    """Approximate Graph Abstraction
    """

    def __init__(self,
                 adata,
                 n_nodes=None,
                 n_neighbors=30,
                 n_pcs=50,
                 n_dcs=10,
                 min_group_size=20,
                 recompute_pca=False,
                 recompute_distances=False,
                 recompute_graph=False,
                 attachedness_measure='connectedness',
                 clusters=None,
                 n_jobs=1):
        super(AGA, self).__init__(adata,
                                  k=n_neighbors,
                                  n_pcs=n_pcs,
                                  n_dcs=n_dcs,
                                  n_jobs=n_jobs,
                                  recompute_pca=recompute_pca,
                                  recompute_distances=recompute_distances,
                                  recompute_graph=recompute_graph)
        self.n_neighbors = n_neighbors
        self.min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * self.X.shape[0])
        self.passed_adata = adata  # just for debugging purposes
        self.choose_largest_segment = True
        self.attachedness_measure = attachedness_measure
        self.clusters = clusters
        self.clusters_precomputed = None
        self.clusters_precomputed_names = None
        self.flavor_develop = 'bi'  # bipartitioning
        if clusters not in {'segments', 'unconstrained_segments'}:
            if clusters not in adata.smp_keys():
                raise ValueError('Did not find {} in adata.smp_keys()! '
                                 'If you do not have any precomputed clusters, pass "segments" for "node_groups" instead'
                                 .format(clusters))
            clusters_array = adata.smp[clusters]
            # transform to a list of index arrays
            self.clusters_precomputed = []
            # TODO: this is not a good solution
            if clusters + '_names' in adata.add:
                self.clusters_precomputed_names = list(adata.add[clusters + '_names'])
            else:
                self.clusters_precomputed_names = []
            from natsort import natsorted
            for cluster_name in natsorted(np.unique(clusters_array)):
                self.clusters_precomputed.append(np.where(cluster_name == clusters_array)[0])
                if clusters + '_names' not in adata.add:
                    self.clusters_precomputed_names.append(cluster_name)
            n_nodes = len(self.clusters_precomputed)
        else:
            if n_nodes is None:
                n_nodes = 1
                logg.hint(
                    'by passing the parameter `n_nodes`, '
                    'choose the number of subgroups to detect')
        self.n_splits = n_nodes - 1

    def splits_segments(self):
        """Detect splits and partition the data into corresponding segments.

        Detect all splits up to `n_nodes`.

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
        """Detect all splits up to `n_nodes`.

        Writes Attributes
        -----------------
        segs : np.ndarray
            List of integer index arrays.
        segs_tips : np.ndarray
            List of indices of the tips of segments.
        """
        logg.info('    abstracted graph will have {} nodes'.format(self.n_splits+1))
        indices_all = np.arange(self.X.shape[0], dtype=int)
        segs = [indices_all]
        if False:  # this is safe, but not compatible with on-the-fly computation
            tips_all = np.array(np.unravel_index(np.argmax(self.Dchosen), self.Dchosen.shape))
        else:
            if self.iroot is not None:
                tip_0 = np.argmax(self.Dchosen[self.iroot])
            else:
                tip_0 = np.argmax(self.Dchosen[0])  # just a random index, here fixed to "0"
            tips_all = np.array([tip_0, np.argmax(self.Dchosen[tip_0])])
        # we keep a list of the tips of each segment
        segs_tips = [tips_all]
        if self.clusters_precomputed_names:
            self.segs_names_original = [', '.join(self.clusters_precomputed_names)]
        segs_undecided = [True]
        segs_adjacency = [[]]
        segs_distances = np.zeros((1, 1))
        segs_adjacency_nodes = [{}]
        # logg.info('    do not consider groups with less than {} points for splitting'
        #           .format(self.min_group_size))
        for ibranch in range(self.n_splits):
            if self.clusters == 'unconstrained_segments':
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
                logg.m('    split', ibranch + 1, v=4)
                stop, segs_distances = self.do_split_constrained(segs, segs_tips,
                                                                 segs_adjacency,
                                                                 segs_adjacency_nodes,
                                                                 segs_distances)
                if stop: break
        # segs, segs_tips, segs_distances, segs_adjacency = self.compute_tree_from_clusters()
        # store as class members
        self.segs = segs
        self.segs_tips = segs_tips
        self.segs_sizes = []
        for iseg, seg in enumerate(self.segs):
            self.segs_sizes.append(len(seg))
        # self.segs_undecided = segs_undecided
        # the following is a bit too much, but this allows easy storage
        realized_distances = []
        for i, neighbors in enumerate(segs_adjacency):
            realized_distances += segs_distances[i][neighbors].tolist()

        median_realized_distances = np.median(realized_distances)
        self.segs_attachedness = np.zeros_like(segs_distances)
        self.segs_attachedness[segs_distances <= median_realized_distances] = 1
        self.segs_attachedness[segs_distances > median_realized_distances] = (
            np.exp(-(segs_distances-median_realized_distances)/median_realized_distances)
            [segs_distances > median_realized_distances])
        self.segs_distances = segs_distances
        self.segs_attachedness_absolute = 1/segs_distances
        if self.attachedness_measure == 'connectedness':
            norm = np.sqrt(np.multiply.outer(self.segs_sizes, self.segs_sizes))
            self.segs_attachedness_absolute /= norm

        minimal_realized_attachedness = MINIMAL_REALIZED_ATTACHEDNESS
        self.segs_adjacency = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
        self.segs_adjacency_absolute = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
        for i, neighbors in enumerate(segs_adjacency):
            clipped_attachedness = self.segs_attachedness[i][neighbors]
            clipped_attachedness[clipped_attachedness < minimal_realized_attachedness] = minimal_realized_attachedness
            self.segs_adjacency[i, neighbors] = clipped_attachedness
            self.segs_attachedness[i, neighbors] = clipped_attachedness
            self.segs_adjacency_absolute[i, neighbors] = self.segs_attachedness_absolute[i, neighbors]
        self.segs_adjacency = self.segs_adjacency.tocsr()

    def do_split_constrained(self, segs, segs_tips,
                             segs_adjacency,
                             segs_adjacency_nodes,
                             segs_distances):

        if max([len(seg) for seg in segs]) < self.min_group_size:
            return True, segs_distances

        def binary_split_largest():
            isegs = np.argsort([len(seg) for seg in segs])[::-1]
            for iseg in isegs:
                seg = segs[iseg]
                logg.m('    splitting group {} with size {}'.format(iseg, len(seg)), v=4)
                jsegs = [jseg for jseg in range(len(segs)) if jseg != iseg]
                dtip = np.zeros(len(seg))
                for jseg in jsegs:
                    if len(segs_tips[jseg]) > 0:
                        jtip = segs_tips[jseg][0]
                        dtip += self.Dchosen[jtip, seg]
                if len(jsegs) > 0: dtip /= len(jsegs)
                itip = segs_tips[iseg][0]
                dtip += self.Dchosen[itip, seg]
                imax = np.argmax(dtip)
                dist_new_itip = dtip[imax]
                new_itip = seg[imax]
                new_seg = self.Dchosen[new_itip, seg] < self.Dchosen[itip, seg]
                ssegs = [seg[new_seg], seg[~new_seg]]
                ssegs_tips = [[new_itip], []]
                sizes = [len(ssegs[0]), len(ssegs[1])]
                if sizes[0] != 0 and sizes[1] != 0: break
            logg.m('    new tip {} with distance {:.6}, constraint was {}'
                   .format(new_itip, dist_new_itip, itip), v=4)
            logg.m('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            if len(segs_tips[iseg]) > 0: ssegs_tips[1] = [segs_tips[iseg][0]]
            return iseg, seg, ssegs, ssegs_tips, sizes

        def new_split(segs_tips):
            # upon initialization, start with no tips
            if len(segs) == 1:
                segs_tips.pop(0)
                segs_tips.append([])
            scores = []
            new_tips = []
            second_tips = []
            third_tips = []
            for iseg, seg in enumerate(segs):
                seg = segs[iseg]
                if len(seg) <= self.min_group_size:
                    scores.append(-1)
                    new_tips.append(0)
                    second_tips.append(0)
                    third_tips.append(0)
                    continue
                jsegs = [jseg for jseg in range(len(segs)) if jseg != iseg]
                dtip_others = np.zeros(len(seg))
                for jseg in jsegs:
                    if len(segs_tips[jseg]) > 0:
                        jtip = segs_tips[jseg][0]
                        dtip_others += self.Dchosen[jtip, seg]
                if len(jsegs) > 0: dtip_others /= len(jsegs)
                dtip = dtip_others
                need_to_compute_another_tip = False
                if len(segs_tips[iseg]) > 0:
                    itip = segs_tips[iseg][0]
                    dtip += self.Dchosen[itip, seg]
                elif len(jsegs) == 0:
                    # just take a random point and the extremum with respect to that
                    # point, the point is fixed to be the first in the segment
                    itip = seg[np.argmax(self.Dchosen[seg[0], seg])]
                    dtip += self.Dchosen[itip, seg]
                else:
                    need_to_compute_another_tip = True
                new_itip = seg[np.argmax(dtip)]
                if need_to_compute_another_tip:
                    itip = seg[np.argmax(self.Dchosen[new_itip, seg])]
                dtip = self.Dchosen[itip, seg] + self.Dchosen[new_itip, seg]
                itip_third = np.argmax(dtip)
                # score = dtip[itip_third] / self.Dchosen[itip, new_itip]
                score = len(seg)
                scores.append(score)
                new_tips.append(new_itip)
                second_tips.append(itip)
                third_tips.append(seg[itip_third])
            iseg = np.argmax(scores)
            new_itip = new_tips[iseg]
            itip = second_tips[iseg]
            third_itip = third_tips[iseg]
            seg = segs[iseg]
            logg.m('... splitting group {} with size {}'.format(iseg, len(seg)), v=4)
            new_seg = self.Dchosen[new_itip, seg] < self.Dchosen[itip, seg]
            size_0 = np.sum(new_seg)
            if False:
                if size_0 > len(seg) - size_0 and len(segs) == 1:
                    new_itip = itip
                    new_seg = ~new_seg
                    size_0 = len(seg) - size_0
                idcs = np.argsort(self.Dchosen[new_itip, seg])
                sorted_dists_from_new_tip = self.Dchosen[new_itip, seg][idcs]
                i = np.argmax(np.diff(sorted_dists_from_new_tip))
                if i <= size_0: new_seg[idcs[i+1:]] = False  # idx starts at zero and this works
            ssegs = [seg[new_seg], seg[~new_seg]]
            ssegs_tips = [[new_itip], []]
            sizes = [len(ssegs[0]), len(ssegs[1])]
            logg.m('    new tip {} with distance {:.6}, constraint was {}'
                   .format(new_itip, 0.0, itip), v=4)
            logg.m('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            logg.m('    the scores where', scores, v=4)
            return iseg, seg, ssegs, ssegs_tips, sizes

        def star_split(segs_tips):
            if len(segs) == 1:
                segs_tips.pop(0)
                segs_tips.append([])
            isegs = np.argsort([len(seg) for seg in segs])[::-1]
            iseg = isegs[0]
            seg = segs[iseg]
            new_tips = [seg[np.argmax(self.Dchosen[seg[0], seg])]]
            dtip_others = self.Dchosen[new_tips[0], seg]
            dists = [np.max(dtip_others)]
            for j in range(10):
                new_tip = seg[np.argmax(dtip_others)]
                if new_tip in new_tips: break
                new_tips.append(new_tip)
                dtip_j = self.Dchosen[new_tips[-1], seg]
                dists.append(np.max(dtip_j))
                dtip_others += dtip_j
            tip_idx_max = np.argmax(dists)
            new_tip = new_tips.pop(tip_idx_max)
            dist_max = dists.pop(tip_idx_max)
            new_seg = np.ones(len(seg), dtype=bool)
            for constraint_tip in new_tips:
                new_seg[self.Dchosen[new_tip, seg] > self.Dchosen[constraint_tip, seg]] = False
            ssegs = [seg[new_seg], seg[~new_seg]]
            ssegs_tips = [[new_tip], new_tips]
            sizes = [len(ssegs[0]), len(ssegs[1])]
            np.set_printoptions(precision=4)
            logg.m('    new tip', new_tip, 'with distance', dist_max,
                   'using constraints {} with distances'
                   .format(new_tips), v=4)
            logg.m('   ', dists, v=4)
            logg.m('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            return iseg, seg, ssegs, ssegs_tips, sizes

        def select_precomputed(segs_tips):
            if len(segs) == 1:
                segs_tips.pop(0)
                segs_tips.append([])
            iseg = 0
            seg = segs[iseg]
            logg.m('    splitting group {} with size {}'.format(iseg, len(seg)), v=4)
            new_tips = [seg[np.argmax(self.Dchosen[seg[0], seg])]]
            dtip_others = self.Dchosen[new_tips[0], seg]
            dists = [np.max(dtip_others)]
            for j in range(10):
                new_tip = seg[np.argmax(dtip_others)]
                if new_tip in new_tips: break
                new_tips.append(new_tip)
                dtip_j = self.Dchosen[new_tips[-1], seg]
                dists.append(np.max(dtip_j))
                dtip_others += dtip_j
            tip_idx_max = np.argmax(dists)
            new_tip = new_tips.pop(tip_idx_max)
            dist_max = dists.pop(tip_idx_max)
            for iclus, clus in enumerate(self.clusters_precomputed):
                if new_tip in set(clus):
                    new_seg = clus
                    clus_name = self.clusters_precomputed_names[iclus]
                    break
            pos_new_seg = np.in1d(seg, new_seg, assume_unique=True)
            ssegs = [new_seg, seg[~pos_new_seg]]
            ssegs_tips = [[new_tip], new_tips]
            sizes = [len(ssegs[0]), len(ssegs[1])]
            np.set_printoptions(precision=4)
            logg.m('    new tip', new_tip, 'with distance', dist_max,
                   'using constraints {} with distances'
                   .format(new_tips), v=4)
            logg.m('   ', dists, v=4)
            logg.m('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            return iseg, seg, ssegs, ssegs_tips, sizes, clus_name

        if self.clusters_precomputed is None:
            iseg, seg, ssegs, ssegs_tips, sizes = binary_split_largest()
            # iseg, seg, ssegs, ssegs_tips, sizes = new_split(segs_tips)
            # iseg, seg, ssegs, ssegs_tips, sizes = star_split(segs_tips)
        else:
            iseg, seg, ssegs, ssegs_tips, sizes, clus_name = select_precomputed(segs_tips)
        trunk = 1
        segs.pop(iseg)
        segs_tips.pop(iseg)
        # insert trunk at same position
        segs.insert(iseg, ssegs[trunk])
        segs_tips.insert(iseg, ssegs_tips[trunk])
        if self.clusters_precomputed_names:
            iseg_name = ' '.join(np.setdiff1d(self.clusters_precomputed_names,
                                              [n for n in self.segs_names_original]
                                              + [clus_name]))
            self.segs_names_original[iseg] = iseg_name
        # append other segments
        segs += [seg for iseg, seg in enumerate(ssegs) if iseg != trunk]
        segs_tips += [seg_tips for iseg, seg_tips in enumerate(ssegs_tips) if iseg != trunk]
        if self.clusters_precomputed_names: self.segs_names_original += [clus_name]
        # correct edges in adjacency matrix
        n_add = len(ssegs) - 1
        new_shape = (segs_distances.shape[0] + n_add, segs_distances.shape[1] + n_add)
        # segs_distances.resize() throws an error!
        segs_distances_help = segs_distances.copy()
        segs_distances = np.zeros((new_shape))
        segs_distances[np.ix_(range(segs_distances_help.shape[0]),
                              range(segs_distances_help.shape[1]))] = segs_distances_help
        segs_distances = self.adjust_adjacency(iseg,
                                               n_add,
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
            if 'bi' == self.flavor_develop:
                score = Dseg[new_tips[0], new_tips[1]]
            elif 'tri' == self.flavor_develop:
                score = dseg[new_tips[2]] / Dseg[new_tips[0], new_tips[1]] * len(seg)
            else:
                raise ValueError('unknown `self.flavor_develop`')
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

    def compute_attachedness(self, jseg, kseg_list, segs, segs_tips,
                             segs_adjacency_nodes):
        distances = []
        median_distances = []
        measure_points_in_jseg = []
        measure_points_in_kseg = []
        if self.attachedness_measure == 'random_walk_approx':
            for kseg in kseg_list:
                reference_point_in_kseg = segs_tips[kseg][0]
                measure_points_in_jseg.append(segs[jseg][np.argmin(self.Dchosen[reference_point_in_kseg, segs[jseg]])])
                reference_point_in_jseg = measure_points_in_jseg[-1]
                measure_points_in_kseg.append(segs[kseg][np.argmin(self.Dchosen[reference_point_in_jseg, segs[kseg]])])
                distances.append(self.Dchosen[measure_points_in_jseg[-1], measure_points_in_kseg[-1]])
                logg.m('   ',
                       jseg, '(tip: {}, clos: {})'.format(segs_tips[jseg][0], measure_points_in_jseg[-1]),
                       kseg, '(tip: {}, clos: {})'.format(segs_tips[kseg][0], measure_points_in_kseg[-1]),
                       '->', distances[-1], v=4)
        elif self.attachedness_measure == 'random_walk':
            for kseg in kseg_list:
                closest_distance = 1e12
                measure_point_in_jseg = 0
                measure_point_in_kseg = 0
                distances_pairs = []
                robust_quantile_jseg = int(0.0*len(segs[jseg]))
                robust_quantile_kseg = int(0.0*len(segs[kseg]))
                for reference_point_in_kseg in segs[kseg]:
                    position_in_jseg = np.argpartition(self.Dchosen[reference_point_in_kseg, segs[jseg]], robust_quantile_jseg)[robust_quantile_jseg]
                    measure_point_in_jseg_test = segs[jseg][position_in_jseg]
                    distances_pairs.append(self.Dchosen[reference_point_in_kseg, measure_point_in_jseg_test])
                    if distances_pairs[-1] < closest_distance:
                        measure_point_in_jseg = measure_point_in_jseg_test
                        measure_point_in_kseg = reference_point_in_kseg
                        closest_distance = distances_pairs[-1]
                measure_points_in_kseg.append(measure_point_in_kseg)
                measure_points_in_jseg.append(measure_point_in_jseg)
                closest_distance = np.partition(distances_pairs, robust_quantile_kseg)[robust_quantile_kseg]
                distances.append(closest_distance)
                median_distance = np.median(self.Dchosen[measure_point_in_kseg, segs[jseg]])
                median_distances.append(median_distance)
                logg.m('   ',
                       jseg, '({})'.format(measure_points_in_jseg[-1]),
                       kseg, '({})'.format(measure_points_in_kseg[-1]),
                       '->', distances[-1], median_distance, v=4)
        elif self.attachedness_measure == 'euclidian_distance_full_pairwise':
            for kseg in kseg_list:
                closest_similarity = 1e12
                measure_point_in_jseg = 0
                measure_point_in_kseg = 0
                for reference_point_in_kseg in segs[kseg]:
                    measure_point_in_jseg_test = segs[jseg][np.argmax(self.Ktilde[reference_point_in_kseg, segs[jseg]])]
                    if self.Ktilde[reference_point_in_kseg, measure_point_in_jseg_test] > closest_similarity:
                        measure_point_in_jseg = measure_point_in_jseg_test
                        measure_point_in_kseg = reference_point_in_kseg
                        closest_similarity = self.Ktilde[reference_point_in_kseg, measure_point_in_jseg_test]
                measure_points_in_kseg.append(measure_point_in_kseg)
                measure_points_in_jseg.append(measure_point_in_jseg)
                closest_distance = 1/closest_similarity
                distances.append(closest_distance)
                logg.m('   ',
                       jseg, '(tip: {}, clos: {})'.format(segs_tips[jseg][0], measure_points_in_jseg[-1]),
                       kseg, '(tip: {}, clos: {})'.format(segs_tips[kseg][0], measure_points_in_kseg[-1]),
                       '->', distances[-1], v=4)
        elif self.attachedness_measure == 'connectedness_brute_force':
            segs_jseg = set(segs[jseg])
            for kseg in kseg_list:
                connectedness = 0
                for reference_point_in_kseg in segs[kseg]:
                    for j in self.Ktilde[reference_point_in_kseg].nonzero()[1]:
                        if j in segs_jseg:
                            connectedness += 1
                distances.append(1./(connectedness+1))
            logg.m(' ', jseg, '-', kseg_list, '->', distances, v=4)
        else:
            raise ValueError('unknown attachedness measure')
        return distances, median_distances, measure_points_in_jseg, measure_points_in_kseg

    def trace_existing_connections(self, jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes, trunk):
        j_connects = segs_adjacency_nodes[jseg].copy()
        connectedness = [0, 0]
        not_trunk = 1 if trunk == 0 else 0
        kseg_trunk = set(segs[kseg_list[trunk]])
        kseg_not_trunk = set(segs[kseg_list[not_trunk]])
        for j_connect, connects in j_connects.items():
            for point_connect, seg_connect in connects:
                if seg_connect == kseg_list[trunk]:
                    # score = 0
                    # if self.Dsq[point_connect, j_connect] > 0:
                    #     score += 1. / (1 + self.Dsq[point_connect, j_connect])  # / (1 + len(segs_adjacency_nodes[jseg]))  # len(segs[jseg])
                    # if self.Dsq[j_connect, point_connect] > 0:
                    #     score += 1. / (1 + self.Dsq[j_connect, point_connect])  # / (1 + len(segs_adjacency_nodes[kseg_list[trunk if in_kseg_trunk else not_trunk]]))  # len(kseg_trunk if in_kseg_trunk else kseg_not_trunk)
                    score = 1
                    in_kseg_trunk = True if point_connect in kseg_trunk else False
                    if in_kseg_trunk:
                        connectedness[trunk] += score
                    else:
                        # elif point_connect in kseg_not_trunk:
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
                        connectedness[not_trunk] += score
        # distances = [1/c if c > 0 else np.inf for c in connectedness]
        distances = [1/(1+c) for c in connectedness]
        logg.m('    ', jseg, '-', kseg_list, '->', distances, v=5)
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
        # treat this in a different loop so we can normalize with surface of segment
        for p, q_list in segs_adjacency_nodes[kseg_loop].items():
            q_list = [q for q, jseg in q_list if jseg == kseg_test]
            for q in q_list:
                # score = 0
                # if self.Dsq[p, q] > 0: score += 1. / (1 + self.Dsq[p, q])  # / (1 + len(segs_adjacency_nodes[kseg_test]))  # len(seg_test)
                # if self.Dsq[q, p] > 0: score += 1. / (1 + self.Dsq[q, p])  # / (1 + len(segs_adjacency_nodes[kseg_loop]))  # len(seg_loop)
                score = 1
                connections += score
        # distance = 1/connections if connections > 0 else np.inf
        distance = 1/(1+connections)
        logg.m('    ', kseg_list[0], '-', kseg_list[1], '->', distance, v=5)
        return distance

    def compute_tree_from_clusters(self):
        segs = self.clusters_precomputed
        segs_tips = [[] for i in range(len(segs))]
        segs_adjacency = [[] for i in range(len(segs))]
        segs_distances = np.ones((len(segs), len(segs)))
        segs_adjacency_nodes = [{} for i in range(len(segs))]
        if not os.path.exists('./distances.npy'):
            for i in range(len(segs)):
                for j in range(i):
                    distance = self.establish_new_connections([i, j], segs, segs_adjacency_nodes)
                    segs_distances[i, j] = distance
                    segs_distances[j, i] = distance
            np.save('./distances.npy', segs_distances)
        else:
            segs_distances = np.load('./distances.npy')
        already_connected = []
        not_connected_points = np.arange(self.X.shape[0], dtype=int)
        not_connected = list(range(len(segs)))
        for count in range(self.n_splits):
            new_tips = [not_connected_points[np.argmax(self.Dchosen[not_connected_points[0], not_connected_points])]]
            dtip_others = self.Dchosen[new_tips[0], not_connected_points]
            dists = [np.max(dtip_others)]
            for j in range(10):
                new_tip = not_connected_points[np.argmax(dtip_others)]
                if new_tip in new_tips: break
                new_tips.append(new_tip)
                dtip_j = self.Dchosen[new_tips[-1], not_connected_points]
                dists.append(np.max(dtip_j))
                dtip_others += dtip_j
            tip_idx_max = np.argmax(dists)
            new_tip = new_tips.pop(tip_idx_max)
            dist_max = dists.pop(tip_idx_max)
            for iseg in not_connected:
                if new_tip in set(segs[iseg]): break
            new_seg = segs[iseg]
            pos_new_seg = np.in1d(not_connected_points, new_seg, assume_unique=True)
            not_connected_points = not_connected_points[~pos_new_seg]
            # adjust adjacency
            do_not_attach_ksegs_with_each_other = False
            continue_after_distance_compute = False
            kseg_list = [iseg]
            for kseg in kseg_list:
                jseg_list = [jseg for jseg in range(len(segs))
                             if jseg != kseg and jseg not in segs_adjacency[kseg]]  # prev_connecting_segments]  # if it's a cluster split, this is allowed?
                idcs = np.argsort(segs_distances[kseg, jseg_list])[::-1]
                for idx in idcs:
                    jseg_min = jseg_list[idx]
                    logg.m('    consider connecting', kseg, 'to', jseg_min, v=4)
                    if jseg_min not in kseg_list:
                        segs_adjacency_sparse = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
                        for i, neighbors in enumerate(segs_adjacency):
                            segs_adjacency_sparse[i, neighbors] = 1
                        G = nx.Graph(segs_adjacency_sparse)
                        paths_all = nx.single_source_dijkstra_path(G, source=kseg)
                        if jseg_min not in paths_all:
                            segs_adjacency[jseg_min].append(kseg)
                            segs_adjacency[kseg].append(jseg_min)
                            logg.m('            attaching new segment', kseg, 'at', jseg_min, v=4)
                            break
                        else:
                            logg.m('        cannot attach new segment', kseg, 'at', jseg_min,
                                      '(would produce cycle)', v=4)
        return segs, segs_tips, segs_distances, segs_adjacency

    def adjust_adjacency(self, iseg, n_add, segs, segs_tips, segs_adjacency,
                         segs_adjacency_nodes, segs_distances, trunk):
        prev_connecting_segments = segs_adjacency[iseg].copy()
        segs_adjacency += [[] for i in range(n_add)]
        segs_adjacency_nodes += [{} for i in range(n_add)]
        kseg_list = list(range(len(segs) - n_add, len(segs))) + [iseg]
        trunk = len(kseg_list) - 1
        if self.attachedness_measure == 'connectedness':
            jseg_list = [jseg for jseg in range(len(segs)) if jseg not in kseg_list]
            for jseg in jseg_list:
                distances = self.trace_existing_connections(jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes, trunk=trunk)
                segs_distances[jseg, kseg_list] = distances
                segs_distances[kseg_list, jseg] = distances
            distance = self.establish_new_connections(kseg_list, segs, segs_adjacency_nodes)
            segs_distances[kseg_list[0], kseg_list[1]] = distance
            segs_distances[kseg_list[1], kseg_list[0]] = distance
        # treat existing connections
        # logg.info('... treat existing connections')
        for jseg in prev_connecting_segments:
            median_distances = []
            if self.attachedness_measure != 'connectedness':
                result = self.compute_attachedness(jseg, kseg_list, segs, segs_tips, segs_adjacency_nodes)
                distances, median_distances, measure_points_in_jseg, measure_points_in_kseg = result
                segs_distances[jseg, kseg_list] = distances
                segs_distances[kseg_list, jseg] = distances
            distances = segs_distances[jseg, kseg_list]
            # in case we do not really have evidence for a connection based on the maximal distances
            if (median_distances
                and ((max(distances) < 0.1 and max(distances) / min(distances) < 2.5)
                     or (min(distances) >= 0.1 and max(distances) / min(distances) < 1.05))
                and min(median_distances) / max(median_distances) < 0.95):
                idx = np.argmin(median_distances)
            else:
                idx = np.argmin(distances)
            kseg_min = kseg_list[idx]
            pos = segs_adjacency[jseg].index(iseg)
            segs_adjacency[jseg][pos] = kseg_min
            pos_2 = segs_adjacency[iseg].index(jseg)
            segs_adjacency[iseg].pop(pos_2)
            segs_adjacency[kseg_min].append(jseg)
            logg.m('    segment {} is now attached to {}'.format(jseg, kseg_min), v=4)
        # in case the segment we split should correspond to two "clusters", we
        # need to check whether the new segments connect to any of the other old
        # segments
        # if not, we add a link between the new segments, if yes, we add two
        # links to connect them at the correct old segments
        # logg.info('... treat new connections')
        do_not_attach_ksegs_with_each_other = False
        continue_after_distance_compute = False
        for kseg in kseg_list:
            jseg_list = [jseg for jseg in range(len(segs))
                         if jseg != kseg and jseg not in segs_adjacency[kseg]]  # prev_connecting_segments]  # if it's a cluster split, this is allowed?
            if self.attachedness_measure != 'connectedness':
                result = self.compute_attachedness(kseg, jseg_list, segs, segs_tips, segs_adjacency_nodes)
                distances, median_distances, measure_points_in_kseg, measure_points_in_jseg = result
                segs_distances[kseg, jseg_list] = distances
                segs_distances[jseg_list, kseg] = distances
            if continue_after_distance_compute: continue
            idx = np.argmin(segs_distances[kseg, jseg_list])
            jseg_min = jseg_list[idx]
            logg.m('    consider connecting', kseg, 'to', jseg_min, v=4)
            if jseg_min not in kseg_list:
                segs_adjacency_sparse = sp.sparse.lil_matrix((len(segs), len(segs)), dtype=float)
                for i, neighbors in enumerate(segs_adjacency):
                    segs_adjacency_sparse[i, neighbors] = 1
                G = nx.Graph(segs_adjacency_sparse)
                paths_all = nx.single_source_dijkstra_path(G, source=kseg)
                if jseg_min not in paths_all:
                    segs_adjacency[jseg_min].append(kseg)
                    segs_adjacency[kseg].append(jseg_min)
                    logg.m('        attaching new segment', kseg, 'at', jseg_min, v=4)
                    # if we split the cluster, we should not attach kseg
                    do_not_attach_ksegs_with_each_other = True
                else:
                    logg.m('        cannot attach new segment', kseg, 'at', jseg_min,
                           '(would produce cycle)', v=4)
                    if kseg != kseg_list[-1]:
                        logg.m('            continue', v=4)
                        continue
                    else:
                        logg.m('            do not add another link', v=4)
                        continue_after_distance_compute = True
            if jseg_min in kseg_list and not do_not_attach_ksegs_with_each_other:
                segs_adjacency[jseg_min].append(kseg)
                segs_adjacency[kseg].append(jseg_min)
                continue_after_distance_compute = True
                logg.m('        attaching new segment', kseg, 'with new segment', jseg_min, v=4)
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
        if 'tri' == self.flavor_develop:
            ssegs = self._do_split_single_wolf17_tri(Dseg, tips)
        elif 'bi' == self.flavor_develop:
            ssegs = self._do_split_single_wolf17_bi(Dseg, tips)
        else:
            raise ValueError('unknown `self.flavor_develop`')
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
            measure_points = np.zeros((3, 3), dtype=int)
            # this is another strategy than for the undecided_cells
            # here it's possible to use the more symmetric procedure
            # shouldn't make much of a difference
            measure_points[0, 1] = ssegs[1][np.argmin(Dseg[reference_point[0]][ssegs[1]])]
            measure_points[1, 0] = ssegs[0][np.argmin(Dseg[reference_point[1]][ssegs[0]])]
            measure_points[0, 2] = ssegs[2][np.argmin(Dseg[reference_point[0]][ssegs[2]])]
            measure_points[2, 0] = ssegs[0][np.argmin(Dseg[reference_point[2]][ssegs[0]])]
            measure_points[1, 2] = ssegs[2][np.argmin(Dseg[reference_point[1]][ssegs[2]])]
            measure_points[2, 1] = ssegs[1][np.argmin(Dseg[reference_point[2]][ssegs[1]])]
            added_dist = np.zeros(3)
            added_dist[0] = Dseg[measure_points[1, 0], measure_points[0, 1]] + Dseg[measure_points[2, 0], measure_points[0, 2]]
            added_dist[1] = Dseg[measure_points[0, 1], measure_points[1, 0]] + Dseg[measure_points[2, 1], measure_points[1, 2]]
            added_dist[2] = Dseg[measure_points[1, 2], measure_points[2, 1]] + Dseg[measure_points[0, 2], measure_points[2, 0]]
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
            #             pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][measure_points[iseg_new, i], 0],
            #                        self.passed_adata.smp['X_diffmap'][seg_reference][measure_points[iseg_new, i], 1], marker='o', c='black')
            #             pl.scatter(self.passed_adata.smp['X_diffmap'][seg_reference][measure_points[i, iseg_new], 0],
            #                        self.passed_adata.smp['X_diffmap'][seg_reference][measure_points[i, iseg_new], 1], marker='x', c='black')
            #     pl.xticks([])
            #     pl.yticks([])
            #     # pl.savefig('./figs/cutting_off_tip={}.png'.format(iseg_new))
            # pl.show()
            # print('trunk', trunk)
        else:
            trunk = 0
            ssegs_adjacency = [[1], [0]]
            reference_point_in_0 = ssegs_tips[0][0]
            measure_point_in_1 = ssegs[1][np.argmin(Dseg[reference_point_in_0][ssegs[1]])]
            reference_point_in_1 = measure_point_in_1  # ssegs_tips[1][0]
            measure_point_in_0 = ssegs[0][np.argmin(Dseg[reference_point_in_1][ssegs[0]])]
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

# Author: Alex Wolf (http://falexwolf.de)

from collections import namedtuple
import numpy as np
import scipy as sp
import pandas as pd
import networkx as nx
import scipy.sparse
from textwrap import indent, dedent
from natsort import natsorted
from .. import logging as logg
from ..data_structs import data_graph
from .. import utils
from .. import settings
from ..plotting import utils as pl_utils


MINIMAL_TREE_ATTACHEDNESS = 0.05

doc_string_base = dedent("""\
    Generate cellular maps of differentiation manifolds with complex
    topologies [Wolf17i]_.

    Approximate graph abstraction (AGA) quantifies the connectivity of
    partitions of a neighborhood graph of single cells, thereby generating a
    much simpler abstracted graph whose nodes label the partitions. Together
    with a random walk-based distance measure, this generates a partial
    coordinatization of data useful for exploring and explaining its
    variation. By default, AGA uses the Louvain algorithm to partition the data,
    which has been suggested for clustering single-cell data by
    [Levine15]_. Also, it extends the random-walk based distance measure
    suggested by [Haghverdi16]_.

    **Note**: In order to compute distances along the graph (pseudotimes), you need
    to provide a root cell, e.g., as in the `example of Nestorowa et al. (2016)
    <https://github.com/theislab/graph_abstraction/blob/master/nestorowa16/nestorowa16.ipynb>`__::

        adata.uns['iroot'] = np.flatnonzero(adata.smp['exp_groups'] == 'Stem')[0]

    You should get good results with the default parameters. Most the parameters
    appear similarly in other tools and are used to generate the graph.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix, optionally with `adata.uns['iroot']`, the index
        of root cell for computing a pseudotime.
    groups : categorical smp/cell annotation or {{'louvain_groups', 'segments'}}, optional (default: 'louvain_groups')
        Criterion to determine the resulting partitions of the single-cell
        graph. 'louvain_groups' uses the louvain algorithm and optimizes
        modularity of the graph, 'segments' uses a bipartioning criterium that
        similar to hierarchical clustering on the graph. You can also pass your
        predefined groups by choosing any sample annotation.
    n_pcs : `int`, optional (default: 50)
        Use n_pcs PCs to compute the euclidean distance matrix, which is the
        basis for generating the graph. Set to 0 if you don't want preprocessing
        with PCA.
    n_neighbors : `int` or `None`, optional (default: `None`)
        Number of nearest neighbors on the knn graph. Often this can be reduced
        down to a value of 4. Defaults to the number of neighbors in a
        precomputed graph. If there is none, defaults to 30.
    n_dcs : `int`, optional (default: 10)
        Number of diffusion components (very similar to eigen vectors of
        adjacency matrix) to use for distance computations.
    resolution : `float`, optional (default: 1.0)
        See :func:`~scanpy.api.louvain`.
    random_state : `int`, optional (default: 0)
        See :func:`~scanpy.api.louvain`.
    tree_detection : {{'iterative_matching', 'min_span_tree'}}, optional (default: 'min_span_tree')
        How to detect a tree structure in the abstracted graph. If choosing
        'min_span_tree', he minimum spanning tree is computed for the abstracted
        graph with inverted weights. If choosing 'iterative_matching', this runs
        a recursive algorithm that greedily attaches partitions (groups) that
        maximize the random-walk based distance measure.
    attachedness_measure : {{'connectedness', 'random_walk'}}, optional (default: 'connectedness')
        How to measure connectedness between groups.
    n_nodes : `int` or `None`, optional (default: `None`)
        Number of nodes in the abstracted graph. Except when choosing
        'segments' for `groups`, for which `n_nodes` defaults to
        `n_nodes=1`, `n_nodes` defaults to the number of groups implied by the
        choice of `groups`.
    recompute_graph : `bool`, optional (default: `False`)
        Recompute single-cell graph.
    recompute_pca : `bool`, optional (default: `False`)
        Recompute PCA.
    recompute_louvain : `bool`, optional (default: `False`)
        When changing the `resolution` parameter, you should set this to True.
    n_jobs : `int` or None (default: `sc.settings.n_jobs`)
        Number of CPUs to use for parallel processing.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy. Otherwise, perform
        computation inplace and return None.

    Returns
    -------
    Returns or updates `adata` depending on `copy` with
    {returns}
    """)


doc_string_returns = dedent("""\
        aga_adjacency_full_attachedness : np.ndarray (adata.uns)
            The full adjacency matrix of the abstracted graph, weights
            correspond to connectedness.
        aga_adjacency_full_confidence : np.ndarray (adata.uns)
            The full adjacency matrix of the abstracted graph, weights
            correspond to confidence in the presence of an edge.
        aga_adjacency_tree_confidence : sc.sparse csr matrix (adata.uns)
            The adjacency matrix of the tree-like subgraph that best explains
            the topology.
        aga_pseudotime : pd.Series (adata.smp, dtype float)
            Pseudotime labels, that is, distance a long the manifold for each
            cell.
    """)


def aga(adata,
        groups='louvain_groups',
        n_pcs=50,
        n_neighbors=None,
        n_dcs=None,
        resolution=None,
        random_state=0,
        attachedness_measure='connectedness',
        tree_detection='min_span_tree',
        tree_based_confidence=True,
        n_nodes=None,
        recompute_pca=False,
        recompute_distances=False,
        recompute_graph=False,
        recompute_louvain=False,
        n_jobs=None,
        copy=False):
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    if tree_detection not in {'iterative_matching', 'min_span_tree'}:
        raise ValueError('`tree_detection` needs to be one of {}'
                         .format({'iterative_matching', 'min_span_tree'}))
    fresh_compute_louvain = False
    if (groups == 'louvain_groups'
        and ('louvain_groups' not in adata.smp_keys()
             # resolution does not match
             or ('louvain_params' in adata.uns
                 and resolution is not None
                 and adata.uns['louvain_params']['resolution'] != resolution)
             or recompute_louvain
             or not data_graph.no_recompute_of_graph_necessary(
            adata,
            recompute_pca=recompute_pca,
            recompute_distances=recompute_distances,
            recompute_graph=recompute_graph,
            n_neighbors=n_neighbors,
            n_dcs=n_dcs))):
        from .louvain import louvain
        louvain(adata,
                resolution=resolution,
                n_neighbors=n_neighbors,
                recompute_pca=recompute_pca,
                recompute_graph=recompute_graph,
                n_pcs=n_pcs,
                n_dcs=n_dcs,
                random_state=random_state)
        fresh_compute_louvain = True
    clusters = groups
    logg.info('running Approximate Graph Abstraction (AGA)', reset=True)
    if ('iroot' not in adata.uns
        and 'xroot' not in adata.uns
        and 'xroot' not in adata.var):
        logg.info('    no root cell found, no computation of pseudotime')
        msg = \
    """To enable computation of pseudotime, pass the index or expression vector
    of a root cell. Either add
        adata.uns['iroot'] = root_cell_index
    or (robust to subsampling)
        adata.var['xroot'] = adata.X[root_cell_index, :]
    where "root_cell_index" is the integer index of the root cell, or
        adata.var['xroot'] = adata[root_cell_name, :].X
    where "root_cell_name" is the name (a string) of the root cell."""
        logg.hint(msg)
    aga = AGA(adata,
              clusters=clusters,
              n_neighbors=n_neighbors,
              n_pcs=n_pcs,
              n_dcs=n_dcs,
              n_jobs=n_jobs,
              tree_based_confidence=tree_based_confidence,
              # we do not need to recompute things both in the louvain
              # call above and here
              recompute_graph=recompute_graph and not fresh_compute_louvain,
              recompute_distances=recompute_distances and not fresh_compute_louvain,
              recompute_pca=recompute_pca and not fresh_compute_louvain,
              n_nodes=n_nodes,
              attachedness_measure=attachedness_measure)
    updated_diffmap = aga.update_diffmap()
    adata.smpm['X_diffmap'] = aga.rbasis[:, 1:]
    adata.smp['X_diffmap0'] = aga.rbasis[:, 0]
    adata.uns['diffmap_evals'] = aga.evals[1:]
    adata.uns['data_graph_distance_local'] = aga.Dsq
    adata.uns['data_graph_norm_weights'] = aga.Ktilde
    if aga.iroot is not None:
        aga.set_pseudotime()  # pseudotimes are random walk distances from root point
        adata.uns['iroot'] = aga.iroot  # update iroot, might have changed when subsampling, for example
        adata.smp['aga_pseudotime'] = aga.pseudotime
    # detect splits and partition the data into segments
    aga.splits_segments()

    if tree_detection == 'min_span_tree':
        min_span_tree = utils.compute_minimum_spanning_tree(
            1./aga.segs_adjacency_full_attachedness)
        min_span_tree.data = 1./min_span_tree.data
        full_confidence, tree_confidence = aga.compute_adjacency_confidence(
            aga.segs_adjacency_full_attachedness, min_span_tree, tree_based_confidence)
    else:
        full_confidence, tree_confidence = aga.segs_adjacency_full_confidence, aga.segs_adjacency_tree_confidence

    y = adata.smp[clusters].cat.categories
    x = np.array(aga.segs_names_original)
    xsorted = np.argsort(x)
    ypos = np.searchsorted(x[xsorted], y)
    indices = xsorted[ypos]

    adata.uns['aga_adjacency_full_attachedness'] = aga.segs_adjacency_full_attachedness[indices, :][:, indices]
    adata.uns['aga_adjacency_full_confidence'] = full_confidence[indices, :][:, indices]
    adata.uns['aga_adjacency_tree_confidence'] = tree_confidence[indices, :][:, indices]
    adata.uns['aga_groups_key'] = clusters
    adata.uns[clusters + '_sizes'] = np.array(aga.segs_sizes)[indices]
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n' + indent(doc_string_returns, '    '))
    return adata if copy else None

aga.__doc__ = doc_string_base.format(returns=doc_string_returns)


def aga_degrees(adata):
    """Compute the degree of each node in the abstracted graph.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    degrees : list
        List of degrees for each node.
    """
    import networkx as nx
    g = nx.Graph(adata.uns['aga_adjacency_full_confidence'])
    degrees = [d for _, d in g.degree(weight='weight')]
    return degrees


def aga_expression_entropies(adata):
    """Compute the median expression entropy for each node-group.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    entropies : list
        Entropies of median expressions for each node.
    """
    from scipy.stats import entropy
    groups_order, groups_masks = utils.select_groups(adata,
                                                     key=adata.uns['aga_groups_key'])
    entropies = []
    for mask in groups_masks:
        X_mask = adata.X[mask]
        x_median = np.median(X_mask, axis=0)
        x_probs = (x_median - np.min(x_median)) / (np.max(x_median) - np.min(x_median))
        entropies.append(entropy(x_probs))
    return entropies


def aga_compare_paths(adata1, adata2,
                      adjacency_key='aga_adjacency_full_confidence'):
    """Compare paths in abstracted graphs in two datasets.

    Compute the fraction of consistent paths between leafs, a measure for the
    topological similarity between graphs.

    By increasing the verbosity to level 4 and 5, the paths that do not agree
    and the paths that agree are written to the output, respectively.

    Parameters
    ----------
    adata1, adata2 : AnnData
        Annotated data matrices to compare.
    adjacency_key : str
        Key for indexing the adjacency matrices to be used in adata1 and adata2.

    Returns
    -------
    OrderedTuple with attributes ``n_steps`` (total number of steps in paths)
    and ``frac_steps`` (fraction of consistent steps), ``n_paths`` and
    ``frac_paths``.
    """
    import networkx as nx
    g1 = nx.Graph(adata1.uns[adjacency_key])
    g2 = nx.Graph(adata2.uns[adjacency_key])
    leaf_nodes1 = [str(x) for x in g1.nodes() if g1.degree(x) == 1]
    logg.msg('leaf nodes in graph 1: {}'.format(leaf_nodes1), v=5, no_indent=True)
    asso_groups1 = utils.identify_groups(adata1.smp['aga_groups'].values,
                                         adata2.smp['aga_groups'].values)
    asso_groups2 = utils.identify_groups(adata2.smp['aga_groups'].values,
                                         adata1.smp['aga_groups'].values)
    orig_names1 = adata1.uns['aga_groups_order_original']
    orig_names2 = adata2.uns['aga_groups_order_original']

    import itertools
    n_steps = 0
    n_agreeing_steps = 0
    n_paths = 0
    n_agreeing_paths = 0
    # loop over all pairs of leaf nodes in the reference adata1
    for (r, s) in itertools.combinations(leaf_nodes1, r=2):
        r2, s2 = asso_groups1[r][0], asso_groups1[s][0]
        orig_names = [orig_names1[int(i)] for i in [r, s]]
        orig_names += [orig_names2[int(i)] for i in [r2, s2]]
        logg.msg('compare shortest paths between leafs ({}, {}) in graph1 and ({}, {}) in graph2:'
               .format(*orig_names), v=4, no_indent=True)
        no_path1 = False
        try:
            path1 = [str(x) for x in nx.shortest_path(g1, int(r), int(s))]
        except nx.NetworkXNoPath:
            no_path1 = True
        no_path2 = False
        try:
            path2 = [str(x) for x in nx.shortest_path(g2, int(r2), int(s2))]
        except nx.NetworkXNoPath:
            no_path2 = True
        if no_path1 and no_path2:
            # consistent behavior
            n_paths += 1
            n_agreeing_paths += 1
            n_steps += 1
            n_agreeing_steps += 1
            continue
        elif no_path1 or no_path2:
            # non-consistent result
            n_paths += 1
            n_steps += 1
            continue
        if len(path1) >= len(path2):
            path_mapped = [asso_groups1[l] for l in path1]
            path_compare = path2
            path_compare_id = 2
            path_compare_orig_names = [[orig_names2[int(s)] for s in l] for l in path_compare]
            path_mapped_orig_names = [[orig_names2[int(s)] for s in l] for l in path_mapped]
        else:
            path_mapped = [asso_groups2[l] for l in path2]
            path_compare = path1
            path_compare_id = 1
            path_compare_orig_names = [[orig_names1[int(s)] for s in l] for l in path_compare]
            path_mapped_orig_names = [[orig_names1[int(s)] for s in l] for l in path_mapped]
        n_agreeing_steps_path = 0
        ip_progress = 0
        for il, l in enumerate(path_compare[:-1]):
            for ip, p in enumerate(path_mapped):
                if ip >= ip_progress and l in p:
                    # check whether we can find the step forward of path_compare in path_mapped
                    if (ip + 1 < len(path_mapped)
                        and
                        path_compare[il + 1] in path_mapped[ip + 1]):
                        # make sure that a step backward leads us to the same value of l
                        # in case we "jumped"
                        logg.msg('found matching step ({} -> {}) at position {} in path{} and position {} in path_mapped'
                               .format(l, path_compare_orig_names[il + 1], il, path_compare_id, ip), v=6)
                        consistent_history = True
                        for iip in range(ip, ip_progress, -1):
                            if l not in path_mapped[iip - 1]:
                                consistent_history = False
                        if consistent_history:
                            # here, we take one step further back (ip_progress - 1); it's implied that this
                            # was ok in the previous step
                            logg.msg('    step(s) backward to position(s) {} in path_mapped are fine, too: valid step'
                                   .format(list(range(ip - 1, ip_progress - 2, -1))), v=6)
                            n_agreeing_steps_path += 1
                            ip_progress = ip + 1
                            break
        n_steps_path = len(path_compare) - 1
        n_agreeing_steps += n_agreeing_steps_path
        n_steps += n_steps_path
        n_paths += 1
        if n_agreeing_steps_path == n_steps_path: n_agreeing_paths += 1

        # only for the output, use original names
        path1_orig_names = [orig_names1[int(s)] for s in path1]
        path2_orig_names = [orig_names2[int(s)] for s in path2]
        logg.msg('      path1 = {},\n'
               'path_mapped = {},\n'
               '      path2 = {},\n'
               '-> n_agreeing_steps = {} / n_steps = {}.'
               .format(path1_orig_names,
                       [list(p) for p in path_mapped_orig_names],
                       path2_orig_names,
                       n_agreeing_steps_path, n_steps_path), v=5, no_indent=True)
    Result = namedtuple('aga_compare_paths_result',
                        ['frac_steps', 'n_steps', 'frac_paths', 'n_paths'])
    return Result(frac_steps=n_agreeing_steps/n_steps if n_steps > 0 else np.nan,
                  n_steps=n_steps if n_steps > 0 else np.nan,
                  frac_paths=n_agreeing_paths/n_paths if n_steps > 0 else np.nan,
                  n_paths=n_paths if n_steps > 0 else np.nan)


def aga_contract_graph(adata, min_group_size=0.01, max_n_contractions=1000, copy=False):
    """Contract the abstracted graph.
    """
    adata = adata.copy() if copy else adata
    if 'aga_adjacency_tree_confidence' not in adata.uns: raise ValueError('run tool aga first!')
    min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * adata.n_smps)
    logg.info('contract graph using `min_group_size={}`'.format(min_group_size))

    def propose_nodes_to_contract(adjacency_tree_confidence, node_groups):
        # nodes with two edges
        n_edges_per_seg = np.sum(adjacency_tree_confidence > 0, axis=1).A1
        for i in range(adjacency_tree_confidence.shape[0]):
            if n_edges_per_seg[i] == 2:
                neighbors = adjacency_tree_confidence[i].nonzero()[1]
                for neighbors_edges in range(1, 20):
                    for n_cnt, n in enumerate(neighbors):
                        if n_edges_per_seg[n] == neighbors_edges:
                            logg.msg('merging node {} into {} (two edges)'
                                   .format(i, n), v=4)
                            return i, n
        # node groups with a very small cell number
        for i in range(adjacency_tree_confidence.shape[0]):
            if node_groups[str(i) == node_groups].size < min_group_size:
                neighbors = adjacency_tree_confidence[i].nonzero()[1]
                neighbor_sizes = [node_groups[str(n) == node_groups].size for n in neighbors]
                n = neighbors[np.argmax(neighbor_sizes)]
                logg.msg('merging node {} into {} '
                       '(smaller than `min_group_size` = {})'
                       .format(i, n, min_group_size), v=4)
                return i, n
        return 0, 0

    def contract_nodes(adjacency_tree_confidence, node_groups):
        for count in range(max_n_contractions):
            i, n = propose_nodes_to_contract(adjacency_tree_confidence, node_groups)
            if i != 0 or n != 0:
                G = nx.Graph(adjacency_tree_confidence)
                G_contracted = nx.contracted_nodes(G, n, i, self_loops=False)
                adjacency_tree_confidence = nx.to_scipy_sparse_matrix(G_contracted)
                node_groups[str(i) == node_groups] = str(n)
                for j in range(i+1, G.size()+1):
                    node_groups[str(j) == node_groups] = str(j-1)
            else:
                break
        return adjacency_tree_confidence, node_groups

    size_before = adata.uns['aga_adjacency_tree_confidence'].shape[0]
    adata.uns['aga_adjacency_tree_confidence'], adata.smp['aga_groups'] = contract_nodes(
        adata.uns['aga_adjacency_tree_confidence'], adata.smp['aga_groups'].values)
    adata.uns['aga_groups_order'] = np.unique(adata.smp['aga_groups'].values)
    for key in ['aga_adjacency_full_confidence', 'aga_groups_original',
                'aga_groups_order_original', 'aga_groups_colors_original']:
        if key in adata.uns: del adata.uns[key]
    logg.info('    contracted graph from {} to {} nodes'
              .format(size_before, adata.uns['aga_adjacency_tree_confidence'].shape[0]))
    logg.msg('removed adata.uns["aga_adjacency_full_confidence"]', v=4)
    return adata if copy else None


class AGA(data_graph.DataGraph):
    """Approximate Graph Abstraction

    This needs to be rewritten in a cleaner way.
    """

    def __init__(self,
                 adata,
                 n_nodes=None,
                 n_neighbors=30,
                 n_pcs=50,
                 n_dcs=10,
                 min_group_size=1,
                 tree_based_confidence=True,
                 minimal_distance_evidence=0.95,
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
        self.minimal_distance_evidence = minimal_distance_evidence
        # the ratio of max(minimal_distances)/min(minimal_distances) has to be smaller than minimal_distance_evidence
        # in order to be considered convincing evidence, otherwise, consider median_distances
        self.min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * self.X.shape[0])
        self.passed_adata = adata  # just for debugging purposes
        self.choose_largest_segment = True
        self.attachedness_measure = attachedness_measure
        self.tree_based_confidence = tree_based_confidence
        self.clusters = clusters
        self.clusters_precomputed = None
        self.clusters_precomputed_names = None
        self.flavor_develop = 'bi'  # bipartitioning
        if clusters not in {'segments', 'unconstrained_segments'}:
            if clusters not in adata.smp_keys():
                raise ValueError('Did not find {} in adata.smp_keys()! '
                                 'If you do not have any precomputed clusters, pass "segments" for "node_groups" instead'
                                 .format(clusters))
            clusters_array = adata.smp[clusters].values
            # transform to a list of index arrays
            self.clusters_precomputed = []
            self.clusters_precomputed_names = list(adata.smp[clusters].cat.categories)
            for cluster_name in self.clusters_precomputed_names:
                self.clusters_precomputed.append(np.where(cluster_name == clusters_array)[0])
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
            Array of dimension (number of segments) × (number of data
            points). Each row stores a mask array that defines a segment.
        segs_tips : np.ndarray
            Array of dimension (number of segments) × 2. Each row stores the
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
                logg.msg('    split', ibranch + 1, v=4)
                stop, segs_distances = self.do_split_constrained(segs, segs_tips,
                                                                 segs_adjacency,
                                                                 segs_adjacency_nodes,
                                                                 segs_distances)
                if stop: break

        # segments
        self.segs = segs
        self.segs_tips = segs_tips
        self.segs_sizes = []
        for iseg, seg in enumerate(self.segs): self.segs_sizes.append(len(seg))

        # the full, unscaled adjacency matrix
        self.segs_adjacency_full_attachedness = 1/segs_distances
        # if self.attachedness_measure == 'connectedness':
        #     norm = np.sqrt(np.multiply.outer(self.segs_sizes, self.segs_sizes))
        #     self.segs_adjacency_full_attachedness /= norm
        self.segs_adjacency_full_confidence, self.segs_adjacency_tree_confidence \
            = self.compute_adjacency_confidence(
                self.segs_adjacency_full_attachedness,
                segs_adjacency,
                self.tree_based_confidence)
        np.fill_diagonal(self.segs_adjacency_full_attachedness, 0)

    def compute_adjacency_confidence(self, full_attachedness, tree_adjacency, tree_based_confidence):
        """Translates the attachedness measure into a confidence measure.
        """
        if sp.sparse.issparse(tree_adjacency):
            tree_adjacency = [tree_adjacency[i].nonzero()[1]
                              for i in range(tree_adjacency.shape[0])]
        segs_distances = 1/full_attachedness
        if not tree_based_confidence:  # inter- and intra-cluster based confidence
            from scipy.stats import norm
            # intra-cluster connections
            total_n = self.k * np.array(self.segs_sizes)  # total number of connections
            a = full_attachedness
            confidence = np.zeros_like(full_attachedness)
            logg.msg('computing confidence', v=5)
            logg.msg('i_name, j_name, connectivity, total_n[i], total_n[j], '
                     'actual, expected, variance, confidence', v=5)
            for i in range(a.shape[0]):
                for j in range(i+1, a.shape[1]):
                    expected = total_n[i] * total_n[j] / np.sum(total_n)**2
                    actual = a[i, j] / np.sum(total_n)
                    variance = expected * (1 - expected) / np.sum(total_n)
                    if actual > expected:
                        confidence[i, j] = 1
                    elif actual < 1e-12:
                        confidence[i, j] = 0
                    else:
                        confidence[i, j] = 2 * norm.cdf(
                            actual, expected, np.sqrt(variance))
                    i_name = self.segs_names_original[i]
                    j_name = self.segs_names_original[j]
                    logg.msg(i_name, j_name, a[i, j], total_n[i], total_n[j],
                             actual, expected, variance, confidence[i, j], v=5)
            full_confidence = confidence + confidence.T
            tree_confidence = self.compute_tree_confidence(
                full_confidence, tree_adjacency)
        else:
            # compute the average tree distances
            tree_distances = []
            for i, neighbors in enumerate(tree_adjacency):
                tree_distances += segs_distances[i][neighbors].tolist()
            median_tree_distances = np.median(tree_distances)
            full_confidence = np.zeros_like(segs_distances)
            full_confidence[segs_distances <= median_tree_distances] = 1
            full_confidence[segs_distances > median_tree_distances] = (
                np.exp(-(segs_distances-median_tree_distances)/median_tree_distances)
                [segs_distances > median_tree_distances])
            np.fill_diagonal(full_confidence, 0)
            tree_confidence = self.compute_tree_confidence(
                full_confidence, tree_adjacency,
                minimal_tree_attachedness=MINIMAL_TREE_ATTACHEDNESS)
        return full_confidence, tree_confidence

    def compute_tree_confidence(self, full_confidence, tree_adjacency, minimal_tree_attachedness=1e-14):
        n = full_confidence.shape[0]
        tree_confidence = sp.sparse.lil_matrix((n, n), dtype=float)
        for i, neighbors in enumerate(tree_adjacency):
            clipped_attachedness = full_confidence[i][neighbors]
            clipped_attachedness[clipped_attachedness < minimal_tree_attachedness] = minimal_tree_attachedness
            tree_confidence[i, neighbors] = clipped_attachedness
            full_confidence[i, neighbors] = clipped_attachedness
        tree_confidence = tree_confidence.tocsr()
        return tree_confidence

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
                logg.msg('    splitting group {} with size {}'.format(iseg, len(seg)), v=4)
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
            logg.msg('    new tip {} with distance {:.6}, constraint was {}'
                   .format(new_itip, dist_new_itip, itip), v=4)
            logg.msg('    new sizes {} and {}'
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
            logg.msg('... splitting group {} with size {}'.format(iseg, len(seg)), v=4)
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
            logg.msg('    new tip {} with distance {:.6}, constraint was {}'
                   .format(new_itip, 0.0, itip), v=4)
            logg.msg('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            logg.msg('    the scores where', scores, v=4)
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
            logg.msg('    new tip', new_tip, 'with distance', dist_max,
                   'using constraints {} with distances'
                   .format(new_tips), v=4)
            logg.msg('   ', dists, v=4)
            logg.msg('    new sizes {} and {}'
                   .format(sizes[0], sizes[1]), v=4)
            return iseg, seg, ssegs, ssegs_tips, sizes

        def select_precomputed(segs_tips):
            if len(segs) == 1:
                segs_tips.pop(0)
                segs_tips.append([])
            iseg = 0
            seg = segs[iseg]
            logg.msg('    splitting group {} with size {}'.format(iseg, len(seg)), v=4)
            new_tips = [seg[np.argmax(self.Dchosen[seg[0], seg])]]
            dtip_others = self.Dchosen[new_tips[0], seg]
            dists = [np.max(dtip_others)]
            # it would be equivalent to just consider one pair of points
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
            logg.msg('    new tip', new_tip, 'with distance', dist_max,
                   'using constraints {} with distances'
                   .format(new_tips), v=4)
            logg.msg('   ', dists, v=4)
            logg.msg('    new sizes {} and {}'
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
            # there is one partition that corresponds to all other partitions...
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
            logg.msg('... group', iseg, 'score', score, 'n_points', len(seg),
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
                logg.msg('   ',
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
                logg.msg('   ',
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
                logg.msg('   ',
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
                # distances.append(1./(connectedness+1))
                distances.append(1./connectedness if connectedness != 0 else np.inf)
            logg.msg(' ', jseg, '-', kseg_list, '->', distances, v=4)
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
        distances = [1/c if c > 0 else np.inf for c in connectedness]
        # distances = [1/(1+c) for c in connectedness]
        logg.msg('    ', jseg, '-', kseg_list, '->', distances, v=5)
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
        # distance = 1/(1+connections)
        distance = 1/connections if connections > 0 else np.inf
        logg.msg('    ', kseg_list[0], '-', kseg_list[1], '->', distance, v=5)
        return distance

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
            # in case we do not have convincing evidence for a connection based on the maximal distances
            if (median_distances
                and ((max(distances) < 0.1 and min(distances) / max(distances) >= 0.4)
                     # all distances are very small, we require significant statistical evidence here
                     or (min(distances) >= 0.1 and min(distances) / max(distances) >= self.minimal_distance_evidence))
                     # distances are larger
                and min(median_distances) / max(median_distances) < self.minimal_distance_evidence):
                     # require median_distances to actually provide better evidence
                logg.msg('        no convincing evidence in minimal distances, consider median distance', v=4)
                idx = np.argmin(median_distances)
            else:
                idx = np.argmin(distances)
            kseg_min = kseg_list[idx]
            pos = segs_adjacency[jseg].index(iseg)
            segs_adjacency[jseg][pos] = kseg_min
            pos_2 = segs_adjacency[iseg].index(jseg)
            segs_adjacency[iseg].pop(pos_2)
            segs_adjacency[kseg_min].append(jseg)
            logg.msg('    group {} is now attached to {}'.format(jseg, kseg_min), v=4)
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
            # candidate for the segment to which we attach would attach the new
            # segment
            jseg_min = jseg_list[idx]
            logg.msg('    consider connecting', kseg, 'to', jseg_min, v=4)
            # if the closest segment is not among the two new segments
            if jseg_min not in kseg_list:
                segs_adjacency_sparse = sp.sparse.lil_matrix(
                    (len(segs), len(segs)), dtype=float)
                for i, neighbors in enumerate(segs_adjacency):
                    segs_adjacency_sparse[i, neighbors] = 1
                G = nx.Graph(segs_adjacency_sparse)
                paths_all = nx.single_source_dijkstra_path(G, source=kseg)
                # we can attach the new segment to an old segment
                if jseg_min not in paths_all:
                    segs_adjacency[jseg_min].append(kseg)
                    segs_adjacency[kseg].append(jseg_min)
                    logg.msg('        attaching new segment',
                           kseg, 'at', jseg_min, v=4)
                    # if we establish the new connection with an old segment
                    # we should not add a new connection to the second new segment
                    do_not_attach_ksegs_with_each_other = True
                # we cannot attach it to an old segment as this
                # would produce a cycle
                else:
                    logg.msg('        cannot attach new segment',
                           kseg, 'at', jseg_min,
                           '(would produce cycle)', v=4)
                    # we still have the other new segment to inspect so it's not
                    # a drama that we couldn't establish a new connection
                    if kseg != kseg_list[-1]:
                        logg.msg('            continue', v=4)
                        continue
                    # we do not add add a new link
                    else:
                        logg.msg('            do not add another link', v=4)
                        continue_after_distance_compute = True
            if jseg_min in kseg_list and not do_not_attach_ksegs_with_each_other:
                segs_adjacency[jseg_min].append(kseg)
                segs_adjacency[kseg].append(jseg_min)
                # we're already done as we found the new connection
                continue_after_distance_compute = True
                logg.msg('        attaching new segment',
                         kseg, 'with new segment', jseg_min, v=4)
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
            logg.msg('... is root itself, never obtain significant correlation', v=4)
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

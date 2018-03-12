from collections import namedtuple
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import minimum_spanning_tree
from textwrap import dedent
from .. import logging as logg
from ..neighbors import Neighbors
from .. import utils
from .. import settings


doc_string_base = dedent("""\
    Generate cellular maps of differentiation manifolds with complex
    topologies [Wolf17i]_.

    Statistical graph abstraction (SGA) quantifies the connectivity of
    partitions of a neighborhood graph of single cells, thereby generating a
    much simpler abstracted graph whose nodes label the partitions. Together
    with a random walk-based distance measure, this generates a partial
    coordinatization of data useful for exploring and explaining its
    variation.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    groups : categorical annotation of observations or 'louvain_groups', optional (default: 'louvain_groups')
        Criterion to determine the resulting partitions of the single-cell
        graph. 'louvain_groups' uses the Louvain algorithm and optimizes
        modularity of the graph. You can also pass your predefined groups by
        choosing any categorical annotation of observations (`adata.obs`).
    tree_based_confidence : `bool`, optional (default: `True`)
        Have high confidence in a connection if its connectivity is
        significantly higher than the median connectivity of the global spanning
        tree of the abstracted graph.
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
        sga_connectivity : np.ndarray (adata.uns)
            The full adjacency matrix of the abstracted graph, weights
            correspond to connectivity.
        sga_confidence : np.ndarray (adata.uns)
            The full adjacency matrix of the abstracted graph, weights
            correspond to confidence in the presence of an edge.
        sga_confidence_tree : sc.sparse csr matrix (adata.uns)
            The adjacency matrix of the tree-like subgraph that best explains
            the topology.
    """)


def sga(adata,
        groups='louvain_groups',
        threshold=0.01,
        tree_based_confidence=False,
        n_jobs=None,
        copy=False):
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    logg.info('running Statistical Graph Abstraction (SGA)', reset=True)
    sga = SGA(adata, groups, threshold=threshold, tree_based_confidence=tree_based_confidence)
    sga.compute()
    adata.uns['sga_connectivity'] = sga.connectivity
    adata.uns['sga_confidence'] = sga.confidence
    adata.uns['sga_confidence_tree'] = sga.confidence_tree
    adata.uns['sga_groups_key'] = groups
    adata.uns[groups + '_sizes'] = np.array(sga.vc.sizes())
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added\n'
        '    \'sga_connectivity\', connectivity adjacency (adata.uns)\n'
        '    \'sga_confidence\', confidence adjacency (adata.uns)\n'
        '    \'sga_confidence_tree\', confidence subtree (adata.uns)')
    return adata if copy else None

sga.__doc__ = doc_string_base.format(returns=doc_string_returns)


class SGA(Neighbors):
    """Statistical Graph Abstraction.
    """

    def __init__(self, adata, groups, threshold=0.05, tree_based_confidence=True):
        super(SGA, self).__init__(adata)
        self._groups = groups
        self.threshold = threshold
        self._tree_based_confidence = tree_based_confidence

    def compute(self):
        self.compute_connectivity()
        self.compute_confidence()

    def compute_connectivity(self):
        import igraph
        ones = self._adata.uns['neighbors_similarities'].copy()
        # graph where edges carry weight 1
        ones.data = np.ones(len(ones.data))
        g = utils.get_igraph_from_adjacency(ones)
        self.vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups].cat.codes.values)
        cg = self.vc.cluster_graph(combine_edges='sum')
        self.connectivity = utils.get_sparse_from_igraph(cg, weight_attr='weight')/2

    def compute_confidence(self):
        """Translates the connectivity measure into a confidence measure.
        """
        pseudo_distance = self.connectivity.copy()
        pseudo_distance.data = 1./pseudo_distance.data
        connectivity_tree = minimum_spanning_tree(pseudo_distance)
        connectivity_tree.data = 1./connectivity_tree.data
        connectivity_tree_indices = [
            connectivity_tree[i].nonzero()[1]
            for i in range(connectivity_tree.shape[0])]
        # inter- and intra-cluster based confidence
        if not self._tree_based_confidence:
            total_n = self.n_neighbors * np.array(self.vc.sizes())
            logg.msg('{:>2} {:>2} {:>4} {:>4} {:>4} '
                     '{:>7} {:>7} {:>7} {:>7}'
                     .format('i', 'j', 'conn', 'n[i]', 'n[j]',
                             'avg', 'thresh', 'var', 'conf'), v=5)
            maximum = self.connectivity.max()
            confidence = self.connectivity.copy()  # initializing
            for i in range(self.connectivity.shape[0]):
                for j in range(i+1, self.connectivity.shape[1]):
                    if self.connectivity[i, j] > 0:
                        minimum = min(total_n[i], total_n[j])
                        average = self.connectivity[i, j] / minimum
                        confidence[i, j] = self.connectivity[i, j] / maximum
                        variance = self.threshold * (1-self.threshold)
                        # if average > self.threshold:
                        #     confidence[i, j] = 1
                        # else:
                        #     confidence[i, j] = norm.cdf(average,
                        #         self.threshold, variance)
                        logg.msg(
                            '{:2} {:2} {:4} {:4} {:4} '
                            '{:7.2} {:7.2} {:7.2} {:7.2}'
                            .format(i, j, int(self.connectivity[i, j]),
                                    total_n[i], total_n[j],
                                    average, self.threshold, variance, confidence[i, j]), v=5)
                        confidence[j, i] = confidence[i, j]
        # tree-based confidence
        else:
            median_connectivity_tree = np.median(connectivity_tree.data)
            confidence = self.connectivity.copy()
            confidence.data[self.connectivity.data >= median_connectivity_tree] = 1
            connectivity_adjusted = self.connectivity.copy()
            connectivity_adjusted.data -= median_connectivity_tree
            connectivity_adjusted.data = np.exp(connectivity_adjusted.data)
            index = self.connectivity.data < median_connectivity_tree
            confidence.data[index] = connectivity_adjusted.data[index]
        confidence_tree = self.compute_confidence_tree(
            confidence, connectivity_tree_indices)
        self.confidence = confidence
        self.confidence_tree = confidence_tree

    def compute_confidence_tree(
            self, confidence, connectivity_tree_indices):
        confidence_tree = sp.sparse.lil_matrix(confidence.shape, dtype=float)
        for i, neighbors in enumerate(connectivity_tree_indices):
            if len(neighbors) > 0:
                confidence_tree[i, neighbors] = confidence[i, neighbors]
        return confidence_tree.tocsr()


def sga_degrees(adata):
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
    g = nx.Graph(adata.uns['sga_adjacency_full_confidence'])
    degrees = [d for _, d in g.degree(weight='weight')]
    return degrees


def sga_expression_entropies(adata):
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
                                                     key=adata.uns['sga_groups_key'])
    entropies = []
    for mask in groups_masks:
        X_mask = adata.X[mask]
        x_median = np.median(X_mask, axis=0)
        x_probs = (x_median - np.min(x_median)) / (np.max(x_median) - np.min(x_median))
        entropies.append(entropy(x_probs))
    return entropies


def sga_compare_paths(adata1, adata2,
                      adjacency_key='sga_adjacency_full_confidence'):
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
    asso_groups1 = utils.identify_groups(adata1.obs['sga_groups'].values,
                                         adata2.obs['sga_groups'].values)
    asso_groups2 = utils.identify_groups(adata2.obs['sga_groups'].values,
                                         adata1.obs['sga_groups'].values)
    orig_names1 = adata1.uns['sga_groups_order_original']
    orig_names2 = adata2.uns['sga_groups_order_original']

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
    Result = namedtuple('sga_compare_paths_result',
                        ['frac_steps', 'n_steps', 'frac_paths', 'n_paths'])
    return Result(frac_steps=n_agreeing_steps/n_steps if n_steps > 0 else np.nan,
                  n_steps=n_steps if n_steps > 0 else np.nan,
                  frac_paths=n_agreeing_paths/n_paths if n_steps > 0 else np.nan,
                  n_paths=n_paths if n_steps > 0 else np.nan)


def sga_contract_graph(adata, min_group_size=0.01, max_n_contractions=1000, copy=False):
    """Contract the abstracted graph.
    """
    adata = adata.copy() if copy else adata
    if 'sga_adjacency_tree_confidence' not in adata.uns: raise ValueError('run tool sga first!')
    min_group_size = min_group_size if min_group_size >= 1 else int(min_group_size * adata.n_obs)
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

    size_before = adata.uns['sga_adjacency_tree_confidence'].shape[0]
    adata.uns['sga_adjacency_tree_confidence'], adata.obs['sga_groups'] = contract_nodes(
        adata.uns['sga_adjacency_tree_confidence'], adata.obs['sga_groups'].values)
    adata.uns['sga_groups_order'] = np.unique(adata.obs['sga_groups'].values)
    for key in ['sga_adjacency_full_confidence', 'sga_groups_original',
                'sga_groups_order_original', 'sga_groups_colors_original']:
        if key in adata.uns: del adata.uns[key]
    logg.info('    contracted graph from {} to {} nodes'
              .format(size_before, adata.uns['sga_adjacency_tree_confidence'].shape[0]))
    logg.msg('removed adata.uns["sga_adjacency_full_confidence"]', v=4)
    return adata if copy else None

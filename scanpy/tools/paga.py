from collections import namedtuple
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import minimum_spanning_tree
from .. import logging as logg
from ..neighbors import Neighbors
from .. import utils
from .. import settings


def paga(
        adata,
        groups='louvain',
        use_rna_velocity=False,
        copy=False):
    """\
    Generate cellular maps of differentiation manifolds with complex
    topologies [Wolf17i]_.

    Partition-based graph abstraction (PAGA) quantifies the connectivities of
    partitions of a neighborhood graph of single cells, thereby generating a
    much simpler abstracted graph whose nodes label the partitions. Together
    with a random walk-based distance measure, this generates a partial
    coordinatization of data useful for exploring and explaining its variation.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    groups : categorical annotation of observations or 'louvain_groups', optional (default: 'louvain_groups')
        Criterion to determine the resulting partitions of the single-cell
        graph. 'louvain_groups' uses the Louvain algorithm and optimizes
        modularity of the graph. You can also pass your predefined groups by
        choosing any categorical annotation of observations (`adata.obs`).
    use_rna_velocity : `bool` (default: `False`)
        Use RNA velocity to orient edges in the abstracted graph and estimate transitions.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy. Otherwise, perform
        computation inplace and return `None`.

    Returns
    -------
    Returns or updates `adata` depending on `copy` with
    connectivities : np.ndarray (adata.uns['connectivities'])
        The full adjacency matrix of the abstracted graph, weights
        correspond to connectivities.
    confidence : np.ndarray (adata.uns['confidence'])
        The full adjacency matrix of the abstracted graph, weights
        correspond to confidence in the presence of an edge.
    confidence_tree : sc.sparse csr matrix (adata.uns['confidence_tree'])
        The adjacency matrix of the tree-like subgraph that best explains
        the topology.
    """
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    logg.info('running partition-based graph abstraction (PAGA)', reset=True)
    paga = PAGA(adata, groups, use_rna_velocity=use_rna_velocity)
    paga.compute()
    # only add if not present
    if 'paga' not in adata.uns:
        adata.uns['paga'] = {}
    if not use_rna_velocity:
        adata.uns['paga']['connectivities'] = paga.connectivities_coarse
        adata.uns['paga']['confidence'] = paga.confidence
        adata.uns['paga']['confidence_tree'] = paga.confidence_tree
        adata.uns[groups + '_sizes'] = np.array(paga.vc.sizes())
    else:
        adata.uns['paga']['transitions_confidence'] = paga.transitions_confidence
        adata.uns['paga']['transitions_ttest'] = paga.transitions_ttest
    adata.uns['paga']['groups'] = groups
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    if use_rna_velocity:
        logg.hint(
            'added\n'
            '    \'paga/transitions_confidence\', confidence adjacency (adata.uns)\n'
            '    \'paga/transitions_ttest\', confidence subtree (adata.uns)')
    else:
        logg.hint(
            'added\n'
            '    \'paga/connectivities\', connectivities adjacency (adata.uns)\n'
            '    \'paga/confidence\', confidence adjacency (adata.uns)\n'
            '    \'paga/confidence_tree\', confidence subtree (adata.uns)')
    return adata if copy else None


class PAGA(Neighbors):

    def __init__(self, adata, groups, use_rna_velocity=False,
                 tree_based_confidence=False):
        super(PAGA, self).__init__(adata)
        self._groups = groups
        self._tree_based_confidence = tree_based_confidence
        self._use_rna_velocity = use_rna_velocity

    def compute(self):
        if self._use_rna_velocity:
            self.compute_transitions_coarse()
        else:
            self.compute_connectivities_coarse()
            self.compute_confidence()

    def compute_connectivities_coarse(self):
        import igraph
        ones = self.connectivities.copy()
        # graph where edges carry weight 1
        ones.data = np.ones(len(ones.data))
        g = utils.get_igraph_from_adjacency(ones)
        self.vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups].cat.codes.values)
        cg = self.vc.cluster_graph(combine_edges='sum')
        self.connectivities_coarse = utils.get_sparse_from_igraph(cg, weight_attr='weight')/2

    def compute_confidence(self):
        """Translates the connectivities_coarse measure into a confidence measure.
        """
        pseudo_distance = self.connectivities_coarse.copy()
        pseudo_distance.data = 1./pseudo_distance.data
        connectivities_coarse_tree = minimum_spanning_tree(pseudo_distance)
        connectivities_coarse_tree.data = 1./connectivities_coarse_tree.data
        connectivities_coarse_tree_indices = [
            connectivities_coarse_tree[i].nonzero()[1]
            for i in range(connectivities_coarse_tree.shape[0])]
        # inter- and intra-cluster based confidence
        if not self._tree_based_confidence:
            total_n = self.n_neighbors * np.array(self.vc.sizes())
            maximum = self.connectivities_coarse.max()
            confidence = self.connectivities_coarse.copy()  # initializing
            for i in range(self.connectivities_coarse.shape[0]):
                for j in range(i+1, self.connectivities_coarse.shape[1]):
                    if self.connectivities_coarse[i, j] > 0:
                        geom_mean = np.sqrt(total_n[i] * total_n[j])
                        confidence[i, j] = self.connectivities_coarse[i, j] / geom_mean
                        confidence[j, i] = confidence[i, j]
        # tree-based confidence
        else:
            median_connectivities_coarse_tree = np.median(connectivities_coarse_tree.data)
            confidence = self.connectivities_coarse.copy()
            confidence.data[self.connectivities_coarse.data >= median_connectivities_coarse_tree] = 1
            connectivities_coarse_adjusted = self.connectivities_coarse.copy()
            connectivities_coarse_adjusted.data -= median_connectivities_coarse_tree
            connectivities_coarse_adjusted.data = np.exp(connectivities_coarse_adjusted.data)
            index = self.connectivities_coarse.data < median_connectivities_coarse_tree
            confidence.data[index] = connectivities_coarse_adjusted.data[index]
        confidence_tree = self.compute_confidence_tree(
            confidence, connectivities_coarse_tree_indices)
        self.confidence = confidence
        self.confidence_tree = confidence_tree

    def compute_confidence_tree(
            self, confidence, connectivities_coarse_tree_indices):
        confidence_tree = sp.sparse.lil_matrix(confidence.shape, dtype=float)
        for i, neighbors in enumerate(connectivities_coarse_tree_indices):
            if len(neighbors) > 0:
                confidence_tree[i, neighbors] = confidence[i, neighbors]
        return confidence_tree.tocsr()

    def compute_transitions_coarse(self):
        # analogous code using networkx
        # membership = adata.obs['clusters'].cat.codes.tolist()
        # partition = defaultdict(list)
        # for n, p in zip(list(range(len(G))), membership):
        #     partition[p].append(n)
        # partition = partition.values()
        # g_abstracted = nx.quotient_graph(g, partition, relabel=True)
        # for some reason, though, edges aren't oriented in the quotient
        # graph...
        import igraph
        g = utils.get_igraph_from_adjacency(
            self._adata.uns['velocyto_transitions'], directed=True)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups].cat.codes.values)
        cg_full = vc.cluster_graph(combine_edges=False)

        g_bool = utils.get_igraph_from_adjacency(
            self._adata.uns['velocyto_transitions'].astype('bool'), directed=True)
        vc_bool = igraph.VertexClustering(
            g_bool, membership=self._adata.obs[self._groups].cat.codes.values)
        cg_bool = vc_bool.cluster_graph(combine_edges='sum')  # collapsed version
        transitions_coarse = utils.get_sparse_from_igraph(cg_bool, weight_attr='weight')
        # translate this into a confidence measure
        # the number of outgoing edges
        # total_n = np.zeros(len(vc.sizes()))
        # # (this is not the convention of standard stochastic matrices)
        # total_outgoing = transitions_coarse.sum(axis=1)
        # for i in range(len(total_n)):
        #     total_n[i] = vc.subgraph(i).ecount()
        #     total_n[i] += total_outgoing[i, 0]
        # use the topology based reference, the velocity one might have very small numbers
        total_n = self.n_neighbors * np.array(vc_bool.sizes())
        transitions_ttest = transitions_coarse.copy()
        transitions_confidence = transitions_coarse.copy()
        from scipy.stats import ttest_1samp
        for i in range(transitions_coarse.shape[0]):
            # no symmetry in transitions_coarse, hence we should not restrict to
            # upper triangle
            neighbors = transitions_coarse[i].nonzero()[1]
            for j in neighbors:
                forward = cg_full.es.select(_source=i, _target=j)['weight']
                backward = cg_full.es.select(_source=j, _target=i)['weight']
                # backward direction: add minus sign
                values = np.array(list(forward) + list(-np.array(backward)))
                # require some minimal number of observations
                if len(values) < 5:
                    transitions_ttest[i, j] = 0
                    transitions_ttest[j, i] = 0
                    transitions_confidence[i, j] = 0
                    transitions_confidence[j, i] = 0
                    continue
                t, prob = ttest_1samp(values, 0.0)
                if t > 0:
                    # number of outgoing edges greater than number of ingoing edges
                    # i.e., transition from i to j
                    transitions_ttest[i, j] = -np.log10(max(prob, 1e-10))
                    transitions_ttest[j, i] = 0
                else:
                    transitions_ttest[j, i] = -np.log10(max(prob, 1e-10))
                    transitions_ttest[i, j] = 0
                # geom_mean
                geom_mean = np.sqrt(total_n[i] * total_n[j])
                diff = (len(forward) - len(backward)) / geom_mean
                if diff > 0:
                    transitions_confidence[i, j] = diff
                    transitions_confidence[j, i] = 0
                else:
                    transitions_confidence[j, i] = -diff
                    transitions_confidence[i, j] = 0
        transitions_ttest.eliminate_zeros()
        transitions_confidence.eliminate_zeros()
        # transpose in order to match convention of stochastic matrices
        # entry ij means transition from j to i
        self.transitions_ttest = transitions_ttest.T
        self.transitions_confidence = transitions_confidence.T


def paga_degrees(adata):
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
    g = nx.Graph(adata.uns['paga']['confidence'])
    degrees = [d for _, d in g.degree(weight='weight')]
    return degrees


def paga_expression_entropies(adata):
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
    groups_order, groups_masks = utils.select_groups(
        adata, key=adata.uns['paga']['groups'])
    entropies = []
    for mask in groups_masks:
        X_mask = adata.X[mask]
        x_median = np.median(X_mask, axis=0)
        x_probs = (x_median - np.min(x_median)) / (np.max(x_median) - np.min(x_median))
        entropies.append(entropy(x_probs))
    return entropies


def paga_compare_paths(adata1, adata2,
                       adjacency_key='confidence', adjacency_key2=None):
    """Compare paths in abstracted graphs in two datasets.

    Compute the fraction of consistent paths between leafs, a measure for the
    topological similarity between graphs.

    By increasing the verbosity to level 4 and 5, the paths that do not agree
    and the paths that agree are written to the output, respectively.

    The PAGA "groups key" needs to be the same in both objects.

    Parameters
    ----------
    adata1, adata2 : AnnData
        Annotated data matrices to compare.
    adjacency_key : str
        Key for indexing the adjacency matrices in `.uns['paga']` to be used in
        adata1 and adata2.
    adjacency_key2 : str, None
        If provided, used for adata2.


    Returns
    -------
    OrderedTuple with attributes ``n_steps`` (total number of steps in paths)
    and ``frac_steps`` (fraction of consistent steps), ``n_paths`` and
    ``frac_paths``.
    """
    import networkx as nx
    g1 = nx.Graph(adata1.uns['paga'][adjacency_key])
    g2 = nx.Graph(adata2.uns['paga'][adjacency_key2 if adjacency_key2 is not None else adjacency_key])
    leaf_nodes1 = [str(x) for x in g1.nodes() if g1.degree(x) == 1]
    logg.msg('leaf nodes in graph 1: {}'.format(leaf_nodes1), v=5, no_indent=True)
    paga_groups = adata1.uns['paga']['groups']
    asso_groups1 = utils.identify_groups(adata1.obs[paga_groups].values,
                                         adata2.obs[paga_groups].values)
    asso_groups2 = utils.identify_groups(adata2.obs[paga_groups].values,
                                         adata1.obs[paga_groups].values)
    orig_names1 = adata1.obs[paga_groups].cat.categories
    orig_names2 = adata2.obs[paga_groups].cat.categories

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
            logg.msg('there are no connecting paths in both graphs', v=5, no_indent=True)
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
    Result = namedtuple('paga_compare_paths_result',
                        ['frac_steps', 'n_steps', 'frac_paths', 'n_paths'])
    return Result(frac_steps=n_agreeing_steps/n_steps if n_steps > 0 else np.nan,
                  n_steps=n_steps if n_steps > 0 else np.nan,
                  frac_paths=n_agreeing_paths/n_paths if n_steps > 0 else np.nan,
                  n_paths=n_paths if n_steps > 0 else np.nan)

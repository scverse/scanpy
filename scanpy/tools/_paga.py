from collections import namedtuple
from typing import List

import numpy as np
import scipy as sp
from scipy.sparse.csgraph import minimum_spanning_tree

from .. import utils
from .. import logging as logg
from ..logging import _settings_verbosity_greater_or_equal_than
from ..neighbors import Neighbors

_AVAIL_MODELS = {'v1.0', 'v1.2'}


def paga(
        adata,
        groups='louvain',
        use_rna_velocity=False,
        model='v1.2',
        copy=False):
    """Mapping out the coarse-grained connectivity structures of complex manifolds [Wolf19]_.

    By quantifying the connectivity of partitions (groups, clusters) of the
    single-cell graph, partition-based graph abstraction (PAGA) generates a much
    simpler abstracted graph (*PAGA graph*) of partitions, in which edge weights
    represent confidence in the presence of connections. By tresholding this
    confidence in :func:`~scanpy.pl.paga`, a much simpler representation of the
    manifold data is obtained, which is nonetheless faithful to the topology of
    the manifold.

    The confidence should be interpreted as the ratio of the actual versus the
    expected value of connetions under the null model of randomly connecting
    partitions. We do not provide a p-value as this null model does not
    precisely capture what one would consider "connected" in real data, hence it
    strongly overestimates the expected value. See an extensive discussion of
    this in [Wolf19]_.

    .. note::
        Note that you can use the result of :func:`~scanpy.pl.paga` in
        :func:`~scanpy.tl.umap` and :func:`~scanpy.tl.draw_graph` via
        `init_pos='paga'` to get single-cell embeddings that are typically more
        faithful to the global topology.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        An annotated data matrix.
    groups : key for categorical in `adata.obs`, optional (default: 'louvain')
        You can pass your predefined groups by choosing any categorical
        annotation of observations (`adata.obs`).
    use_rna_velocity : `bool` (default: `False`)
        Use RNA velocity to orient edges in the abstracted graph and estimate
        transitions. Requires that `adata.uns` contains a directed single-cell
        graph with key `['velocity_graph']`. This feature might be subject
        to change in the future.
    model : {'v1.2', 'v1.0'}, optional (default: 'v1.2')
        The PAGA connectivity model.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy. Otherwise, perform
        computation inplace and return `None`.

    Returns
    -------
    **connectivities** : :class:`numpy.ndarray` (adata.uns['connectivities'])
        The full adjacency matrix of the abstracted graph, weights correspond to
        confidence in the connectivities of partitions.
    **connectivities_tree** : :class:`scipy.sparse.csr_matrix` (adata.uns['connectivities_tree'])
        The adjacency matrix of the tree-like subgraph that best explains
        the topology.

    Notes
    -----
    Together with a random walk-based distance measure
    (e.g. :func:`scanpy.tl.dpt`) this generates a partial coordinatization of
    data useful for exploring and explaining its variation.

    See Also
    --------
    pl.paga
    pl.paga_path
    pl.paga_compare
    """
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    adata = adata.copy() if copy else adata
    utils.sanitize_anndata(adata)
    logg.info('running PAGA', reset=True)
    paga = PAGA(adata, groups, model=model)
    # only add if not present
    if 'paga' not in adata.uns:
        adata.uns['paga'] = {}
    if not use_rna_velocity:
        paga.compute_connectivities()
        adata.uns['paga']['connectivities'] = paga.connectivities
        adata.uns['paga']['connectivities_tree'] = paga.connectivities_tree
        # adata.uns['paga']['expected_n_edges_random'] = paga.expected_n_edges_random
        adata.uns[groups + '_sizes'] = np.array(paga.ns)
    else:
        paga.compute_transitions()
        adata.uns['paga']['transitions_confidence'] = paga.transitions_confidence
        # adata.uns['paga']['transitions_ttest'] = paga.transitions_ttest
    adata.uns['paga']['groups'] = groups
    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    if use_rna_velocity:
        logg.hint(
            'added\n'
            '    \'paga/transitions_confidence\', connectivities adjacency (adata.uns)')
        #    '    \'paga/transitions_ttest\', t-test on transitions (adata.uns)')
    else:
        logg.hint(
            'added\n'
            '    \'paga/connectivities\', connectivities adjacency (adata.uns)\n'
            '    \'paga/connectivities_tree\', connectivities subtree (adata.uns)')
    return adata if copy else None


class PAGA():

    def __init__(self, adata, groups, model='v1.2'):
        self._adata = adata
        self._neighbors = Neighbors(adata)
        self._groups_key = groups
        self._model = model

    def compute_connectivities(self):
        if self._model == 'v1.2':
            return self._compute_connectivities_v1_2()
        elif self._model == 'v1.0':
            return self._compute_connectivities_v1_0()
        else:
            raise ValueError(
                '`model` {} needs to be one of {}.'
                .format(self._model, _AVAIL_MODELS))

    def _compute_connectivities_v1_2(self):
        import igraph
        ones = self._neighbors.distances.copy()
        ones.data = np.ones(len(ones.data))
        # should be directed if we deal with distances
        g = utils.get_igraph_from_adjacency(ones, directed=True)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups_key].cat.codes.values)
        ns = vc.sizes()
        n = sum(ns)
        es_inner_cluster = [vc.subgraph(i).ecount() for i in range(len(ns))]
        cg = vc.cluster_graph(combine_edges='sum')
        inter_es = utils.get_sparse_from_igraph(cg, weight_attr='weight')
        es = np.array(es_inner_cluster) + inter_es.sum(axis=1).A1
        inter_es = inter_es + inter_es.T  # \epsilon_i + \epsilon_j
        connectivities = inter_es.copy()
        expected_n_edges = inter_es.copy()
        inter_es = inter_es.tocoo()
        for i, j, v in zip(inter_es.row, inter_es.col, inter_es.data):
            expected_random_null = (es[i]*ns[j] + es[j]*ns[i])/(n - 1)
            if expected_random_null != 0:
                scaled_value = v / expected_random_null
            else:
                scaled_value = 1
            if scaled_value > 1:
                scaled_value = 1
            connectivities[i, j] = scaled_value
            expected_n_edges[i, j] = expected_random_null
        # set attributes
        self.ns = ns
        self.expected_n_edges_random = expected_n_edges
        self.connectivities = connectivities
        self.connectivities_tree = self._get_connectivities_tree_v1_2()
        return inter_es.tocsr(), connectivities

    def _compute_connectivities_v1_0(self):
        import igraph
        ones = self._neighbors.connectivities.copy()
        ones.data = np.ones(len(ones.data))
        g = utils.get_igraph_from_adjacency(ones)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups_key].cat.codes.values)
        ns = vc.sizes()
        cg = vc.cluster_graph(combine_edges='sum')
        inter_es = utils.get_sparse_from_igraph(cg, weight_attr='weight')/2
        connectivities = inter_es.copy()
        inter_es = inter_es.tocoo()
        n_neighbors_sq = self._neighbors.n_neighbors**2
        for i, j, v in zip(inter_es.row, inter_es.col, inter_es.data):
            # have n_neighbors**2 inside sqrt for backwards compat
            geom_mean_approx_knn = np.sqrt(
                n_neighbors_sq * ns[i] * ns[j])
            if geom_mean_approx_knn != 0:
                scaled_value = v / geom_mean_approx_knn
            else:
                scaled_value = 1
            connectivities[i, j] = scaled_value
        # set attributes
        self.ns = ns
        self.connectivities = connectivities
        self.connectivities_tree = self._get_connectivities_tree_v1_0(inter_es)
        return inter_es.tocsr(), connectivities

    def _get_connectivities_tree_v1_2(self):
        inverse_connectivities = self.connectivities.copy()
        inverse_connectivities.data = 1./inverse_connectivities.data
        connectivities_tree = minimum_spanning_tree(inverse_connectivities)
        connectivities_tree_indices = [
            connectivities_tree[i].nonzero()[1]
            for i in range(connectivities_tree.shape[0])]
        connectivities_tree = sp.sparse.lil_matrix(self.connectivities.shape, dtype=float)
        for i, neighbors in enumerate(connectivities_tree_indices):
            if len(neighbors) > 0:
                connectivities_tree[i, neighbors] = self.connectivities[i, neighbors]
        return connectivities_tree.tocsr()

    def _get_connectivities_tree_v1_0(self, inter_es):
        inverse_inter_es = inter_es.copy()
        inverse_inter_es.data = 1./inverse_inter_es.data
        connectivities_tree = minimum_spanning_tree(inverse_inter_es)
        connectivities_tree_indices = [
            connectivities_tree[i].nonzero()[1]
            for i in range(connectivities_tree.shape[0])]
        connectivities_tree = sp.sparse.lil_matrix(inter_es.shape, dtype=float)
        for i, neighbors in enumerate(connectivities_tree_indices):
            if len(neighbors) > 0:
                connectivities_tree[i, neighbors] = self.connectivities[i, neighbors]
        return connectivities_tree.tocsr()

    def compute_transitions(self):
        vkey = 'velocity_graph'
        if vkey not in self._adata.uns:
            if 'velocyto_transitions' in self._adata.uns:
                self._adata.uns[vkey] = self._adata.uns['velocyto_transitions']
                logg.msg('The key \'velocyto_transitions\' has been changed to \'velocity_graph\'.', v=5, no_indent=True)
            else:
                raise ValueError(
                    'The passed AnnData needs to have an `uns` annotation '
                    'with key \'velocity_graph\' - a sparse matrix from RNA velocity.')
        if self._adata.uns[vkey].shape != (self._adata.n_obs, self._adata.n_obs):
            raise ValueError(
                'The passed \'velocity_graph\' have shape {} but shoud have shape {}'
                .format(self._adata.uns[vkey].shape,
                        (self._adata.n_obs, self._adata.n_obs)))
        # restore this at some point
        # if 'expected_n_edges_random' not in self._adata.uns['paga']:
        #     raise ValueError(
        #         'Before running PAGA with `use_rna_velocity=True`, run it with `False`.')
        import igraph
        g = utils.get_igraph_from_adjacency(
            self._adata.uns[vkey].astype('bool'), directed=True)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups_key].cat.codes.values)
        # set combine_edges to False if you want self loops
        cg_full = vc.cluster_graph(combine_edges='sum')
        transitions = utils.get_sparse_from_igraph(cg_full, weight_attr='weight')
        transitions = transitions - transitions.T
        transitions_conf = transitions.copy()
        transitions = transitions.tocoo()
        total_n = self._neighbors.n_neighbors * np.array(vc.sizes())
        # total_n_sum = sum(total_n)
        # expected_n_edges_random = self._adata.uns['paga']['expected_n_edges_random']
        for i, j, v in zip(transitions.row, transitions.col, transitions.data):
            # if expected_n_edges_random[i, j] != 0:
            #     # factor 0.5 because of asymmetry
            #     reference = 0.5 * expected_n_edges_random[i, j]
            # else:
            #     # approximate
            #     reference = self._neighbors.n_neighbors * total_n[i] * total_n[j] / total_n_sum
            reference = np.sqrt(total_n[i] * total_n[j])
            transitions_conf[i, j] = 0 if v < 0 else v / reference
        transitions_conf.eliminate_zeros()
        # transpose in order to match convention of stochastic matrices
        # entry ij means transition from j to i
        self.transitions_confidence = transitions_conf.T

    def compute_transitions_old(self):
        import igraph
        g = utils.get_igraph_from_adjacency(
            self._adata.uns['velocyto_transitions'], directed=True)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups_key].cat.codes.values)
        # this stores all single-cell edges in the cluster graph
        cg_full = vc.cluster_graph(combine_edges=False)
        # this is the boolean version that simply counts edges in the clustered graph
        g_bool = utils.get_igraph_from_adjacency(
            self._adata.uns['velocyto_transitions'].astype('bool'), directed=True)
        vc_bool = igraph.VertexClustering(
            g_bool, membership=self._adata.obs[self._groups_key].cat.codes.values)
        cg_bool = vc_bool.cluster_graph(combine_edges='sum')  # collapsed version
        transitions = utils.get_sparse_from_igraph(cg_bool, weight_attr='weight')
        total_n = self._neighbors.n_neighbors * np.array(vc_bool.sizes())
        transitions_ttest = transitions.copy()
        transitions_confidence = transitions.copy()
        from scipy.stats import ttest_1samp
        for i in range(transitions.shape[0]):
            neighbors = transitions[i].nonzero()[1]
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


def paga_degrees(adata) -> List[int]:
    """Compute the degree of each node in the abstracted graph.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    List of degrees for each node.
    """
    import networkx as nx
    g = nx.Graph(adata.uns['paga']['connectivities'])
    degrees = [d for _, d in g.degree(weight='weight')]
    return degrees


def paga_expression_entropies(adata) -> List[float]:
    """Compute the median expression entropy for each node-group.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    Entropies of median expressions for each node.
    """
    from scipy.stats import entropy
    groups_order, groups_masks = utils.select_groups(
        adata, key=adata.uns['paga']['groups'])
    entropies = []
    for mask in groups_masks:
        X_mask = adata.X[mask].todense()
        x_median = np.nanmedian(X_mask, axis=1,overwrite_input=True)
        x_probs = (x_median - np.nanmin(x_median)) / (np.nanmax(x_median) - np.nanmin(x_median))
        entropies.append(entropy(x_probs))
    return entropies


def paga_compare_paths(adata1, adata2,
                       adjacency_key='connectivities', adjacency_key2=None):
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

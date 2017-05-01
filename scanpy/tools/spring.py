# coding: utf-8
# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Force-directed drawing of Data Graph

By default, the Fruchterman-Reingold algorithm is used. This is simple
force-directed graph drawing.

References
----------
- General: https://en.wikipedia.org/wiki/Force-directed_graph_drawing
- Suggested for drawing knn-graphs in the context of single-cell
  transcriptomics: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
"""

import numpy as np
from .. import settings as sett
from .. import utils

step_size = 10


def spring(adata, k=4, n_comps=2, n_steps=12, rep=None, copy=False):
    u"""
    Visualize data using the force-directed Fruchterman-Reingold algorithm.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata.smp['X_pca']: np.ndarray
            Result of preprocessing with PCA: observations × variables.
            If it exists, spring will use this instead of adata.X.
    k : int
        Number of nearest neighbors in graph.
    n_comps : int
        Number of principal components in preprocessing PCA.
    n_steps : int
        Number of steps in optimizing the spring representation.
    rep : float of None (default: None)
        Repulsion = strength of springs. If None the distance is set to
        1/sqrt(n) where n is the number of cells. Increase this value to move
        nodes farther apart.

    Returns
    -------
    X_spring : np.ndarray
         Force-directed graph drawing representation of the data with shape
         n_variables x n_comps.
    """
    adata = adata.copy() if copy else adata
    sett.m(0, 'draw knn graph')
    if 'X_pca' in adata.smp:
        X = adata.smp['X_pca']
        sett.m(0, '--> using X_pca for building graph')
    else:
        X = adata.X
        sett.m(0, '--> using X for building graph')
    D = utils.comp_distance(X, metric='euclidean')
    # deterimine the distance of the k nearest neighbors
    indices = np.zeros((D.shape[0], k), dtype=np.int_)
    for irow, row in enumerate(D):
        # the last item is already in its sorted position as
        # argpartition puts the (k-1)th element - starting to count from
        # zero - in its sorted position
        idcs = np.argpartition(row, k-1)[:k]
        indices[irow] = idcs
    # compute adjacency matrix
    # make this float, as we might put a weight matrix here
    Adj = np.zeros(D.shape, dtype=float)
    for irow, row in enumerate(indices):
        Adj[irow, row] = 1
        # symmetrize as in DPT
        # for j in row:
        #    if irow not in indices[j]:
        #        Adj[j,irow] = 1
    # just sample initial positions, the rest is done by the plotting tool
    np.random.seed(1)
    Y = np.asarray(np.random.random((Adj.shape[0], 2)), dtype=Adj.dtype)
    sett.m(0, 'is very slow, will be sped up soon')
    for istep in 1 + np.arange(n_steps, dtype=int):
        sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
        Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size, rep=rep)
    adata.smp['X_spring'] = Y
    return adata if copy else None


def fruchterman_reingold_layout_networkX(Adj, Yinit=None):
    """
    Wrapper around networkX function.

    Is equivalent to what is hard-coded in this module up to a shift of the
    coordinates.
    """
    try:
        import networkx as nx
        print(nx.__file__)
    except ImportError:
        msg = ('importing package networkx for graph drawing failed\n'
               '--> install via "pip install networkx"')
        raise ImportError(msg)
    G = nx.DiGraph(Adj)
    posinit = dict(zip(G, Yinit))
    pos = nx.fruchterman_reingold_layout(G, pos=posinit)
    # map to numpy array as in other visualization tools
    Y = np.zeros((Adj.shape[0], 2))
    for k, y in pos.items():
        Y[k] = y
    return Y

# ------------------------------------------------------------------------------
# The following is largely adapted from NetworkX graph drawing, see copyright
# ------------------------------------------------------------------------------
#
#    Copyright (C) 2004-2016 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg <aric.hagberg@gmail.com>,
#          Dan Schult <dschult@colgate.edu>

def fruchterman_reingold_layout(W,
                                rep=None,
                                Yinit=None,
                                fixed=None,
                                iterations=50,
                                scale=1.0,
                                dim=2):
    """
    Position nodes using Fruchterman-Reingold force-directed algorithm.

    Parameters
    ----------
    W : np.ndarray
        Weight matrix (weighted adjacency matrix).
    rep : float (default=None)
        Repulsion = optimal distance between nodes = the strength of springs.
        If None the distance is set to 1/sqrt(n) where n is the number of nodes.
        Increase this value to move nodes farther apart.
    Yinit : np.ndarray or None  optional (default=None)
        Initial positions as an array with shape (W.shape[1] x 2). If None, then
        use random initial positions.
    iterations : int  optional (default=50)
        Number of iterations of spring-force relaxation
    scale : float (default=1.0)
        Scale factor for positions. The nodes are positioned
        in a box of size [0, scale] x [0, scale].
    dim : int
        Dimension of layout

    Returns
    -------
    Y : np.ndarray
        An array with positions.
    """
    if rep is None:
        nnodes, _ = W.shape
        rep = 1.0 / np.sqrt(nnodes)
    sett.m(0, 'using repulsion rep =', rep)
    if len(W) >= 500:
        Y = _fruchterman_reingold_sparse(W, rep, Yinit, fixed, iterations, dim)
    else:
        Y = _fruchterman_reingold_dense(W, rep, Yinit, fixed, iterations, dim)
    if fixed is None:
        Y = rescale_layout(Y, scale=scale)
    return Y

def _fruchterman_reingold_sparse(A, rep, Y=None, fixed=None,
                                 iterations=50, dim=2):
    try:
        from scipy.sparse import spdiags, coo_matrix
    except ImportError:
        msg = "_sparse_fruchterman_reingold() scipy numpy: http://scipy.org/ "
        raise ImportError(msg)
    # make sure we have a List of Lists representation (LInked List format)
    A = (coo_matrix(A)).tolil()
    nnodes, _ = A.shape
    if Y is None:
        # random initial positions
        Y = np.asarray(np.random.random((nnodes, dim)), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        Y = Y.astype(A.dtype)
    # the initial "temperature"  is about .1 of domain area
    # this is the largest step allowed in the dynamics.
    t = 0.1 * max(max(Y.T[0]) - min(Y.T[0]), max(Y.T[1]) - min(Y.T[1]))
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations+1)

    displacement = np.zeros((dim, nnodes))
    for iteration in range(iterations):
        displacement *= 0
        # loop over rows
        for i in range(A.shape[0]):
            # difference between this row's node position and all others
            delta = (Y[i] - Y).T
            # distance between points
            distance = np.sqrt((delta**2).sum(axis=0))
            # enforce minimum distance of 0.01
            distance = np.where(distance < 0.01, 0.01, distance)
            # the adjacency matrix row
            Ai = np.asarray(A.getrowview(i).toarray())
            # displacement "force"
            displacement[:, i] +=\
                (delta * (rep * rep / distance**2 - Ai * distance / rep)).sum(axis=1)
        # update positions
        length = np.sqrt((displacement**2).sum(axis=0))
        length = np.where(length < 0.01, 0.1, length)
        Y += (displacement * t / length).T
        # cool temperature
        t -= dt
    return Y

def _fruchterman_reingold_dense(A, rep, Y=None, fixed=None,
                                iterations=50, dim=2):
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # make sure we have an array instead of a matrix
    A = np.asarray(A)

    if Y is None:
        # random initial positions
        Y = np.asarray(np.random.random((nnodes, dim)), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        Y = Y.astype(A.dtype)

    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    # We need to calculate this in case our fixed positions force our domain
    # to be much bigger than 1x1
    t = max(max(Y.T[0]) - min(Y.T[0]), max(Y.T[1]) - min(Y.T[1]))*0.1
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t/float(iterations+1)
    delta = np.zeros((Y.shape[0], Y.shape[0], Y.shape[1]), dtype=A.dtype)
    # the inscrutable (but fast) version
    # this is still O(V^2)
    # could use multilevel methods to speed this up significantly
    for iteration in range(iterations):
        # matrix of difference between points
        for i in range(Y.shape[1]):
            delta[:, :, i] = Y[:, i, None] - Y[:, i]
        # distance between points
        distance = np.sqrt((delta**2).sum(axis=-1))
        # enforce minimum distance of 0.01
        distance = np.where(distance < 0.01, 0.01, distance)
        # displacement "force"
        displacement = np.transpose(np.transpose(delta) *
                                    (rep * rep / distance**2 - A * distance / rep)
                                    ).sum(axis=1)
        # update positions
        length = np.sqrt((displacement**2).sum(axis=1))
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = np.transpose(np.transpose(displacement) * t / length)
        if fixed is not None:
            # don't change positions of fixed nodes
            delta_pos[fixed] = 0.0
        Y += delta_pos
        # cool temperature
        t -= dt
    return Y

def rescale_layout(Y, scale=1):
    """Return scaled position array to (-scale, scale) in all axes.

    The function acts on NumPy arrays which hold position information.
    Each position is one row of the array. The dimension of the space
    equals the number of columns. Each coordinate in one column.

    To rescale, the mean (center) is subtracted from each axis separately.
    Then all values are scaled so that the largest magnitude value
    from all axes equals `scale` (thus, the aspect ratio is preserved).
    The resulting NumPy Array is returned (order of rows unchanged).

    Parameters
    ----------
    Y : numpy array
        positions to be scaled. Each row is a position.
    scale : number (default: 1)
        The size of the resulting extent in all directions.

    Returns
    -------
    Y : numpy array
        Scaled positions. Each row is a position.
    """
    # Find max length over all dimensions
    lim = 0  # max coordinate for all axes
    for i in range(Y.shape[1]):
        Y[:, i] -= Y[:, i].mean()
        lim = max(Y[:, i].max(), lim)
    # rescale to (-scale, scale) in all directions, preserves aspect
    if lim > 0:
        for i in range(Y.shape[1]):
            Y[:, i] *= scale / lim

    return Y

# coding: utf-8
# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Draw the Data Graph

By default, the Fruchterman-Reingold algorithm is used. This is simple
force-directed graph drawing.

References
----------
- General: https://en.wikipedia.org/wiki/Force-directed_graph_drawing
- Suggested for drawing knn-graphs in the context of single-cell
  transcriptomics: Weinreb et al., bioRxiv doi:10.1101/090332 (2016) 
"""

from __future__ import absolute_import
from collections import OrderedDict as odict
import numpy as np
from ..compat.matplotlib import pyplot as pl
from ..tools.pca import pca
from .. import settings as sett
from .. import plotting as plott
from .. import utils

step_size = 10

def drawg(adata, k=4, n_comps=2):
    u"""
    Visualize data using graph drawing algorithms.

    In particular the force-directed Fruchterman-Reingold algorithm.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata['X_pca']: np.ndarray
            Result of preprocessing with PCA: observations × variables.
            If it exists, drawg will use this instead of adata.X.
    k : int
        Number of nearest neighbors in graph.
    n_comps : int
        Number of principal components in preprocessing PCA.

    Returns
    -------
    d : dict containing
        Y : np.ndarray
            Fruchterman Reingold representation of data.
    """
    sett.m(0,'draw knn graph')
    if 'X_pca' in adata:
        X = adata['X_pca']
        sett.m(0, '--> using X_pca for building graph')
    else:
        X = adata.X
        sett.m(0, '--> using X for building graph')
    D = utils.comp_distance(X, metric='euclidean')
    # deterimine the distance of the k nearest neighbors
    indices = np.zeros((D.shape[0],k),dtype=np.int_)
    for irow, row in enumerate(D):
        # the last item is already in its sorted position as
        # argpartition puts the (k-1)th element - starting to count from
        # zero - in its sorted position
        idcs = np.argpartition(row,k-1)[:k]
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
    n_steps = 12
    for istep in 1 + np.arange(n_steps, dtype=int):
        sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
        Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size)
    return {'type': 'drawg', 'Y': Y, 'Adj': Adj, 'istep': istep}

def plot(ddrawg, adata,
         add_steps=0,
         smp=None,
         names=None,
         comps='1,2',
         cont=None,
         layout='2d',
         legendloc='lower right',
         cmap=None,
         right_margin=0.75):
    """
    Scatter plots.

    Parameters
    ----------
    ddrawg : dict
        Dict returned by diffmap tool.
    adata : AnnData
        Annotated data matrix.
    add_steps : int
        Steps to iterate graph drawing algorithm.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: continuous: viridis/ categorical: finite palette)
         String denoting matplotlib color map.
    right_margin : float (default: 0.2)
         Adjust how far the plotting panel extends to the right.
    """
    params = locals(); del params['adata']; del params['ddrawg']
    Y = ddrawg['Y']

    if params['add_steps'] == 0:
        del params['add_steps']
        sett.m(0, 'set parameter add_steps > 0 to iterate. ' 
               'the current step is', ddrawg['istep'],
               '\n--> append, for example, "--plotparams add_steps 1", for a single step')
        istep = ddrawg['istep']
        _plot(ddrawg, adata, istep, **params)
    else:
        Adj = ddrawg['Adj']
        istep = ddrawg['istep']
        # TODO: don't save the adjacency matrix!!!
        import scanpy as sc
        sc.write(ddrawg['writekey']+'_step{:02}'.format(istep), ddrawg)
        # compute the next steps
        istep_init = istep + 1
        add_steps = params['add_steps']
        del params['add_steps']
        for istep in istep_init + np.arange(add_steps, dtype=int):
            sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
            Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size)
            sett.mt(0, 'finished computation')
            _plot({'Y': Y}, adata, istep, **params)
        # save state of Y to outfile
        ddrawg['Y'] = Y 
        ddrawg['istep'] = istep
        sc.write(ddrawg['writekey'], ddrawg)

def _plot(dplot, adata, istep=0, **params):
    from .. import plotting as plott
    plott.plot_tool(dplot, adata,
                    # defined in plotting
                    subtitles=['Fruchterman-Reingold step: ' + str(istep)],
                    component_name='FR',
                    **params)

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

def fruchterman_reingold_layout(W, rep=None,
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
    if len(W) >= 500:
        Y = _fruchterman_reingold_sparse(W, rep, Yinit, fixed, iterations, dim)
    else:
        Y = _fruchterman_reingold_dense(W, rep, Yinit, fixed, iterations, dim)
    if fixed is None:
        Y = rescale_layout(Y, scale=scale)
    return Y

def _fruchterman_reingold_sparse(A, rep, Y=None, fixed=None,
                                 iterations=50, dim=2):
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Sparse version
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

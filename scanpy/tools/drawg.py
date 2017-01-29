# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Draw the Data Graph
===================

Simple force-directed graph drawing. In particular, the Fruchterman-Reingold
algorithm. See https://en.wikipedia.org/wiki/Force-directed_graph_drawing for
some history.

Suggested for drawing knn graphs in the context of single-cell transcriptomics 
by Weinreb et al., bioRxiv doi:10.1101/090332 (2016.
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

def drawg(ddata, k=4, n_components=2):
    """
    Visualize data using graph drawing algorithms.

    In particular the force-directed Fruchterman-Reingold algorithm.

    Parameters
    ----------
    ddata : dict containing
        X : np.ndarray
            Data array, rows store observations, columns covariates.
    k : int
        Number of nearest neighbors in graph.
    n_components : int
        Number of principal components in preprocessing PCA.

    Returns
    -------
    d : dict containing
        Y : np.ndarray
            Fruchterman Reingold representation of data.
    """
    sett.m(0,'draw knn graph')
    X = ddata['X']
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
        Adj[irow,row] = 1
        # symmetrize as in DPT
        # for j in row:
        #    if irow not in indices[j]:
        #        Adj[j,irow] = 1
    
    # just sample initial positions, the rest is done by the plotting tool
    np.random.seed(1)
    Y = np.asarray(np.random.random((Adj.shape[0], 2)), dtype=Adj.dtype)
    # compute first step 
    istep = 1
    sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
    Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size)
    sett.mt(0, 'finished')
    return {'type': 'drawg', 'Y': Y, 'Adj': Adj, 'istep': istep}

def plot(ddrawg, ddata,
         nrsteps=0,
         layout='2d',
         legendloc='lower right',
         cmap='jet',
         adjust_right=0.75): # consider changing to 'viridis'
    """
    Plot the results of a DPT analysis.

    Parameters
    ----------
    ddrawg : dict
        Dict returned by diffmap tool.
    ddata : dict
        Data dictionary.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: jet)
         String denoting matplotlib color map. 
    """
    params = locals(); del params['ddata']; del params['ddrawg']
    Y = ddrawg['Y']

    if params['nrsteps'] == 0:
        sett.m(0, 'set parameter nrsteps > 0 to iterate. ' 
               'the current step is', ddrawg['istep'],
               '\n--> append, for example, "--plotparams nrsteps 1", for a single step')
        istep = ddrawg['istep']
        _plot(Y, ddata, params, istep)
        if sett.autoshow:
            pl.show()
        pl.savefig(sett.figdir+ddrawg['writekey']
                   +'_step{:02}'.format(istep)+'.'+sett.extf)
    else:
        Adj = ddrawg['Adj']
        istep = ddrawg['istep']
        # save a copy to retrieve later
        # TODO: don't save the adjacency matrix!!!
        import scanpy as sc
        sc.write(ddrawg['writekey']+'_step{:02}'.format(istep), ddrawg)
        # compute the next steps
        istep_init = istep + 1
        for istep in istep_init + np.arange(params['nrsteps'], dtype=int):
            sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
            Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size)
            sett.mt(0, 'finished computation')
            _plot(Y, ddata, params, istep)
            if sett.autoshow:
                pl.show(block=False)
                sett.mt(0, 'finished plotting')
            pl.savefig(sett.figdir+ddrawg['writekey']
                       +'_step{:02}'.format(istep)+'.'+sett.extf)
        # save state of Y to outfile
        ddrawg['Y'] = Y 
        ddrawg['istep'] = istep
        sc.write(ddrawg['writekey'], ddrawg)
        #
        if sett.autoshow:
            pl.show()

def _plot(Y, ddata, params, istep):
    # highlights
    highlights = []
    if False:
        if 'highlights' in ddata:
            highlights = ddata['highlights']
    # base figure
    axs = plott.scatter(Y,
                        subtitles=['draw graph: step ' + str(istep)],
                        component_name='FR',
                        layout=params['layout'],
                        c='grey',
                        highlights=highlights,
                        cmap=params['cmap'])
    # annotated groups
    if 'groupmasks' in ddata:
        for igroup, group in enumerate(ddata['groupmasks']):
            plott.group(axs[0], igroup, ddata, Y, params['layout'])
        axs[0].legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
        # right margin
        pl.subplots_adjust(right=params['adjust_right'])

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

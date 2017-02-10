# coding: utf-8
# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
t-SNE

References
----------
This module automatically choose from three t-SNE versions from
- sklearn.manifold.TSNE
- Dmitry Ulyanov (multicore, fastest)
  https://github.com/DmitryUlyanov/Multicore-TSNE
  install via 'pip install psutil cffi', get code from github
- Laurens van der Maaten (slowest, oldest), slow fall back option
  https://lvdmaaten.github.io/tsne/
  Copyright 2008 Laurens van der Maaten, Tilburg University.
"""

from collections import OrderedDict as odict
import numpy as np
from ..compat.matplotlib import pyplot as pl
from ..tools.pca import pca
from .. import settings as sett
from .. import plotting as plott
from .. import utils

def tsne(adata, nr_pcs=50, perplexity=30):
    u"""
    Visualize data using t-SNE as of van der Maaten & Hinton (2008).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata['Xpca']: np.ndarray
            Result of preprocessing with PCA: observations × variables.
            If there it exists, tsne will use this instead of adata.X.
    nr_pcs : int
        Number of principal components in preprocessing PCA.

    Parameters as used in sklearn.manifold.TSNE:
    perplexity : float, optional (default: 30) 
        Perplexity.

    Returns
    -------
    dtsne : dict containing
        Y : np.ndarray
            tSNE representation of the data.
    """
    params = locals(); del params['adata']
    sett.m(0,'perform tSNE')
    sett.m(0,'--> mind that this is not deterministic!')
    # preprocessing by PCA
    if 'Xpca' in adata:
        X = adata['Xpca']
        sett.m(0, 'using Xpca for tSNE')
    else:
        if params['nr_pcs'] > 0 and adata.X.shape[1] > params['nr_pcs']:
            sett.m(0, 'preprocess using PCA with', params['nr_pcs'], 'PCs')
            sett.m(0, '--> avoid this by setting nr_pcs = 0')
            dpca = pca(adata, nr_comps=params['nr_pcs'])
            X = dpca['Y']
        else:
            X = adata.X
    # params for sklearn
    params_sklearn = {k: v for k, v in params.items() if not k=='nr_pcs'}
    params_sklearn['verbose'] = sett.verbosity
    # deal with different tSNE implementations
    try:
        from MulticoreTSNE import MulticoreTSNE as TSNE
        tsne = TSNE(n_jobs=4, **params_sklearn)
        sett.m(0,'--> perform tSNE using MulticoreTSNE')
        Y = tsne.fit_transform(X)
    except ImportError:
        try:
            from sklearn.manifold import TSNE
            tsne = TSNE(**params_sklearn)
            sett.m(1,'--> perform tSNE using sklearn!')
            sett.m(1,'--> can be sped up by installing\n' 
                     '    https://github.com/DmitryUlyanov/Multicore-TSNE')
            Y = tsne.fit_transform(X)            
        except ImportError:
            sett.m(0,'--> perform tSNE using slow/unreliable original\n' 
                     '    code by L. van der Maaten!?\n'
                     '--> consider installing sklearn\n'
                     '    using "conda/pip install scikit-learn"')
            Y = _tsne_vandermaaten(X, 2, params['perplexity'])
    return {'type': 'tsne', 'Y': Y}

def plot(dplot, adata,
         smp='',
         comps='1,2',
         layout='2d',
         legendloc='lower right',
         cmap='jet',
         adjust_right=0.75):
    """
    Plot the results of a DPT analysis.

    Parameters
    ----------
    dplot : dict
        Dict returned by plotting tool.
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first anntotated group)
        Sample annotation for coloring, possible are all keys in adata.smp_keys(),
        or gene names.
    comps : str, optional (default: "1,2")
         String in the form "comp1,comp2,comp3".
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: "jet")
         String denoting matplotlib color map.
    adjust_right : float, optional (default: 0.75)
         Increase to increase the right margin.
    """
    from .. import plotting as plott
    plott.plot_tool(dplot, adata,
                    smp,
                    comps,
                    layout,
                    legendloc,
                    cmap,
                    adjust_right,
                    # defined in plotting
                    subtitles=['tSNE'],
                    component_name='tSNE')

def _tsne_vandermaaten(X = np.array([]), no_dims = 2, perplexity = 30.0):
    """
    Runs t-SNE on the dataset in the NxD array X to reduce its dimensionality to
    no_dims dimensions.  The syntaxis of the function is Y = tsne.tsne(X,
    no_dims, perplexity), where X is an NxD NumPy array.
    """

    # Initialize variables
    (n, d) = X.shape;
    max_iter = 1000;
    initial_momentum = 0.5;
    final_momentum = 0.8;
    eta = 500;
    min_gain = 0.01;
    Y = np.random.randn(n, no_dims);
    dY = np.zeros((n, no_dims));
    iY = np.zeros((n, no_dims));
    gains = np.ones((n, no_dims));

    # Compute P-values
    P = _x2p_vandermaaten(X, 1e-5, perplexity);
    P = P + np.transpose(P);
    P = P / np.sum(P);
    P = P * 4;                                    # early exaggeration
    P = np.maximum(P, 1e-12);

    # Run iterations
    for iter in range(max_iter):

        # Compute pairwise affinities
        sum_Y = np.sum(np.square(Y), 1);
        num = 1 / (1 + np.add(np.add(-2 * np.dot(Y, Y.T), sum_Y).T, sum_Y));
        num[range(n), range(n)] = 0;
        Q = num / np.sum(num);
        Q = np.maximum(Q, 1e-12);

        # Compute gradient
        PQ = P - Q;
        for i in range(n):
            dY[i,:] = np.sum(np.tile(PQ[:,i] * num[:,i], (no_dims, 1)).T * (Y[i,:] - Y), 0);

        # Perform the update
        if iter < 20:
            momentum = initial_momentum
        else:
            momentum = final_momentum
        gains = (gains + 0.2) * ((dY > 0) != (iY > 0)) + (gains * 0.8) * ((dY > 0) == (iY > 0));
        gains[gains < min_gain] = min_gain;
        iY = momentum * iY - eta * (gains * dY);
        Y = Y + iY;
        Y = Y - np.tile(np.mean(Y, 0), (n, 1));

        # Compute current value of cost function
        if (iter + 1) % 10 == 0:
            C = np.sum(P * np.log(P / Q));
            sett.m(0,"Iteration " + str(iter + 1) + ": error is " + str(C))

        # Stop lying about P-values
        if iter == 100:
            P = P / 4;

    # Return solution
    return Y;
        
def _Hbeta_vandermaaten(D = np.array([]), beta = 1.0):
    """
    Compute the perplexity and the P-row for a specific value of the
    precision of a Gaussian distribution.
    """

    # Compute P-row and corresponding perplexity
    P = np.exp(-D.copy() * beta);
    sumP = sum(P);
    H = np.log(sumP) + beta * np.sum(D * P) / sumP;
    P = P / sumP;
    return H, P;

def _x2p_vandermaaten(X = np.array([]), tol = 1e-5, perplexity = 30.0):
    """
    Performs a binary search to get P-values in such a way that each
    conditional Gaussian has the same perplexity.
    """

    # Initialize some variables
    sett.m(0,"Computing pairwise distances...")
    (n, d) = X.shape;
    sum_X = np.sum(np.square(X), 1);
    D = np.add(np.add(-2 * np.dot(X, X.T), sum_X).T, sum_X);
    P = np.zeros((n, n));
    beta = np.ones((n, 1));
    logU = np.log(perplexity);

    # Loop over all datapoints
    for i in range(n):

        # Print progress
        if i % 500 == 0:
            sett.m(0,"Computing P-values for point ", i, " of ", n, "...")

        # Compute the Gaussian kernel and entropy for the current precision
        betamin = -np.inf;
        betamax =  np.inf;
        Di = D[i, np.concatenate((np.r_[0:i], np.r_[i+1:n]))];
        (H, thisP) = _Hbeta_vandermaaten(Di, beta[i]);

        # Evaluate whether the perplexity is within tolerance
        Hdiff = H - logU;
        tries = 0;
        while np.abs(Hdiff) > tol and tries < 50:

            # If not, increase or decrease precision
            if Hdiff > 0:
                betamin = beta[i].copy();
                if betamax == np.inf or betamax == -np.inf:
                    beta[i] = beta[i] * 2;
                else:
                    beta[i] = (beta[i] + betamax) / 2;
            else:
                betamax = beta[i].copy();
                if betamin == np.inf or betamin == -np.inf:
                    beta[i] = beta[i] / 2;
                else:
                    beta[i] = (beta[i] + betamin) / 2;

            # Recompute the values
            (H, thisP) = _Hbeta_vandermaaten(Di, beta[i]);
            Hdiff = H - logU;
            tries = tries + 1;

        # Set the final row of P
        P[i, np.concatenate((np.r_[0:i], np.r_[i+1:n]))] = thisP;

    # Return final P-matrix
    sett.m(0,"Mean value of sigma: ", np.mean(np.sqrt(1 / beta)))
    return P;

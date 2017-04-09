# coding: utf-8
# Author: F. Alex Wolf (http://falexwolf.de)
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

import numpy as np
from ..tools.pca import pca
from .. import settings as sett
from .. import plotting as plott
from .. import utils

def tsne(adata, random_state=0, n_pcs=50, perplexity=30, n_cpus=1):
    u"""
    Visualize data using t-SNE as of van der Maaten & Hinton (2008).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, optionally with metadata:
        adata['X_pca']: np.ndarray
            Result of preprocessing with PCA: observations × variables.
            If it exists, tsne will use this instead of adata.X.
    random_state : unsigned int or -1, optional (default: 0)
        Change to use different intial states for the optimization, if -1, use
        default behavior of implementation (sklearn uses np.random.seed,
        Multicore-TSNE produces a new plot at every call).
    n_pcs : int, optional (default: 50)
        Number of principal components in preprocessing PCA.
    perplexity : float, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    n_cpus : int
        Use the multicore implementation, if it is installed.

    Adds annotation
    ---------------
    X_tsne : np.ndarray of shape n_samples x 2
        Array that stores the tSNE representation of the data. Analogous
        to X_pca, X_diffmap and X_spring.
    """
    sett.mt(0,'compute tSNE')
    # preprocessing by PCA
    if 'X_pca' in adata and adata['X_pca'].shape[1] >= n_pcs:
        X = adata['X_pca']
        sett.m(0, 'using X_pca for tSNE')
    else:
        if n_pcs > 0 and adata.X.shape[1] > n_pcs:
            sett.m(0, 'preprocess using PCA with', n_pcs, 'PCs')
            sett.m(0, '--> avoid this by setting n_pcs = 0')
            X = pca(adata.X, random_state=random_state, n_comps=n_pcs)
            adata['X_pca'] = X
        else:
            X = adata.X
    # params for sklearn
    params_sklearn = {'perplexity' : perplexity,
                      'random_state': None if random_state == -1 else random_state,
                      'verbose': sett.verbosity,
                      'learning_rate': 200,
                      'early_exaggeration': 12}
    # deal with different tSNE implementations
    if n_cpus > 1:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE
            tsne = TSNE(n_jobs=2, **params_sklearn)
            sett.m(0,'... compute tSNE using MulticoreTSNE')
            Y = tsne.fit_transform(X.astype(np.float64))
        except ImportError:
            print('--> did not find package MulticoreTSNE: install it from\n'
                  '    https://github.com/DmitryUlyanov/Multicore-TSNE')
            sys.exit()
    else:    
        try:
            from sklearn.manifold import TSNE
            tsne = TSNE(**params_sklearn)
            sett.m(0,'--> can be sped up using the option `n_cpus`')
            Y = tsne.fit_transform(X)
        except ImportError:
            sett.m(0,'--> perform tSNE using slow and unreliable original\n'
                     '    implementation by L. van der Maaten\n'
                     '--> consider installing sklearn\n'
                     '    using "conda/pip install scikit-learn"')
            Y = _tsne_vandermaaten(X, 2, params['perplexity'])
    # update AnnData instance
    adata['X_tsne'] = Y
    sett.mt(0, 'finished tSNE')
    return adata

def plot_tsne(adata,
         smp=None,
         names=None,
         comps=None,
         cont=None,
         layout='2d',
         legendloc='right margin',
         cmap=None,
         pal=None,
         right_margin=None,
         size=3,
         titles=None):
    """
    Scatter plots.

    Parameters
    ----------
    dplot : dict
        Dict returned by plotting tool.
    adata : AnnData
        Annotated data matrix.
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
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    from .. import plotting as plott
    smps = plott.scatter(adata,
                    basis='tsne',
                    smp=smp,
                    names=names,
                    comps=comps,
                    cont=cont,
                    layout=layout,
                    legendloc=legendloc,
                    cmap=cmap,
                    pal=pal,
                    right_margin=right_margin,
                    size=size,
                    titles=titles)
    writekey = sett.basekey + '_tsne'
    writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
    plott.savefig(writekey)
    if not sett.savefigs and sett.autoshow:
        from ..compat.matplotlib import pyplot as pl
        pl.show()

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

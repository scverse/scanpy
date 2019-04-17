import sys
from typing import List, Tuple

import numba
import pandas as pd
import numpy as np
from numpy import linalg as la
from scipy.sparse import issparse
from anndata import AnnData

from .. import logging as logg

def design_mat(model: pd.DataFrame, batch_levels: List[str]) -> pd.DataFrame:
    """
    Computes a simple design matrix.

    At the moment, only includes the categorical annotations passed with the 'key' argument
    to the combat function

    Parameters
    --------
    model
        Contains the batch annotation
    batch_levels
        Levels of the batch annotation

    Returns
    --------
    The design matrix for the regression problem
    """
    import patsy

    design = patsy.dmatrix(
        "~ 0 + C(batch, levels={})".format(batch_levels),
        model,
        return_type="dataframe",
    )
    model = model.drop(["batch"], axis=1)
    logg.info("Found {} batches\n".format(design.shape[1]))
    other_cols = [c for i, c in enumerate(model.columns)]
    factor_matrix = model[other_cols]
    design = pd.concat((design, factor_matrix), axis=1)
    if other_cols:
        logg.info("Found {} categorical variables:".format(len(other_cols)))
        logg.info("\t" + ", ".join(other_cols) + '\n')

    return design


def stand_data(
    model: pd.DataFrame,
    data: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Standardizes the data per gene.

    The aim here is to make mean and variance be comparable across batches.

    Parameters
    --------
    model
        Contains the batch annotation
    data
        Contains the Data

    Returns
    --------
    s_data : pandas.DataFrame
        Standardized Data
    design : pandas.DataFrame
        Batch assignment as one-hot encodings
    var_pooled : numpy.ndarray
        Pooled variance per gene
    stand_mean : numpy.ndarray
        Gene-wise mean
    """

    # compute the design matrix
    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))
    drop_cols = [cname for cname, inter in ((model == 1).all()).iteritems() if inter]
    model = model[[c for c in model.columns if c not in drop_cols]]
    design = design_mat(model, batch_levels)

    # compute pooled variance estimator
    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T)
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch, :])
    var_pooled = (data - np.dot(design, B_hat).T)**2
    var_pooled = np.dot(var_pooled, np.ones((int(n_array), 1)) / int(n_array))

    # Compute the means
    if np.sum(var_pooled == 0) > 0:
        print(
            'Found {} genes with zero variance.'
            .format(np.sum(var_pooled == 0))
        )
    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, int(n_array))))
    tmp = np.array(design.copy())
    tmp[:, :n_batch] = 0
    stand_mean += np.dot(tmp, B_hat).T

    # need to be a bit careful with the zero variance genes
    # just set the zero variance genes to zero in the standardized data
    s_data = np.where(var_pooled == 0, 0, (
        (data - stand_mean) /
        np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array))))
    ))
    s_data = pd.DataFrame(s_data, index=data.index, columns=data.columns)

    return s_data, design, var_pooled, stand_mean


def combat(adata: AnnData, key: str = 'batch', inplace: bool = True):
    """ComBat function for batch effect correction [Johnson07]_ [Leek12]_ [Pedersen12]_.

    Corrects for batch effects by fitting linear models, gains statistical power
    via an EB framework where information is borrowed across genes. This uses the
    implementation of `ComBat <https://github.com/brentp/combat.py>`__ [Pedersen12]_.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix
    key: `str`, optional (default: `"batch"`)
        Key to a categorical annotation from adata.obs that will be used for batch effect removal
    inplace: bool, optional (default: `True`)
        Wether to replace adata.X or to return the corrected data

    Returns
    -------
    Depending on the value of inplace, either returns an updated AnnData object
    or modifies the passed one.
    """

    # check the input
    if key not in adata.obs.keys():
        raise ValueError('Could not find the key {!r} in adata.obs'.format(key))

    # only works on dense matrices so far
    if issparse(adata.X):
        X = adata.X.A.T
    else:
        X = adata.X.T
    data = pd.DataFrame(
        data=X,
        index=adata.var_names,
        columns=adata.obs_names,
    )

    # construct a pandas series of the batch annotation
    batch = pd.Series(adata.obs[key])
    model = pd.DataFrame({'batch': batch})
    batch_items = model.groupby("batch").groups.items()
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    # standardize across genes using a pooled variance estimator
    logg.info("Standardizing Data across genes.\n")
    s_data, design, var_pooled, stand_mean = stand_data(model, data)

    # fitting the parameters on the standardized data
    logg.info("Fitting L/S model and finding priors\n")
    batch_design = design[design.columns[:n_batch]]
    # first estimate of the additive batch effect
    gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)
    delta_hat = []

    # first estimate for the multiplicative batch effect
    for i, batch_idxs in enumerate(batch_info):
        delta_hat.append(s_data[batch_idxs].var(axis=1))

    # empirically fix the prior hyperparameters
    gamma_bar = gamma_hat.mean(axis=1)
    t2 = gamma_hat.var(axis=1)
    # a_prior and b_prior are the priors on lambda and theta from Johnson and Li (2006)
    a_prior = list(map(aprior, delta_hat))
    b_prior = list(map(bprior, delta_hat))

    logg.info("Finding parametric adjustments\n")
    # gamma star and delta star will be our empirical bayes (EB) estimators
    # for the additive and multiplicative batch effect per batch and cell
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info):
        # temp stores our estimates for the batch effect parameters.
        # temp[0] is the additive batch effect
        # temp[1] is the multiplicative batch effect
        gamma, delta = _it_sol(
            s_data[batch_idxs].values,
            gamma_hat[i],
            delta_hat[i].values,
            gamma_bar[i],
            t2[i],
            a_prior[i],
            b_prior[i],
        )

        gamma_star.append(gamma)
        delta_star.append(delta)

    logg.info("Adjusting data\n")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)

    # we now apply the parametric adjustment to the standardized data from above
    # loop over all batches in the data
    for j, batch_idxs in enumerate(batch_info):

        # we basically substract the additive batch effect, rescale by the ratio
        # of multiplicative batch effect to pooled variance and add the overall gene
        # wise mean
        dsq = np.sqrt(delta_star[j,:])
        dsq = dsq.reshape((len(dsq), 1))
        denom =  np.dot(dsq, np.ones((1, n_batches[j])))
        numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.loc[batch_idxs], gamma_star).T)
        bayesdata[batch_idxs] = numer / denom

    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    bayesdata = bayesdata * np.dot(vpsq, np.ones((1, int(n_array)))) + stand_mean

    # put back into the adata object or return
    if inplace:
        adata.X = bayesdata.values.transpose()
    else:
        return bayesdata.values.transpose()


@numba.jit
def _it_sol(s_data, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001) -> Tuple[float, float]:
    """
    Iteratively compute the conditional posterior means for gamma and delta.

    gamma is an estimator for the additive batch effect, deltat is an estimator
    for the multiplicative batch effect. We use an EB framework to estimate these
    two. Analytical expressions exist for both parameters, which however depend on each other.
    We therefore iteratively evalutate these two expressions until convergence is reached.

    Parameters
    --------
    s_data : pd.DataFrame
        Contains the standardized Data
    g_hat : float
        Initial guess for gamma
    d_hat : float
        Initial guess for delta
    g_bar, t_2, a, b : float
        Hyperparameters
    conv: float, optional (default: `0.0001`)
        convergence criterium

    Returns:
    --------
    gamma : float
        estimated value for gamma
    delta : float
        estimated value for delta
    """

    n = (1 - np.isnan(s_data)).sum(axis=1)
    g_old = g_hat.copy()
    d_old = d_hat.copy()

    change = 1
    count = 0

    # we place a normally distributed prior on gamma and and inverse gamma prior on delta
    # in the loop, gamma and delta are updated together. they depend on each other. we iterate until convergence.
    while change > conv:
        g_new = (t2*n*g_hat + d_old*g_bar) / (t2*n + d_old)
        sum2 = s_data - g_new.reshape((g_new.shape[0], 1)) @ np.ones((1, s_data.shape[1]))
        sum2 = sum2 ** 2
        sum2 = sum2.sum(axis=1)
        d_new = (0.5*sum2 + b) / (n/2.0 + a-1.0)

        change = max((abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max())
        g_old = g_new  # .copy()
        d_old = d_new  # .copy()
        count = count + 1

    return g_new, d_new


def aprior(delta_hat):
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (2*s2 + m**2) / s2


def bprior(delta_hat):
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (m*s2 + m**3) / s2

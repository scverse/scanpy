from typing import Optional, Dict
from warnings import warn

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse

from scanpy import logging as logg

from scanpy._utils import view_to_actual, check_nonnegative_integers
from scanpy.get import _get_obs_rep, _set_obs_rep

from scanpy.preprocessing._pca import pca


def _pearson_residuals(X, theta, clip, check_values, copy=False):

    X = X.copy() if copy else X

    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        raise ValueError('Pearson residuals require theta > 0')
    # prepare clipping
    if clip is None:
        n = X.shape[0]
        clip = np.sqrt(n)
    if clip < 0:
        raise ValueError("Pearson residuals require `clip>=0` or `clip=None`.")

    if check_values and not check_nonnegative_integers(X):
        warn(
            "`normalize_pearson_residuals()` expects raw count data, but non-integers were found.",
            UserWarning,
        )

    if issparse(X):
        sums_genes = np.sum(X, axis=0)
        sums_cells = np.sum(X, axis=1)
        sum_total = np.sum(sums_genes).squeeze()
    else:
        sums_genes = np.sum(X, axis=0, keepdims=True)
        sums_cells = np.sum(X, axis=1, keepdims=True)
        sum_total = np.sum(sums_genes)

    mu = np.array(sums_cells @ sums_genes / sum_total)
    diff = np.array(X - mu)
    residuals = diff / np.sqrt(mu + mu ** 2 / theta)

    # clip
    residuals = np.clip(residuals, a_min=-clip, a_max=clip)

    return residuals


def normalize_pearson_residuals(
    adata: AnnData,
    theta: float = 100,
    clip: Optional[float] = None,
    layer: Optional[str] = None,
    copy: bool = False,
    check_values: bool = True,
    inplace: bool = True,
) -> Optional[Dict[str, np.ndarray]]:
    """\
    Computes analytic Pearson residuals, based on [Lause20]_.

    Assuming a negative binomial offset model with overdispersion
    theta shared across genes, computes Pearson residuals. By default, residuals
    are clipped to sqrt(n) and overdispersion theta=100 is used.

    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    theta
        The NB overdispersion parameter theta. Higher values correspond to
        less overdispersion (var = mean + mean^2/theta), and `theta=np.Inf`
        corresponds to a Poisson model.
    clip
        Determines if and how residuals are clipped:

            * If `None`, residuals are clipped to the interval \
            [-sqrt(n), sqrt(n)], where n is the number of cells in the dataset (default behavior).
            * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
            `clip=np.Inf` for no clipping.

    layer
        Layer to normalize instead of `X`. If `None`, `X` is normalized.
    copy
        Whether to modify copied input object. Not compatible with
        `inplace=False`.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to True.
    inplace
        Whether to update `adata` or return dictionary with normalized copies
        of `adata.X` and `adata.layers`.

    Returns
    -------
    Returns dictionary with Pearson residuals and settings
    or updates `adata` with normalized version of the original
    `adata.X` and `adata.layers`, depending on `inplace`.

    """

    if copy:
        if not inplace:
            raise ValueError("`copy=True` cannot be used with `inplace=False`.")
        adata = adata.copy()

    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer)
    computed_on = layer if layer else 'adata.X'

    msg = f'computing analytic Pearson residuals on {computed_on}'
    start = logg.info(msg)

    residuals = _pearson_residuals(X, theta, clip, check_values, copy=~inplace)
    settings_dict = dict(theta=theta, clip=clip, computed_on=computed_on)

    if inplace:
        _set_obs_rep(adata, residuals, layer=layer)
        adata.uns['pearson_residuals_normalization'] = settings_dict
    else:
        results_dict = dict(X=residuals, **settings_dict)

    logg.info('    finished ({time_passed})', time=start)

    if copy:
        return adata
    elif not inplace:
        return results_dict


def normalize_pearson_residuals_pca(
    adata: AnnData,
    theta: float = 100,
    clip: Optional[float] = None,
    n_comps_pca: Optional[int] = 50,
    random_state_pca: Optional[float] = 0,
    use_highly_variable: bool = True,
    kwargs_pca: Optional[dict] = {},
    check_values: bool = True,
    inplace: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    Applies Pearson residual normalization and PCA, based on [Lause20]_.

    Operates on the subset of highly variable genes in `adata.var['highly_variable']` by default.

    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    use_highly_variable
        Whether to use the gene selection in `adata.var['highly_variable']` to
        subset the data before normalizing (default) or proceed on the full
        dataset.
    theta
        This is the NB overdispersion parameter theta for Pearson residual
        computations. Higher values correspond to less overdispersion
        (var = mean + mean^2/theta), and `theta=np.Inf` corresponds to a
        Poisson model.
    clip
        This determines how Pearson residuals are clipped:

            * If `None`, residuals are clipped to the interval \
            [-sqrt(n), sqrt(n)], where n is the number of cells in the dataset (default behavior).
            * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
            `clip=np.Inf` for no clipping.

    n_comps_pca
        Number of principal components to compute.
    random_state_pca
        Change to use different initial states for the optimization.
    kwargs_pca
        Dictionary of further keyword arguments passed on to `scanpy.pp.pca()`.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to True.
    inplace
        Whether to place results in `adata` or return them.


    Returns
    -------
    If `inplace=False`, returns the Pearson residual-based PCA results
    (`adata_pca`).
    If `inplace=True`, updates `adata` with the following fields:

    `.uns['pearson_residuals_normalization']['pearson_residuals_df']`
         The hvg-subset, normalized by Pearson residuals
    `.uns['pearson_residuals_normalization']['theta']`
         The used value of the overdisperion parameter theta
    `.uns['pearson_residuals_normalization']['clip']`
         The used value of the clipping parameter

    `.obsm['X_pca']`
        PCA representation of data after gene selection and Pearson residual
        normalization.
    `.uns['pca']['PCs']`
         The principal components containing the loadings.
    `.uns['pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pca']['variance']`
         Explained variance, equivalent to the eigenvalues of the
         covariance matrix.

    """

    if use_highly_variable and 'highly_variable' in adata.var_keys():
        # TODO: are these copies needed?
        adata_pca = adata[:, adata.var['highly_variable']].copy()
    else:
        # TODO: are these copies needed?
        adata_pca = adata.copy()

    normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps_pca, random_state=random_state_pca, **kwargs_pca)

    if inplace:
        norm_settings = adata_pca.uns['pearson_residuals_normalization']
        norm_dict = dict(**norm_settings, pearson_residuals_df=adata_pca.to_df())
        pca_settings = adata_pca.uns['pca']
        pca_dict = dict(**pca_settings, PCs=adata_pca.varm['PCs'])
        adata.uns['pca'] = pca_dict
        adata.uns['pearson_residuals_normalization'] = norm_dict
        adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
        return None
    else:
        return adata_pca

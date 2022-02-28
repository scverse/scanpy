from typing import Optional, Dict
from warnings import warn

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse

from scanpy import logging as logg

from scanpy._utils import view_to_actual, check_nonnegative_integers
from scanpy.get import _get_obs_rep, _set_obs_rep
from scanpy._utils import _doc_params
from scanpy.preprocessing._pca import pca
from scanpy.experimental._docs import (
    doc_adata,
    doc_dist_params,
    doc_layer,
    doc_copy,
    doc_inplace,
)


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
    residuals = diff / np.sqrt(mu + mu**2 / theta)

    # clip
    residuals = np.clip(residuals, a_min=-clip, a_max=clip)

    return residuals


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    layer=doc_layer,
    inplace=doc_inplace,
    copy=doc_copy,
)
def normalize_pearson_residuals(
    adata: AnnData,
    *,
    theta: float = 100,
    clip: Optional[float] = None,
    check_values: bool = True,
    layer: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Dict[str, np.ndarray]]:
    """\
    Applies analytic Pearson residual normalization, based on [Lause21]_.

    The residuals are based on a negative binomial offset model with overdispersion
    `theta` shared across genes. By default, residuals are clipped to sqrt(n) and
    overdispersion `theta=100` is used.

    Expects raw count input.

    Params
    ------
    {adata}
    {dist_params}
    {layer}
    {inplace}
    {copy}

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
    *,
    theta: float = 100,
    clip: Optional[float] = None,
    n_comps: Optional[int] = 50,
    random_state: Optional[float] = 0,
    kwargs_pca: Optional[dict] = {},
    use_highly_variable: Optional[bool] = None,
    check_values: bool = True,
    inplace: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    Applies analytic Pearson residual normalization and PCA, based on [Lause21]_.

    The residuals are based on a negative binomial offset model with overdispersion
    `theta` shared across genes. By default, residuals are clipped to sqrt(n),
    overdispersion `theta=100` is used, and PCA is run with 50 components.

    Operates on the subset of highly variable genes in `adata.var['highly_variable']`
    by default. Expects raw count input.


    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    theta
        The negative binomial overdispersion parameter theta for Pearson residuals.
        Higher values correspond to less overdispersion (var = mean + mean^2/theta),
        and `theta=np.Inf` corresponds to a Poisson model.
    clip
        Determines if and how residuals are clipped:

            * If `None`, residuals are clipped to the interval [-sqrt(n), sqrt(n)], \
            where n is the number of cells in the dataset (default behavior).
            * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
            `clip=np.Inf` for no clipping.

    n_comps
        Number of principal components to compute for the PCA step.
    random_state
        Change to use different initial states for the optimization of the PCA step.
    kwargs_pca
        Dictionary of further keyword arguments passed on to `scanpy.pp.pca()`.
    use_highly_variable
        Whether to use the gene selection in `adata.var['highly_variable']` to subset
        the data before normalizing (default) or proceed on the full dataset.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to
        True.
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
        PCA representation of data after gene selection (if applicable) and Pearson
        residual normalization.
    `.varm['PCs']`
         The principal components containing the loadings. When `inplace=True` and
         `use_highly_variable=True`, this will contain empty rows for the genes not
         selected.
    `.uns['pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pca']['variance']`
         Explained variance, equivalent to the eigenvalues of the covariance matrix.

    """

    # check if HVG selection is there if user wants to use it
    if use_highly_variable and 'highly_variable' not in adata.var_keys():
        raise ValueError(
            "You passed `use_highly_variable=True`, but no HVG selection was found "
            "(e.g., there was no 'highly_variable' column in adata.var).'"
        )

    # default behavior: if there is a HVG selection, we will use it
    if use_highly_variable is None and 'highly_variable' in adata.var_keys():
        use_highly_variable = True

    if use_highly_variable:
        adata_sub = adata[:, adata.var['highly_variable']].copy()
        adata_pca = AnnData(
            adata_sub.X.copy(), obs=adata_sub.obs[[]], var=adata_sub.var[[]]
        )
    else:
        adata_pca = AnnData(adata.X.copy(), obs=adata.obs[[]], var=adata.var[[]])

    normalize_pearson_residuals(
        adata_pca, theta=theta, clip=clip, check_values=check_values
    )
    pca(adata_pca, n_comps=n_comps, random_state=random_state, **kwargs_pca)

    if inplace:
        norm_settings = adata_pca.uns['pearson_residuals_normalization']
        norm_dict = dict(**norm_settings, pearson_residuals_df=adata_pca.to_df())
        if use_highly_variable:
            adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm['PCs'][adata.var['highly_variable']] = adata_pca.varm['PCs']
        else:
            adata.varm['PCs'] = adata_pca.varm['PCs']
        adata.uns['pca'] = adata_pca.uns['pca']
        adata.uns['pearson_residuals_normalization'] = norm_dict
        adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
        return None
    else:
        return adata_pca

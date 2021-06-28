from typing import Optional, Union, Iterable, Dict
from warnings import warn

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs

from scanpy import logging as logg
from scanpy._compat import Literal

from scanpy.preprocessing._utils import view_to_actual, check_nonnegative_integers
from scanpy.get import _get_obs_rep, _set_obs_rep

from ._pca import pca


def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    if issubclass(X.dtype.type, (int, np.integer)):
        X = X.astype(np.float32)  # TODO: Check if float64 should be used
    counts = np.asarray(counts)  # dask doesn't do medians
    after = np.median(counts[counts > 0], axis=0) if after is None else after
    counts += counts == 0
    counts = counts / after
    if issparse(X):
        sparsefuncs.inplace_row_scale(X, 1 / counts)
    else:
        np.divide(X, counts[:, None], out=X)
    return X


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

    if check_values and (check_nonnegative_integers(X) is False):
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
    Computes analytic Pearson residuals, assuming a negative binomial offset
    model with overdispersion theta shared across genes. By default, residuals
    are clipped to sqrt(n) and overdispersion theta=100 is used.

    Based on "Analytic Pearson residuals for normalization of single-cell
    RNA-seq UMI data", bioRxiv, [Lause20]_.

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

        * If `None`, residuals are clipped to the interval
        [-sqrt(n), sqrt(n)], where n is the number of cells in the dataset
        (default behavior).
        * If any scalar c, residuals are clipped to the interval [-c, c]. Set
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

    msg = 'computing analytic Pearson residuals on %s' % computed_on
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
    Applies PCA based on Pearson residual normalization. Operates on the
    subset of highly variable genes in `adata.var['highly_variable']` by
    default.

    This workflow is based on "Analytic Pearson residuals for normalization of
    single-cell RNA-seq UMI data", bioRxiv, [Lause20]_.


    Parameters
    ----------
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
        This determines if and how Pearson residuals are clipped:

        * If `None`, residuals are clipped to the interval
        [-sqrt(n), sqrt(n)], where n is the number of cells in the dataset
        (default behavior).
        * If any scalar c, residuals are clipped to the interval [-c, c]. Set
        `clip=np.Inf` for no clipping.

    n_comps_pca
        Number of principal components to compute.
    random_state_pca
        Change to use different initial states for the optimization.
    kwargs_pca
        Dictionary of further keyword arguments passed on to `sc.pp.pca()`.
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

    `.obsm['X_pearson_residuals_pca']`
        PCA representation of data after gene selection and Pearson residual
        normalization.
    `.uns['pearson_residuals_pca']['PCs']`
         The principal components containing the loadings.
    `.uns['pearson_residuals_pca']['variance_ratio']`
         Ratio of explained variance.
    `.uns['pearson_residuals_pca']['variance']`
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
        adata.uns['pearson_residuals_pca'] = pca_dict
        adata.uns['pearson_residuals_normalization'] = norm_dict
        adata.obsm['X_pearson_residuals_pca'] = adata_pca.obsm['X_pca']
        return None
    else:
        adata_pca.obsm['X_pearson_residuals_pca'] = adata_pca.obsm['X_pca'].copy()
        adata_pca.uns['pearson_residuals_pca'] = adata_pca.uns['pca'].copy()
        del adata_pca.obsm['X_pca']
        del adata_pca.uns['pca']
        return adata_pca


def normalize_total(
    adata: AnnData,
    target_sum: Optional[float] = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    layers: Union[Literal['all'], Iterable[str]] = None,
    layer_norm: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Dict[str, np.ndarray]]:
    """\
    Normalize counts per cell.

    Normalize each cell by total counts over all genes,
    so that every cell has the same total count after normalization.
    If choosing `target_sum=1e6`, this is CPM normalization.

    If `exclude_highly_expressed=True`, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes [Weinreb17]_.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    target_sum
        If `None`, after normalization, each observation (cell) has a total
        count equal to the median of total counts for observations (cells)
        before normalization.
    exclude_highly_expressed
        Exclude (very) highly expressed genes for the computation of the
        normalization factor (size factor) for each cell. A gene is considered
        highly expressed, if it has more than `max_fraction` of the total counts
        in at least one cell. The not-excluded genes will sum up to
        `target_sum`.
    max_fraction
        If `exclude_highly_expressed=True`, consider cells as highly expressed
        that have more counts than `max_fraction` of the original total counts
        in at least one cell.
    key_added
        Name of the field in `adata.obs` where the normalization factor is
        stored.
    layer
        Layer to normalize instead of `X`. If `None`, `X` is normalized.
    inplace
        Whether to update `adata` or return dictionary with normalized copies of
        `adata.X` and `adata.layers`.
    copy
        Whether to modify copied input object. Not compatible with inplace=False.

    Returns
    -------
    Returns dictionary with normalized copies of `adata.X` and `adata.layers`
    or updates `adata` with normalized version of the original
    `adata.X` and `adata.layers`, depending on `inplace`.

    Example
    --------
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> sc.settings.verbosity = 2
    >>> np.set_printoptions(precision=2)
    >>> adata = AnnData(np.array([
    ...    [3, 3, 3, 6, 6],
    ...    [1, 1, 1, 2, 2],
    ...    [1, 22, 1, 2, 2],
    ... ]))
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)['X']
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(
    ...     adata, target_sum=1, exclude_highly_expressed=True,
    ...     max_fraction=0.2, inplace=False
    ... )['X']
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)
    """
    if copy:
        if not inplace:
            raise ValueError("`copy=True` cannot be used with `inplace=False`.")
        adata = adata.copy()

    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('Choose max_fraction between 0 and 1.')

    # Deprecated features
    if layers is not None:
        warn(
            FutureWarning(
                "The `layers` argument is deprecated. Instead, specify individual "
                "layers to normalize with `layer`."
            )
        )
    if layer_norm is not None:
        warn(
            FutureWarning(
                "The `layer_norm` argument is deprecated. Specify the target size "
                "factor directly with `target_sum`."
            )
        )

    if layers == 'all':
        layers = adata.layers.keys()
    elif isinstance(layers, str):
        raise ValueError(
            f"`layers` needs to be a list of strings or 'all', not {layers!r}"
        )

    view_to_actual(adata)

    X = _get_obs_rep(adata, layer=layer)

    gene_subset = None
    msg = 'normalizing counts per cell'
    if exclude_highly_expressed:
        counts_per_cell = X.sum(1)  # original counts per cell
        counts_per_cell = np.ravel(counts_per_cell)

        # at least one cell as more than max_fraction of counts per cell

        gene_subset = (X > counts_per_cell[:, None] * max_fraction).sum(0)
        gene_subset = np.ravel(gene_subset) == 0

        msg += (
            ' The following highly-expressed genes are not considered during '
            f'normalization factor computation:\n{adata.var_names[~gene_subset].tolist()}'
        )
        counts_per_cell = X[:, gene_subset].sum(1)
    else:
        counts_per_cell = X.sum(1)
    start = logg.info(msg)
    counts_per_cell = np.ravel(counts_per_cell)

    cell_subset = counts_per_cell > 0
    if not np.all(cell_subset):
        warn(UserWarning('Some cells have zero counts'))

    if inplace:
        if key_added is not None:
            adata.obs[key_added] = counts_per_cell
        _set_obs_rep(
            adata, _normalize_data(X, counts_per_cell, target_sum), layer=layer
        )
    else:
        # not recarray because need to support sparse
        dat = dict(
            X=_normalize_data(X, counts_per_cell, target_sum, copy=True),
            norm_factor=counts_per_cell,
        )

    # Deprecated features
    if layer_norm == 'after':
        after = target_sum
    elif layer_norm == 'X':
        after = np.median(counts_per_cell[cell_subset])
    elif layer_norm is None:
        after = None
    else:
        raise ValueError('layer_norm should be "after", "X" or None')

    for layer_to_norm in layers if layers is not None else ():
        res = normalize_total(
            adata, layer=layer_to_norm, target_sum=after, inplace=inplace
        )
        if not inplace:
            dat[layer_to_norm] = res["X"]

    logg.info(
        '    finished ({time_passed})',
        time=start,
    )
    if key_added is not None:
        logg.debug(
            f'and added {key_added!r}, counts per cell before normalization (adata.obs)'
        )

    if copy:
        return adata
    elif not inplace:
        return dat

import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata import AnnData


from .. import logging as logg
from .._settings import settings, Verbosity
from .._utils import sanitize_anndata, check_nonnegative_integers, view_to_actual
from scanpy.get import _get_obs_rep, _set_obs_rep
from .._compat import Literal
from ._utils import _get_mean_var
from ._distributed import materialize_as_ndarray
from ._simple import filter_genes


def _highly_variable_genes_seurat_v3(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: int = 2000,
    batch_key: Optional[str] = None,
    check_values: bool = True,
    span: float = 0.3,
    subset: bool = False,
    inplace: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    See `highly_variable_genes`.

    For further implemenation details see https://www.overleaf.com/read/ckptrbgzzzpg

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pd.DataFrame`) or
    updates `.var` with the following fields:

    highly_variable : bool
        boolean indicator of highly-variable genes.
    **means**
        means per gene.
    **variances**
        variance per gene.
    **variances_norm**
        normalized variance per gene, averaged in the case of multiple batches.
    highly_variable_rank : float
        Rank of the gene according to normalized variance, median rank in the case of multiple batches.
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as HVG.
    """

    try:
        from skmisc.loess import loess
    except ImportError:
        raise ImportError(
            'Please install skmisc package via `pip install --user scikit-misc'
        )
    df = pd.DataFrame(index=adata.var_names)
    X = adata.layers[layer] if layer is not None else adata.X

    if check_values and not check_nonnegative_integers(X):
        warnings.warn(
            "`flavor='seurat_v3'` expects raw count data, but non-integers were found.",
            UserWarning,
        )

    df['means'], df['variances'] = _get_mean_var(X)

    if batch_key is None:
        batch_info = pd.Categorical(np.zeros(adata.shape[0], dtype=int))
    else:
        batch_info = adata.obs[batch_key].values

    norm_gene_vars = []
    for b in np.unique(batch_info):
        X_batch = X[batch_info == b]

        mean, var = _get_mean_var(X_batch)
        not_const = var > 0
        estimat_var = np.zeros(X.shape[1], dtype=np.float64)

        y = np.log10(var[not_const])
        x = np.log10(mean[not_const])
        model = loess(x, y, span=span, degree=2)
        model.fit()
        estimat_var[not_const] = model.outputs.fitted_values
        reg_std = np.sqrt(10 ** estimat_var)

        batch_counts = X_batch.astype(np.float64).copy()
        # clip large values as in Seurat
        N = X_batch.shape[0]
        vmax = np.sqrt(N)
        clip_val = reg_std * vmax + mean
        if sp_sparse.issparse(batch_counts):
            batch_counts = sp_sparse.csr_matrix(batch_counts)
            mask = batch_counts.data > clip_val[batch_counts.indices]
            batch_counts.data[mask] = clip_val[batch_counts.indices[mask]]

            squared_batch_counts_sum = np.array(batch_counts.power(2).sum(axis=0))
            batch_counts_sum = np.array(batch_counts.sum(axis=0))
        else:
            clip_val_broad = np.broadcast_to(clip_val, batch_counts.shape)
            np.putmask(
                batch_counts,
                batch_counts > clip_val_broad,
                clip_val_broad,
            )

            squared_batch_counts_sum = np.square(batch_counts).sum(axis=0)
            batch_counts_sum = batch_counts.sum(axis=0)

        norm_gene_var = (1 / ((N - 1) * np.square(reg_std))) * (
            (N * np.square(mean))
            + squared_batch_counts_sum
            - 2 * batch_counts_sum * mean
        )
        norm_gene_vars.append(norm_gene_var.reshape(1, -1))

    norm_gene_vars = np.concatenate(norm_gene_vars, axis=0)
    # argsort twice gives ranks, small rank means most variable
    ranked_norm_gene_vars = np.argsort(np.argsort(-norm_gene_vars, axis=1), axis=1)

    # this is done in SelectIntegrationFeatures() in Seurat v3
    ranked_norm_gene_vars = ranked_norm_gene_vars.astype(np.float32)
    num_batches_high_var = np.sum(
        (ranked_norm_gene_vars < n_top_genes).astype(int), axis=0
    )
    ranked_norm_gene_vars[ranked_norm_gene_vars >= n_top_genes] = np.nan
    ma_ranked = np.ma.masked_invalid(ranked_norm_gene_vars)
    median_ranked = np.ma.median(ma_ranked, axis=0).filled(np.nan)

    df['highly_variable_nbatches'] = num_batches_high_var
    df['highly_variable_rank'] = median_ranked
    df['variances_norm'] = np.mean(norm_gene_vars, axis=0)

    sorted_index = (
        df[['highly_variable_rank', 'highly_variable_nbatches']]
        .sort_values(
            ['highly_variable_rank', 'highly_variable_nbatches'],
            ascending=[True, False],
            na_position='last',
        )
        .index
    )
    df['highly_variable'] = False
    df.loc[sorted_index[: int(n_top_genes)], 'highly_variable'] = True

    if inplace or subset:
        adata.uns['hvg'] = {'flavor': 'seurat_v3'}
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'highly_variable_rank\', float vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'variances\', float vector (adata.var)\n'
            '    \'variances_norm\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = df['highly_variable'].values
        adata.var['highly_variable_rank'] = df['highly_variable_rank'].values
        adata.var['means'] = df['means'].values
        adata.var['variances'] = df['variances'].values
        adata.var['variances_norm'] = df['variances_norm'].values.astype(
            'float64', copy=False
        )
        if batch_key is not None:
            adata.var['highly_variable_nbatches'] = df[
                'highly_variable_nbatches'
            ].values
        if subset:
            adata._inplace_subset_var(df['highly_variable'].values)
    else:
        if batch_key is None:
            df = df.drop(['highly_variable_nbatches'], axis=1)
        return df


def _highly_variable_pearson_residuals(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: int = 1000,
    batch_key: Optional[str] = None,
    theta: float = 100,
    clip: Optional[float] = None,
    chunksize: int = 100,
    check_values: bool = True,
    subset: bool = False,
    inplace: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    See `highly_variable_genes`.

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pd.DataFrame`)
    or updates `.var` with the following fields:

    highly_variable : bool
        boolean indicator of highly-variable genes.
    means : float
        means per gene.
    variances : float
        variances per gene.
    residual_variances : float
        Pearson residual variance per gene. Averaged in the case of multiple
        batches.
    highly_variable_rank : float
        Rank of the gene according to residual variance, median rank in the
        case of multiple batches. NaN for non-HVGs.
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are
        detected as HVG.
    highly_variable_intersection : bool
        If batch_key is given, this denotes the genes that are highly variable
        in all batches.
    """

    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer)
    computed_on = layer if layer else 'adata.X'

    # Check for raw counts
    if check_values and (check_nonnegative_integers(X) is False):
        warnings.warn(
            "`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
            UserWarning,
        )
    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        raise ValueError('Pearson residuals require theta > 0')
    # prepare clipping

    if batch_key is None:
        batch_info = np.zeros(adata.shape[0], dtype=int)
    else:
        batch_info = adata.obs[batch_key].values
    n_batches = len(np.unique(batch_info))

    # Get pearson residuals for each batch separately
    residual_gene_vars = []
    for batch in np.unique(batch_info):

        adata_subset = adata[batch_info == batch]

        # Filter out zero genes
        with settings.verbosity.override(Verbosity.error):
            nonzero_genes = filter_genes(adata_subset, min_cells=1, inplace=False)[0]
        adata_subset = adata_subset[:, nonzero_genes]

        if layer is not None:
            X_batch = adata_subset.layers[layer]
        else:
            X_batch = adata_subset.X

        # Prepare clipping
        if clip is None:
            n = X_batch.shape[0]
            clip = np.sqrt(n)
        if clip < 0:
            raise ValueError("Pearson residuals require `clip>=0` or `clip=None`.")

        if sp_sparse.issparse(X_batch):
            sums_genes = np.sum(X_batch, axis=0)
            sums_cells = np.sum(X_batch, axis=1)
            sum_total = np.sum(sums_genes).squeeze()
        else:
            sums_genes = np.sum(X_batch, axis=0, keepdims=True)
            sums_cells = np.sum(X_batch, axis=1, keepdims=True)
            sum_total = np.sum(sums_genes)

        # Compute pearson residuals in chunks
        residual_gene_var = np.empty((X_batch.shape[1]))
        for start in np.arange(0, X_batch.shape[1], chunksize):
            stop = start + chunksize
            mu = np.array(sums_cells @ sums_genes[:, start:stop] / sum_total)
            X_dense = X_batch[:, start:stop].toarray()
            residuals = (X_dense - mu) / np.sqrt(mu + mu ** 2 / theta)
            residuals = np.clip(residuals, a_min=-clip, a_max=clip)
            residual_gene_var[start:stop] = np.var(residuals, axis=0)

        # Add 0 values for genes that were filtered out
        zero_gene_var = np.zeros(np.sum(~nonzero_genes))
        residual_gene_var = np.concatenate((residual_gene_var, zero_gene_var))
        # Order as before filtering
        idxs = np.concatenate((np.where(nonzero_genes)[0], np.where(~nonzero_genes)[0]))
        residual_gene_var = residual_gene_var[np.argsort(idxs)]
        residual_gene_vars.append(residual_gene_var.reshape(1, -1))

    residual_gene_vars = np.concatenate(residual_gene_vars, axis=0)

    # Get cutoffs and define hvgs per batch
    residual_gene_vars_sorted = np.sort(residual_gene_vars, axis=1)
    cutoffs_per_batch = residual_gene_vars_sorted[:, -n_top_genes]
    highly_variable_per_batch = np.greater_equal(
        residual_gene_vars.T, cutoffs_per_batch
    ).T

    # Merge hvgs across batches
    highly_variable_nbatches = np.sum(highly_variable_per_batch, axis=0)
    highly_variable_intersection = highly_variable_nbatches == n_batches

    # Get rank per gene within each batch
    # argsort twice gives ranks, small rank means most variable
    ranks_residual_var = np.argsort(np.argsort(-residual_gene_vars, axis=1), axis=1)
    ranks_residual_var = ranks_residual_var.astype(np.float32)
    ranks_residual_var[ranks_residual_var >= n_top_genes] = np.nan
    ranks_masked_array = np.ma.masked_invalid(ranks_residual_var)
    # Median rank across batches,
    # ignoring batches in which gene was not selected
    medianrank_residual_var = np.ma.median(ranks_masked_array, axis=0).filled(np.nan)

    means, variances = materialize_as_ndarray(_get_mean_var(X))
    df = pd.DataFrame.from_dict(
        dict(
            means=means,
            variances=variances,
            residual_variances=np.mean(residual_gene_vars, axis=0).astype(
                np.float32, copy=False
            ),
            highly_variable_rank=medianrank_residual_var,
            highly_variable_nbatches=highly_variable_nbatches.astype(np.int64),
            highly_variable_intersection=highly_variable_intersection,
        )
    )
    df = df.set_index(adata.var_names)

    # Sort genes by how often they selected as hvg within each batch and
    # break ties with median rank of residual variance across batches
    df.sort_values(
        ['highly_variable_nbatches', 'highly_variable_rank'],
        ascending=[False, True],
        na_position='last',
        inplace=True,
    )
    df['highly_variable'] = False
    df.highly_variable.iloc[:n_top_genes] = True
    # TODO: following line raises a pandas warning
    # (also for flavor = seurat and cellranger..)
    df = df.loc[adata.var_names]

    if inplace:
        adata.uns['hvg'] = {'flavor': 'pearson_residuals', 'computed_on': computed_on}
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'highly_variable_rank\', float vector (adata.var)\n'
            '    \'highly_variable_nbatches\', int vector (adata.var)\n'
            '    \'highly_variable_intersection\', boolean vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'variances\', float vector (adata.var)\n'
            '    \'residual_variances\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = df['highly_variable'].values
        adata.var['highly_variable_rank'] = df['highly_variable_rank'].values
        adata.var['means'] = df['means'].values
        adata.var['variances'] = df['variances'].values
        adata.var['residual_variances'] = df['residual_variances']

        if batch_key is not None:
            adata.var['highly_variable_nbatches'] = df[
                'highly_variable_nbatches'
            ].values
            adata.var['highly_variable_intersection'] = df[
                'highly_variable_intersection'
            ].values
        if subset:
            adata._inplace_subset_var(df['highly_variable'].values)
    else:
        if batch_key is None:
            df = df.drop(
                ['highly_variable_nbatches', 'highly_variable_intersection'], axis=1
            )
        if subset:
            df = df.iloc[df.highly_variable.values, :]

        return df


def _highly_variable_genes_single_batch(
    adata: AnnData,
    layer: Optional[str] = None,
    min_disp: Optional[float] = 0.5,
    max_disp: Optional[float] = np.inf,
    min_mean: Optional[float] = 0.0125,
    max_mean: Optional[float] = 3,
    n_top_genes: Optional[int] = None,
    n_bins: int = 20,
    flavor: Literal['seurat', 'cell_ranger'] = 'seurat',
) -> pd.DataFrame:
    """\
    See `highly_variable_genes`.

    Returns
    -------
    A DataFrame that contains the columns
    `highly_variable`, `means`, `dispersions`, and `dispersions_norm`.
    """
    X = adata.layers[layer] if layer is not None else adata.X
    if flavor == 'seurat':
        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            X *= np.log(adata.uns['log1p']['base'])
        X = np.expm1(X)

    mean, var = materialize_as_ndarray(_get_mean_var(X))
    # now actually compute the dispersion
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    if flavor == 'seurat':  # logarithmized mean as in Seurat
        dispersion[dispersion == 0] = np.nan
        dispersion = np.log(dispersion)
        mean = np.log1p(mean)
    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['means'] = mean
    df['dispersions'] = dispersion
    if flavor == 'seurat':
        df['mean_bin'] = pd.cut(df['means'], bins=n_bins)
        disp_grouped = df.groupby('mean_bin')['dispersions']
        disp_mean_bin = disp_grouped.mean()
        disp_std_bin = disp_grouped.std(ddof=1)
        # retrieve those genes that have nan std, these are the ones where
        # only a single gene fell in the bin and implicitly set them to have
        # a normalized disperion of 1
        one_gene_per_bin = disp_std_bin.isnull()
        gen_indices = np.where(one_gene_per_bin[df['mean_bin'].values])[0].tolist()
        if len(gen_indices) > 0:
            logg.debug(
                f'Gene indices {gen_indices} fell into a single bin: their '
                'normalized dispersion was set to 1.\n    '
                'Decreasing `n_bins` will likely avoid this effect.'
            )
        # Circumvent pandas 0.23 bug. Both sides of the assignment have dtype==float32,
        # but there’s still a dtype error without “.value”.
        disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[
            one_gene_per_bin.values
        ].values
        disp_mean_bin[one_gene_per_bin.values] = 0
        # actually do the normalization
        df['dispersions_norm'] = (
            df['dispersions'].values  # use values here as index differs
            - disp_mean_bin[df['mean_bin'].values].values
        ) / disp_std_bin[df['mean_bin'].values].values
    elif flavor == 'cell_ranger':
        from statsmodels import robust

        df['mean_bin'] = pd.cut(
            df['means'],
            np.r_[-np.inf, np.percentile(df['means'], np.arange(10, 105, 5)), np.inf],
        )
        disp_grouped = df.groupby('mean_bin')['dispersions']
        disp_median_bin = disp_grouped.median()
        # the next line raises the warning: "Mean of empty slice"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            disp_mad_bin = disp_grouped.apply(robust.mad)
            df['dispersions_norm'] = (
                df['dispersions'].values - disp_median_bin[df['mean_bin'].values].values
            ) / disp_mad_bin[df['mean_bin'].values].values
    else:
        raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
    dispersion_norm = df['dispersions_norm'].values
    if n_top_genes is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[
            ::-1
        ].sort()  # interestingly, np.argpartition is slightly slower
        if n_top_genes > adata.n_vars:
            logg.info('`n_top_genes` > `adata.n_var`, returning all genes.')
            n_top_genes = adata.n_vars
        disp_cut_off = dispersion_norm[n_top_genes - 1]
        gene_subset = np.nan_to_num(df['dispersions_norm'].values) >= disp_cut_off
        logg.debug(
            f'the {n_top_genes} top genes correspond to a '
            f'normalized dispersion cutoff of {disp_cut_off}'
        )
    else:
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce(
            (
                mean > min_mean,
                mean < max_mean,
                dispersion_norm > min_disp,
                dispersion_norm < max_disp,
            )
        )

    df['highly_variable'] = gene_subset
    return df


def highly_variable_genes(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: Optional[int] = None,
    min_disp: Optional[float] = 0.5,
    max_disp: Optional[float] = np.inf,
    min_mean: Optional[float] = 0.0125,
    max_mean: Optional[float] = 3,
    span: Optional[float] = 0.3,
    n_bins: int = 20,
    theta: float = 100,
    clip: Optional[float] = None,
    chunksize: int = 1000,
    flavor: Literal[
        'seurat', 'cell_ranger', 'seurat_v3', 'pearson_residuals'
    ] = 'seurat',
    subset: bool = False,
    inplace: bool = True,
    batch_key: Optional[str] = None,
    check_values: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    Annotate highly variable genes [Satija15]_ [Zheng17]_ [Stuart19]_.

    Expects logarithmized data, except when `flavor='seurat_v3'` or
    `flavor='pearson_residuals'`, in which count data is expected.

    Depending on `flavor`, this reproduces the R-implementations of Seurat
    [Satija15]_, Cell Ranger [Zheng17]_, and Seurat v3 [Stuart19]_, or uses
    analytical Peason residuals [Lause20]_.

    For the dispersion-based methods ([Satija15]_ and [Zheng17]_), the normalized
    dispersion is obtained by scaling with the mean and standard deviation of
    the dispersions for genes falling into a given bin for mean expression of
    genes. This means that for each bin of mean expression, highly variable
    genes are selected.

    For [Stuart19]_, a normalized variance for each gene is computed. First, the data
    are standardized (i.e., z-score normalization per feature) with a regularized
    standard deviation. Next, the normalized variance is computed as the variance
    of each gene after the transformation. Genes are ranked by the normalized variance.

    For [Lause20]_, Pearson residuals of a negative binomial offset model (with
    overdispersion theta shared across genes) are computed. By default, overdispersion
    theta=100 is used and residuals are clipped to sqrt(n). Finally, genes are ranked
    by residual variance.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    layer
        If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.
    n_top_genes
        Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'` or
        `flavor='pearson_residuals'`.
    min_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'` or
        `flavor='pearson_residuals'`.
    max_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'` or
        `flavor='pearson_residuals'`.
    min_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'` or
        `flavor='pearson_residuals'`.
    max_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'` or
        `flavor='pearson_residuals'`.
    span
        The fraction of the data (cells) used when estimating the variance in the loess
        model fit if `flavor='seurat_v3'`.
    n_bins
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`. Ignored if
        `flavor='pearson_residuals'`.
    theta
        If `flavor='pearson_residuals'`, this is the NB overdispersion parameter theta.
        Higher values correspond to less overdispersion (var = mean + mean^2/theta), and
        `theta=np.Inf` corresponds to a Poisson model.
    clip
        If `flavor='pearson_residuals'`, this determines if and how residuals are clipped:

        * If `None`, residuals are clipped to the interval [-sqrt(n), sqrt(n)],
        where n is the number of cells in the dataset (default behavior).
        * If any scalar c, residuals are clipped to the interval [-c, c]. Set
        `clip=np.Inf` for no clipping.

    chunksize
        If `flavor='pearson_residuals'`, this dertermines how many genes are processed at
        once while computing the residual variance. Choosing a smaller value will reduce
        the required memory.
    flavor
        Choose the flavor for identifying highly variable genes. For the dispersion
        based methods in their default workflows, Seurat passes the cutoffs whereas
        Cell Ranger passes `n_top_genes`.
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    batch_key
        If specified, highly-variable genes are selected within each batch separately and merged.
        This simple process avoids the selection of batch-specific genes and acts as a
        lightweight batch correction method. For all flavors, genes are first sorted
        by how many batches they are a HVG. For dispersion-based flavors ties are broken
        by normalized dispersion. If `flavor = 'seurat_v3'`, ties are broken by the median
        (across batches) rank based on within-batch normalized variance. If
        `flavor='pearson_residuals'`, ties are broken by the median rank (across batches)
        based on within-batch residual variance.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to True.
        Only used if `flavor='seurat_v3'` or `flavor='pearson_residuals'`.


    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pandas.DataFrame`) or
    updates `.var` with the following fields

    highly_variable : bool
        boolean indicator of highly-variable genes
    **means**
        means per gene
    **dispersions**
        For dispersion-based flavors, dispersions per gene
    **dispersions_norm**
        For dispersion-based flavors, normalized dispersions per gene
    **variances**
        For `flavor='seurat_v3'` and `flavor='pearson_residuals'`, variance per gene
    **variances_norm**
        For `flavor='seurat_v3'`, normalized variance per gene, averaged in
        the case of multiple batches
    **residual_variances**
        For `flavor='pearson_residuals'`, residual variance per gene. Averaged in the case of
        multiple batches.
    highly_variable_rank : float
        For `flavor='seurat_v3'`, rank of the gene according to normalized
        variance, median rank in the case of multiple batches
        For `flavor='pearson_residuals'`, rank of the gene according to residual
        variance, median rank in the case of multiple batches
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as HVG
    highly_variable_intersection : bool
        If batch_key is given, this denotes the genes that are highly variable in all batches

    Notes
    -----
    This function replaces :func:`~scanpy.pp.filter_genes_dispersion`.
    """

    if n_top_genes is not None and not all(
        m is None for m in [min_disp, max_disp, min_mean, max_mean]
    ):
        logg.info('If you pass `n_top_genes`, all cutoffs are ignored.')

    start = logg.info('extracting highly variable genes')

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`pp.highly_variable_genes` expects an `AnnData` argument, '
            'pass `inplace=False` if you want to return a `pd.DataFrame`.'
        )

    if flavor == 'seurat_v3':
        return _highly_variable_genes_seurat_v3(
            adata,
            layer=layer,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            check_values=check_values,
            span=span,
            subset=subset,
            inplace=inplace,
        )
    if flavor == 'pearson_residuals':
        if n_top_genes is None:
            raise ValueError(
                "`pp.highly_variable_genes` requires the argument `n_top_genes`"
                " for `flavor='pearson_residuals'`"
            )
        return _highly_variable_pearson_residuals(
            adata,
            layer=layer,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            theta=theta,
            clip=clip,
            chunksize=chunksize,
            subset=subset,
            check_values=check_values,
            inplace=inplace,
        )

    if batch_key is None:
        df = _highly_variable_genes_single_batch(
            adata,
            layer=layer,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
            n_top_genes=n_top_genes,
            n_bins=n_bins,
            flavor=flavor,
        )
    else:
        sanitize_anndata(adata)
        batches = adata.obs[batch_key].cat.categories
        df = []
        gene_list = adata.var_names
        for batch in batches:
            adata_subset = adata[adata.obs[batch_key] == batch]

            # Filter to genes that are in the dataset
            with settings.verbosity.override(Verbosity.error):
                filt = filter_genes(adata_subset, min_cells=1, inplace=False)[0]

            adata_subset = adata_subset[:, filt]

            hvg = _highly_variable_genes_single_batch(
                adata_subset,
                min_disp=min_disp,
                max_disp=max_disp,
                min_mean=min_mean,
                max_mean=max_mean,
                n_top_genes=n_top_genes,
                n_bins=n_bins,
                flavor=flavor,
            )

            # Add 0 values for genes that were filtered out
            missing_hvg = pd.DataFrame(
                np.zeros((np.sum(~filt), len(hvg.columns))),
                columns=hvg.columns,
            )
            missing_hvg['highly_variable'] = missing_hvg['highly_variable'].astype(bool)
            missing_hvg['gene'] = gene_list[~filt]
            hvg['gene'] = adata_subset.var_names.values
            hvg = hvg.append(missing_hvg, ignore_index=True)

            # Order as before filtering
            idxs = np.concatenate((np.where(filt)[0], np.where(~filt)[0]))
            hvg = hvg.loc[np.argsort(idxs)]

            df.append(hvg)

        df = pd.concat(df, axis=0)
        df['highly_variable'] = df['highly_variable'].astype(int)
        df = df.groupby('gene').agg(
            dict(
                means=np.nanmean,
                dispersions=np.nanmean,
                dispersions_norm=np.nanmean,
                highly_variable=np.nansum,
            )
        )
        df.rename(
            columns=dict(highly_variable='highly_variable_nbatches'), inplace=True
        )
        df['highly_variable_intersection'] = df['highly_variable_nbatches'] == len(
            batches
        )

        if n_top_genes is not None:
            # sort genes by how often they selected as hvg within each batch and
            # break ties with normalized dispersion across batches
            df.sort_values(
                ['highly_variable_nbatches', 'dispersions_norm'],
                ascending=False,
                na_position='last',
                inplace=True,
            )
            df['highly_variable'] = False
            df.highly_variable.iloc[:n_top_genes] = True
            df = df.loc[adata.var_names]
        else:
            df = df.loc[adata.var_names]
            dispersion_norm = df.dispersions_norm.values
            dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
            gene_subset = np.logical_and.reduce(
                (
                    df.means > min_mean,
                    df.means < max_mean,
                    df.dispersions_norm > min_disp,
                    df.dispersions_norm < max_disp,
                )
            )
            df['highly_variable'] = gene_subset

    logg.info('    finished', time=start)

    if inplace or subset:
        adata.uns['hvg'] = {'flavor': flavor}
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'dispersions\', float vector (adata.var)\n'
            '    \'dispersions_norm\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = df['highly_variable'].values
        adata.var['means'] = df['means'].values
        adata.var['dispersions'] = df['dispersions'].values
        adata.var['dispersions_norm'] = df['dispersions_norm'].values.astype(
            'float32', copy=False
        )
        if batch_key is not None:
            adata.var['highly_variable_nbatches'] = df[
                'highly_variable_nbatches'
            ].values
            adata.var['highly_variable_intersection'] = df[
                'highly_variable_intersection'
            ].values
        if subset:
            adata._inplace_subset_var(df['highly_variable'].values)
    else:
        return df

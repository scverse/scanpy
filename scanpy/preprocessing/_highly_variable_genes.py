import warnings
from typing import Optional

import numpy as np
import pandas as pd
from numba import jit
import scipy.sparse as sp_sparse
from anndata import AnnData


from .. import logging as logg
from .._settings import settings, Verbosity
from .._utils import sanitize_anndata
from .._compat import Literal
from ._utils import _get_mean_var
from ._distributed import materialize_as_ndarray
from ._simple import filter_genes


def _highly_variable_genes_single_batch(
    adata: AnnData,
    min_disp: Optional[float] = None,
    max_disp: Optional[float] = None,
    min_mean: Optional[float] = None,
    max_mean: Optional[float] = None,
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

    if n_top_genes is not None and not all(
        m is None for m in [min_disp, max_disp, min_mean, max_mean]
    ):
        logg.info('If you pass `n_top_genes`, all cutoffs are ignored.')

    if min_disp is None:
        min_disp = 0.5
    if min_mean is None:
        min_mean = 0.0125
    if max_mean is None:
        max_mean = 3
    if max_disp is None:
        max_disp = np.inf

    X = adata.X
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
    dispersion_norm = df['dispersions_norm'].values.astype('float32')
    if n_top_genes is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[
            ::-1
        ].sort()  # interestingly, np.argpartition is slightly slower
        if n_top_genes > adata.n_vars:
            logg.info(f'`n_top_genes` > `adata.n_var`, returning all genes.')
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
    min_disp: Optional[float] = None,
    max_disp: Optional[float] = None,
    min_mean: Optional[float] = None,
    max_mean: Optional[float] = None,
    n_top_genes: Optional[int] = None,
    n_bins: int = 20,
    flavor: Literal['seurat', 'cell_ranger'] = 'seurat',
    subset: bool = False,
    inplace: bool = True,
    batch_key: Optional[str] = None,
) -> Optional[np.recarray]:
    """\
    Annotate highly variable genes [Satija15]_ [Zheng17]_.

    Expects logarithmized data.

    Depending on `flavor`, this reproduces the R-implementations of Seurat
    [Satija15]_ and Cell Ranger [Zheng17]_.

    The normalized dispersion is obtained by scaling with the mean and standard
    deviation of the dispersions for genes falling into a given bin for mean
    expression of genes. This means that for each bin of mean expression, highly
    variable genes are selected.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    min_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    n_top_genes
        Number of highly-variable genes to keep.
    n_bins
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    flavor
        Choose the flavor for computing normalized dispersion. In their default
        workflows, Seurat passes the cutoffs whereas Cell Ranger passes
        `n_top_genes`.
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    batch_key
        If specified, highly-variable genes are selected within each batch separately and merged.
        This simple process avoids the selection of batch-specific genes and acts as a
        lightweight batch correction method.

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~numpy.recarray`) or
    updates `.var` with the following fields

    highly_variable : bool
        boolean indicator of highly-variable genes
    **means**
        means per gene
    **dispersions**
        dispersions per gene
    **dispersions_norm**
        normalized dispersions per gene
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as HVG
    highly_variable_intersection : bool
        If batch_key is given, this denotes the genes that are highly variable in all batches

    Notes
    -----
    This function replaces :func:`~scanpy.pp.filter_genes_dispersion`.
    """

    start = logg.info('extracting highly variable genes')

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`pp.highly_variable_genes` expects an `AnnData` argument, '
            'pass `inplace=False` if you want to return a `np.recarray`.'
        )

    if batch_key is None:
        df = _highly_variable_genes_single_batch(
            adata,
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
                np.zeros((np.sum(~filt), len(hvg.columns))), columns=hvg.columns,
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
            df.loc[:n_top_genes, 'highly_variable'] = True
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
        arrays = [
            df['highly_variable'].values,
            df['means'].values,
            df['dispersions'].values,
            df['dispersions_norm'].values.astype('float32', copy=False),
        ]
        dtypes = [
            ('highly_variable', np.bool_),
            ('means', 'float32'),
            ('dispersions', 'float32'),
            ('dispersions_norm', 'float32'),
        ]
        if batch_key is not None:
            arrays.extend(
                [
                    df['highly_variable_nbatches'].values,
                    df['highly_variable_intersection'].values,
                ]
            )
            dtypes.append(
                [
                    ('highly_variable_nbatches', int),
                    ('highly_variable_intersection', np.bool_),
                ]
            )
        return np.rec.fromarrays(arrays, dtype=dtypes)


def highly_variable_genes_seurat_v3(
    adata: AnnData,
    n_top_genes: int = 2000,
    batch_key: Optional[str] = None,
    lowess_frac: Optional[float] = 0.15,
    subset: bool = False,
    inplace: bool = True,
    use_lowess: bool = False,
):
    """\
    Annotate highly variable genes [Stuart19]_.

    Expects raw count data.

    The major difference in this implementation is the use of lowess insted of loess.

    For further details of the sparse arithmetic see https://www.overleaf.com/read/ckptrbgzzzpg

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    n_top_genes
        Number of highly-variable genes to keep.
    batch_key
        If specified, highly-variable genes are selected within each batch separately and merged.
        This simple process avoids the selection of batch-specific genes and acts as a
        lightweight batch correction method.
    lowess_frac
        The fraction of the data (cells) used when estimating the variance in the lowess model fit.
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    use_lowess
        Whether to use statsmodels lowess implementation

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~numpy.recarray`) or
    updates `.var` with the following fields

    highly_variable : bool
        boolean indicator of highly-variable genes
    **variances_norm**
        normalized variance per gene, averaged in the case of multiple batches
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as HVG
    highly_variable_intersection : bool
        If batch_key is given, this denotes the genes that are highly variable in all batches
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    if batch_key is None:
        batch_info = pd.Categorical(np.zeros(adata.shape[0], dtype=int))
    else:
        batch_info = adata.obs[batch_key]

    norm_gene_vars = []
    for b in np.unique(batch_info):

        mean, var = _get_mean_var(adata[batch_info == b].X)
        not_const = var > 0
        estimat_var = np.zeros(adata.shape[1])

        y = np.log10(var[not_const])
        x = np.log10(mean[not_const])
        if use_lowess is True:
            estimat_var[not_const] = lowess(y, x, frac=lowess_frac, return_sorted=False)
        else:
            estimat_var[not_const] = lowess_quadratic(y, x)
        reg_std = np.sqrt(10 ** estimat_var)

        batch_counts = adata[batch_info == b].X.astype(np.float64).copy()
        # clip large values as in Seurat
        N = np.sum(batch_info == b)
        vmax = np.sqrt(N)
        clip_val = reg_std * vmax + mean
        if sp_sparse.issparse(batch_counts):
            batch_counts = sp_sparse.csr_matrix(batch_counts)
            mask = batch_counts.data > clip_val[batch_counts.indices]
            batch_counts.data[mask] = clip_val[batch_counts.indices[mask]]
        else:
            clip_val_broad = np.broadcast_to(clip_val, batch_counts.shape)
            np.putmask(
                batch_counts, batch_counts > clip_val_broad, clip_val_broad,
            )

        if sp_sparse.issparse(batch_counts):
            squared_batch_counts_sum = np.array(batch_counts.power(2).sum(axis=0))
            batch_counts_sum = np.array(batch_counts.sum(axis=0))
        else:
            squared_batch_counts_sum = np.square(batch_counts).sum(axis=0)
            batch_counts_sum = batch_counts.sum(axis=0)

        norm_gene_var = (1 / ((N - 1) * np.square(reg_std))) * (
            (N * np.square(mean))
            + squared_batch_counts_sum
            - 2 * batch_counts_sum * mean
        )
        norm_gene_vars.append(norm_gene_var.reshape(1, -1))

    norm_gene_vars = np.concatenate(norm_gene_vars, axis=0)
    # argsort twice gives ranks
    ranked_norm_gene_vars = np.argsort(np.argsort(norm_gene_vars, axis=1), axis=1)
    median_ranked = np.median(ranked_norm_gene_vars, axis=0)

    num_batches_high_var = np.sum(
        ranked_norm_gene_vars >= (adata.X.shape[1] - n_top_genes), axis=0
    )
    df = pd.DataFrame(index=np.array(adata.var_names))
    df["highly_variable_nbatches"] = num_batches_high_var
    df["highly_variable_median_rank"] = median_ranked
    df["variances_norm"] = np.mean(norm_gene_vars, axis=0)
    df["means"] = mean
    df["variances"] = var

    df.sort_values(
        ["highly_variable_nbatches", "highly_variable_median_rank"],
        ascending=False,
        na_position="last",
        inplace=True,
    )
    df["highly_variable"] = False
    df.loc[: int(n_top_genes), "highly_variable"] = True
    df = df.loc[adata.var_names]

    adata.var["highly_variable"] = df["highly_variable"].values
    if batch_key is not None:
        batches = adata.obs[batch_key].cat.categories
        adata.var["highly_variable_nbatches"] = df["highly_variable_nbatches"].values
        adata.var["highly_variable_intersection"] = df[
            "highly_variable_nbatches"
        ] == len(batches)
    adata.var["highly_variable_median_rank"] = df["highly_variable_median_rank"].values

    if inplace or subset:
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'variances\', float vector (adata.var)\n'
            '    \'variances_norm\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = df['highly_variable'].values
        adata.var['means'] = df['means'].values
        adata.var['variances'] = df['variances'].values
        adata.var['variances_norm'] = df['variances_norm'].values.astype(
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
        arrays = [
            df['highly_variable'].values,
            df['means'].values,
            df['variances'].values,
            df['variances_norm'].values.astype('float32', copy=False),
        ]
        dtypes = [
            ('highly_variable', np.bool_),
            ('means', 'float32'),
            ('variances', 'float32'),
            ('variances_norm', 'float32'),
        ]
        if batch_key is not None:
            arrays.extend(
                [
                    df['highly_variable_nbatches'].values,
                    df['highly_variable_intersection'].values,
                ]
            )
            dtypes.append(
                [
                    ('highly_variable_nbatches', int),
                    ('highly_variable_intersection', np.bool_),
                ]
            )
        return np.rec.fromarrays(arrays, dtype=dtypes)


@jit(nopython=True)
def lowess_quadratic(y, x, alpha=0.3):
    """\
    Lowess smoother: Locally weighted regression.

    This is an adaptation of
     https://xavierbourretsicotte.github.io/loess.html#Implementation-in-python-(using-bell-shaped-kernel)

    with additional use of quadratic polynomial.

    Parameters
    ----------
    alpha
        Smoothing span. A larger value for alpha will result in a
        smoother curve.
    iter
        The number of robustifying iterations

    Returns
    -------
    Estimated values from lowess fit
    """
    n = len(x)
    r = int(np.ceil(alpha * n))
    h = 1 / np.array([np.sort(np.abs(x - x[i]))[r] for i in range(n)]).reshape(1, -1)
    w = np.abs((x.reshape(-1, 1) - x.reshape(1, -1))) * h
    out_w = np.empty_like(w)
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            if w[i, j] < 0:
                out_w[i, j] = 0
            elif w[i, j] > 1:
                out_w[i, j] = 1
            else:
                out_w[i, j] = w[i, j]
    w = out_w
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    x_sq = np.square(x)
    for i in range(n):
        weights = w[:, i]
        b = np.array(
            [np.sum(weights * y), np.sum(weights * y * x), np.sum(weights * y * x_sq),]
        )
        A = np.array(
            [
                [np.sum(weights), np.sum(weights * x), np.sum(weights * x_sq)],
                [
                    np.sum(weights * x),
                    np.sum(weights * x * x),
                    np.sum(weights * x_sq * x),
                ],
                [
                    np.sum(weights * x_sq),
                    np.sum(weights * x * x_sq),
                    np.sum(weights * x_sq * x_sq),
                ],
            ]
        )
        beta = np.linalg.solve(A, b)
        yest[i] = beta[0] + beta[1] * x[i] + beta[2] * x_sq[i]

    return yest

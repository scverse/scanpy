import warnings
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .. import logging as logg
from ._distributed import materialize_as_ndarray
from ._utils import _get_mean_var
from ..utils import sanitize_anndata


def _highly_variable_genes_single_batch(
    adata,
    min_disp=None, max_disp=None,
    min_mean=None, max_mean=None,
    n_top_genes=None,
    n_bins=20,
    flavor='seurat',
) -> pd.DataFrame:
    """Internal function for annotating highly variable genes [Satija15]_ [Zheng17]_.

    Expects logarithmized data.

    Depending on `flavor`, this reproduces the R-implementations of Seurat
    [Satija15]_ and Cell Ranger [Zheng17]_.

    The normalized dispersion is obtained by scaling with the mean and standard
    deviation of the dispersions for genes falling into a given bin for mean
    expression of genes. This means that for each bin of mean expression, highly
    variable genes are selected.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_mean : `float`, optional (default: 0.0125)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_mean : `float`, optional (default: 3)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    min_disp : `float`, optional (default: 0.5)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_disp : `float`, optional (default: `None`)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    n_top_genes : `int` or `None`, optional (default: `None`)
        Number of highly-variable genes to keep.
    n_bins : `int`, optional (default: 20)
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    flavor : `{'seurat', 'cell_ranger'}`, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion. In their default
        workflows, Seurat passes the cutoffs whereas Cell Ranger passes
        `n_top_genes`.

    Returns
    -------
    highly_variable_data_frame
        A DataFrame that contains colums highly_variable, means, dispersions and dispersions_norm.
    """

    if n_top_genes is not None and not all([
            min_disp is None, max_disp is None, min_mean is None, max_mean is None]):
        logg.info('If you pass `n_top_genes`, all cutoffs are ignored.')

    if min_disp is None: min_disp = 0.5
    if min_mean is None: min_mean = 0.0125
    if max_mean is None: max_mean = 3
    if max_disp is None: max_disp = np.inf

    X = np.expm1(adata.X) if flavor == 'seurat' else adata.X
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
            logg.msg(
                'Gene indices {} fell into a single bin: their '
                'normalized dispersion was set to 1.\n    '
                'Decreasing `n_bins` will likely avoid this effect.'
                .format(gen_indices),
                v=4
            )
        # Circumvent pandas 0.23 bug. Both sides of the assignment have dtype==float32,
        # but there’s still a dtype error without “.value”.
        disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[one_gene_per_bin.values].values
        disp_mean_bin[one_gene_per_bin.values] = 0
        # actually do the normalization
        df['dispersions_norm'] = (
            (
                df['dispersions'].values  # use values here as index differs
                - disp_mean_bin[df['mean_bin'].values].values
            ) / disp_std_bin[df['mean_bin'].values].values
        )
    elif flavor == 'cell_ranger':
        from statsmodels import robust
        df['mean_bin'] = pd.cut(df['means'], np.r_[
            -np.inf,
            np.percentile(df['means'], np.arange(10, 105, 5)),
            np.inf
        ])
        disp_grouped = df.groupby('mean_bin')['dispersions']
        disp_median_bin = disp_grouped.median()
        # the next line raises the warning: "Mean of empty slice"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            disp_mad_bin = disp_grouped.apply(robust.mad)
        df['dispersions_norm'] = (
            np.abs(
                df['dispersions'].values
                - disp_median_bin[df['mean_bin'].values].values
            ) / disp_mad_bin[df['mean_bin'].values].values
        )
    else:
        raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
    dispersion_norm = df['dispersions_norm'].values.astype('float32')
    if n_top_genes is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_genes-1]
        gene_subset = np.nan_to_num(df['dispersions_norm'].values) >= disp_cut_off
        logg.msg(
            'the {} top genes correspond to a normalized dispersion cutoff of {}'
            .format(n_top_genes, disp_cut_off),
            v=5,
        )
    else:
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce((
            mean > min_mean, mean < max_mean,
            dispersion_norm > min_disp,
            dispersion_norm < max_disp,
        ))

    df['highly_variable'] = gene_subset
    return df

def highly_variable_genes(
    adata,
    min_disp=None, max_disp=None,
    min_mean=None, max_mean=None,
    n_top_genes=None,
    n_bins=20,
    flavor='seurat',
    subset=False,
    inplace=True,
    batch_key=None,
) -> Optional[np.recarray]:
    """Annotate highly variable genes [Satija15]_ [Zheng17]_.

    Expects logarithmized data.

    Depending on `flavor`, this reproduces the R-implementations of Seurat
    [Satija15]_ and Cell Ranger [Zheng17]_.

    The normalized dispersion is obtained by scaling with the mean and standard
    deviation of the dispersions for genes falling into a given bin for mean
    expression of genes. This means that for each bin of mean expression, highly
    variable genes are selected.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_mean : `float`, optional (default: 0.0125)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_mean : `float`, optional (default: 3)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    min_disp : `float`, optional (default: 0.5)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    max_disp : `float`, optional (default: `None`)
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored.
    n_top_genes : `int` or `None`, optional (default: `None`)
        Number of highly-variable genes to keep.
    n_bins : `int`, optional (default: 20)
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    flavor : `{'seurat', 'cell_ranger'}`, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion. In their default
        workflows, Seurat passes the cutoffs whereas Cell Ranger passes
        `n_top_genes`.
    subset : `bool`, optional (default: `False`)
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace : `bool`, optional (default: `True`)
        Whether to place calculated metrics in `.var` or return them.
    batch_key : `str`, optional (default: `None`)
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

    logg.msg('extracting highly variable genes', r=True, v=4)

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`pp.highly_variable_genes` expects an `AnnData` argument, '
            'pass `inplace=False` if you want to return a `np.recarray`.')

    if batch_key is None:
        df = _highly_variable_genes_single_batch(adata,
                                                 min_disp=min_disp, max_disp=max_disp,
                                                 min_mean=min_mean, max_mean=max_mean,
                                                 n_top_genes=n_top_genes,
                                                 n_bins=n_bins,
                                                 flavor=flavor)
    else:
        sanitize_anndata(adata)
        batches = adata.obs[batch_key].cat.categories
        df = []
        for batch in batches:
            adata_subset = adata[adata.obs[batch_key] == batch]
            hvg = _highly_variable_genes_single_batch(adata_subset,
                                                      min_disp=min_disp, max_disp=max_disp,
                                                      min_mean=min_mean, max_mean=max_mean,
                                                      n_top_genes=n_top_genes,
                                                      n_bins=n_bins,
                                                      flavor=flavor)
            hvg['gene'] = adata.var_names.values
            df.append(hvg)

        df = pd.concat(df, axis=0)
        df['highly_variable'] = df['highly_variable'].astype(int)
        df = df.groupby('gene').agg({'means': np.nanmean,
                                     'dispersions': np.nanmean,
                                     'dispersions_norm': np.nanmean,
                                     'highly_variable': np.nansum})
        df.rename(columns={'highly_variable': 'highly_variable_nbatches'}, inplace=True)
        df['highly_variable_intersection'] = df['highly_variable_nbatches'] == len(batches)

        if n_top_genes is not None:
            # sort genes by how often they selected as hvg within each batch and
            # break ties with normalized dispersion across batches
            df.sort_values(['highly_variable_nbatches', 'dispersions_norm'],
                           ascending=False, na_position='last', inplace=True)
            df['highly_variable'] = False
            df.loc[:n_top_genes, 'highly_variable'] = True
            df = df.loc[adata.var_names]
        else:
            df = df.loc[adata.var_names]
            dispersion_norm = df.dispersions_norm.values
            dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
            gene_subset = np.logical_and.reduce((
                df.means > min_mean, df.means < max_mean,
                df.dispersions_norm > min_disp,
                df.dispersions_norm < max_disp,
            ))
            df['highly_variable'] = gene_subset

    logg.msg('    finished', time=True, v=4)

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
        adata.var['dispersions_norm'] = df['dispersions_norm'].values.astype('float32', copy=False)
        if batch_key is not None:
            adata.var['highly_variable_nbatches'] = df['highly_variable_nbatches'].values
            adata.var['highly_variable_intersection'] = df['highly_variable_intersection'].values
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
            arrays.extend([df['highly_variable_nbatches'].values,
                           df['highly_variable_intersection'].values])
            dtypes.append([('highly_variable_nbatches', int),
                           ('highly_variable_intersection', np.bool_)])
        return np.rec.fromarrays(arrays, dtype=dtypes)

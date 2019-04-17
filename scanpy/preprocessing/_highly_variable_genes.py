import warnings
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .. import logging as logg
from ._distributed import materialize_as_ndarray
from ._utils import _get_mean_var


def highly_variable_genes(
    adata,
    min_disp=None, max_disp=None,
    min_mean=None, max_mean=None,
    n_top_genes=None,
    n_bins=20,
    flavor='seurat',
    subset=False,
    inplace=True
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

    Notes
    -----
    This function replaces :func:`~scanpy.pp.filter_genes_dispersion`.
    """
    logg.msg('extracting highly variable genes', r=True, v=4)

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`pp.highly_variable_genes` expects an `AnnData` argument, '
            'pass `inplace=False` if you want to return a `np.recarray`.')

    if n_top_genes is not None and not all([
            min_disp is None, max_disp is None, min_mean is None, max_mean is None]):
        logg.info('If you pass `n_top_genes`, all cutoffs are ignored.')
    if min_disp is None: min_disp = 0.5
    if min_mean is None: min_mean = 0.0125
    if max_mean is None: max_mean = 3

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
    df['mean'] = mean
    df['dispersion'] = dispersion
    if flavor == 'seurat':
        df['mean_bin'] = pd.cut(df['mean'], bins=n_bins)
        disp_grouped = df.groupby('mean_bin')['dispersion']
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
        df['dispersion_norm'] = (
            (
                df['dispersion'].values  # use values here as index differs
                - disp_mean_bin[df['mean_bin'].values].values
            ) / disp_std_bin[df['mean_bin'].values].values
        )
    elif flavor == 'cell_ranger':
        from statsmodels import robust
        df['mean_bin'] = pd.cut(df['mean'], np.r_[
            -np.inf,
            np.percentile(df['mean'], np.arange(10, 105, 5)),
            np.inf
        ])
        disp_grouped = df.groupby('mean_bin')['dispersion']
        disp_median_bin = disp_grouped.median()
        # the next line raises the warning: "Mean of empty slice"
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            disp_mad_bin = disp_grouped.apply(robust.mad)
        df['dispersion_norm'] = (
            np.abs(
                df['dispersion'].values
                - disp_median_bin[df['mean_bin'].values].values
            ) / disp_mad_bin[df['mean_bin'].values].values
        )
    else:
        raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
    dispersion_norm = df['dispersion_norm'].values.astype('float32')
    if n_top_genes is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_genes-1]
        gene_subset = np.nan_to_num(df['dispersion_norm'].values) >= disp_cut_off
        logg.msg(
            'the {} top genes correspond to a normalized dispersion cutoff of'
            .format(n_top_genes, disp_cut_off),
            v=5,
        )
    else:
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce((
            mean > min_mean, mean < max_mean,
            dispersion_norm > min_disp,
            dispersion_norm < max_disp,
        ))

    logg.msg('    finished', time=True, v=4)

    if inplace or subset:
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'dispersions\', float vector (adata.var)\n'
            '    \'dispersions_norm\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = gene_subset
        adata.var['means'] = df['mean'].values
        adata.var['dispersions'] = df['dispersion'].values
        adata.var['dispersions_norm'] = df['dispersion_norm'].values.astype('float32', copy=False)
        if subset:
            adata._inplace_subset_var(gene_subset)
    else:
        arrays = (
             gene_subset,
             df['mean'].values,
             df['dispersion'].values,
             df['dispersion_norm'].values.astype('float32', copy=False)
        )
        dtypes = [
            ('highly_variable', np.bool_),
            ('means', 'float32'),
            ('dispersions', 'float32'),
            ('dispersions_norm', 'float32'),
        ]
        return np.rec.fromarrays(arrays, dtype=dtypes)

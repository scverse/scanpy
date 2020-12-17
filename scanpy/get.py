"""This module contains helper functions for accessing data."""
from typing import Optional, Iterable, Tuple, Union

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix, issparse

from anndata import AnnData

# --------------------------------------------------------------------------------
# Plotting data helpers
# --------------------------------------------------------------------------------


# TODO: implement diffxpy method, make singledispatch
def rank_genes_groups_df(
    adata: AnnData,
    group: Union[str, Iterable[str]],
    *,
    key: str = "rank_genes_groups",
    pval_cutoff: Optional[float] = None,
    log2fc_min: Optional[float] = None,
    log2fc_max: Optional[float] = None,
    gene_symbols: Optional[str] = None,
) -> pd.DataFrame:
    """\
    :func:`scanpy.tl.rank_genes_groups` results in the form of a
    :class:`~pandas.DataFrame`.

    Params
    ------
    adata
        Object to get results from.
    group
        Which group (as in :func:`scanpy.tl.rank_genes_groups`'s `groupby`
        argument) to return results from. Can be a list. All groups are
        returned if groups is `None`.
    key
        Key differential expression groups were stored under.
    pval_cutoff
        Minimum adjusted pval to return.
    log2fc_min
        Minumum logfc to return.
    log2fc_max
        Maximum logfc to return.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. Specifying
        this will add that column to the returned dataframe.

    Example
    -------
    >>> import scanpy as sc
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(pbmc, groupby="louvain", use_raw=True)
    >>> dedf = sc.get.rank_genes_groups_df(pbmc, group="0")
    """
    if isinstance(group, str):
        group = [group]
    if group is None:
        group = list(adata.uns[key]['names'].dtype.names)
    colnames = ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']

    d = [pd.DataFrame(adata.uns[key][c])[group] for c in colnames]
    d = pd.concat(d, axis=1, names=[None, 'group'], keys=colnames)
    d = d.stack(level=1).reset_index()
    d['group'] = pd.Categorical(d['group'], categories=group)
    d = d.sort_values(['group', 'level_0']).drop(columns='level_0')

    if pval_cutoff is not None:
        d = d[d["pvals_adj"] < pval_cutoff]
    if log2fc_min is not None:
        d = d[d["logfoldchanges"] > log2fc_min]
    if log2fc_max is not None:
        d = d[d["logfoldchanges"] < log2fc_max]
    if gene_symbols is not None:
        d = d.join(adata.var[gene_symbols], on="names")

    for pts, name in {'pts': 'pct_nz_group', 'pts_rest': 'pct_nz_reference'}.items():
        if pts in adata.uns[key]:
            pts_df = (
                adata.uns[key][pts][group]
                .rename_axis(index='names')
                .reset_index()
                .melt(id_vars='names', var_name='group', value_name=name)
            )
            d = d.merge(pts_df)

    # remove group column for backward compat if len(group) == 1
    if len(group) == 1:
        d.drop(columns='group', inplace=True)

    return d.reset_index(drop=True)


def obs_df(
    adata: AnnData,
    keys: Iterable[str] = (),
    obsm_keys: Iterable[Tuple[str, int]] = (),
    *,
    layer: str = None,
    gene_symbols: str = None,
    use_raw: bool = False,
) -> pd.DataFrame:
    """\
    Return values for observations in adata.

    Params
    ------
    adata
        AnnData object to get values from.
    keys
        Keys from either `.var_names`, `.var[gene_symbols]`, or `.obs.columns`.
    obsm_keys
        Tuple of `(key from obsm, column index of obsm[key])`.
    layer
        Layer of `adata` to use as expression values.
    gene_symbols
        Column of `adata.var` to search for `keys` in.
    use_raw
        Whether to get expression values from `adata.raw`.

    Returns
    -------
    A dataframe with `adata.obs_names` as index, and values specified by `keys`
    and `obsm_keys`.

    Examples
    --------
    Getting value for plotting:

    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> plotdf = sc.get.obs_df(
            pbmc,
            keys=["CD8B", "n_genes"],
            obsm_keys=[("X_umap", 0), ("X_umap", 1)]
        )
    >>> plotdf.plot.scatter("X_umap0", "X_umap1", c="CD8B")

    Calculating mean expression for marker genes by cluster:

    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> marker_genes = ['CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ']
    >>> genedf = sc.get.obs_df(
            pbmc,
            keys=["louvain", *marker_genes]
        )
    >>> grouped = genedf.groupby("louvain")
    >>> mean, var = grouped.mean(), grouped.var()
    """
    if use_raw:
        assert (
            layer is None
        ), "Cannot specify use_raw=True and a layer at the same time."
        if gene_symbols is not None:
            gene_names = pd.Series(
                adata.raw.var_names, index=adata.raw.var[gene_symbols]
            )
        else:
            gene_names = pd.Series(adata.raw.var_names, index=adata.raw.var_names)
    else:
        if gene_symbols is not None:
            gene_names = pd.Series(adata.var_names, index=adata.var[gene_symbols])
        else:
            gene_names = pd.Series(adata.var_names, index=adata.var_names)
    obs_names = []
    var_names = []
    var_symbol = []
    not_found = []
    for key in keys:
        if key in adata.obs.columns:
            obs_names.append(key)
        elif key in gene_names.index:
            var_names.append(gene_names[key])
            var_symbol.append(key)
        else:
            not_found.append(key)
    if len(not_found) > 0:
        if use_raw:
            if gene_symbols is None:
                gene_error = "`adata.raw.var_names`"
            else:
                gene_error = "gene_symbols column `adata.raw.var[{}].values`".format(
                    gene_symbols
                )
        else:
            if gene_symbols is None:
                gene_error = "`adata.var_names`"
            else:
                gene_error = "gene_symbols column `adata.var[{}].values`".format(
                    gene_symbols
                )
        raise KeyError(
            f"Could not find keys '{not_found}' in columns of `adata.obs` or in"
            f" {gene_error}."
        )

    # Make df
    df = pd.DataFrame(index=adata.obs.index)

    # add var values
    if len(var_names) > 0:
        X = _get_obs_rep(adata, layer=layer, use_raw=use_raw)
        if use_raw:
            var_idx = adata.raw.var_names.get_indexer(var_names)
        else:
            var_idx = adata.var_names.get_indexer(var_names)

        # for backed AnnData is important that the indices are ordered
        if adata.isbacked:
            var_order = np.argsort(var_idx)
            matrix = X[:, var_idx[var_order]][:, np.argsort(var_order)]
        else:
            matrix = X[:, var_idx]

        from scipy.sparse import issparse

        if issparse(matrix):
            matrix = matrix.toarray()
        df = df.join(pd.DataFrame(matrix, columns=var_symbol, index=adata.obs.index))

    # add obs values
    if len(obs_names) > 0:
        df = df.join(adata.obs[obs_names])

    # reorder columns to given order
    df = df[keys]
    for k, idx in obsm_keys:
        added_k = f"{k}-{idx}"
        val = adata.obsm[k]
        if isinstance(val, np.ndarray):
            df[added_k] = np.ravel(val[:, idx])
        elif isinstance(val, spmatrix):
            df[added_k] = np.ravel(val[:, idx].toarray())
        elif isinstance(val, pd.DataFrame):
            df[added_k] = val.loc[:, idx]
    return df


def var_df(
    adata: AnnData,
    keys: Iterable[str] = (),
    varm_keys: Iterable[Tuple[str, int]] = (),
    *,
    layer: str = None,
) -> pd.DataFrame:
    """\
    Return values for observations in adata.

    Params
    ------
    adata
        AnnData object to get values from.
    keys
        Keys from either `.obs_names`, or `.var.columns`.
    varm_keys
        Tuple of `(key from varm, column index of varm[key])`.
    layer
        Layer of `adata` to use as expression values.

    Returns
    -------
    A dataframe with `adata.var_names` as index, and values specified by `keys`
    and `varm_keys`.
    """
    # Argument handling
    obs_names = []
    var_names = []
    not_found = []
    for key in keys:
        if key in adata.obs_names:
            obs_names.append(key)
        elif key in adata.var.columns:
            var_names.append(key)
        else:
            not_found.append(key)
    if len(not_found) > 0:
        raise KeyError(
            f"Could not find keys '{not_found}' in columns of `adata.var` or"
            " in `adata.obs_names`."
        )

    # initialize df
    df = pd.DataFrame(index=adata.var.index)

    # add obs values
    if len(obs_names) > 0:
        X = _get_obs_rep(adata, layer=layer)
        obs_idx = adata.obs_names.get_indexer(obs_names)

        # for backed AnnData is important that the indices are ordered
        if adata.isbacked:
            obs_order = np.argsort(obs_idx)
            matrix = X[obs_idx[obs_order], :][np.argsort(obs_order)]
        else:
            matrix = X[obs_idx, :]
        from scipy.sparse import issparse

        if issparse(matrix):
            matrix = matrix.toarray()

        df = df.join(pd.DataFrame(matrix.T, columns=obs_names, index=adata.var.index))

    # add obs values
    if len(var_names) > 0:
        df = df.join(adata.var[var_names])

    # reorder columns to given order
    df = df[keys]

    for k, idx in varm_keys:
        added_k = f"{k}-{idx}"
        val = adata.varm[k]
        if isinstance(val, np.ndarray):
            df[added_k] = np.ravel(val[:, idx])
        elif isinstance(val, spmatrix):
            df[added_k] = np.ravel(val[:, idx].toarray())
        elif isinstance(val, pd.DataFrame):
            df[added_k] = val.loc[:, idx]
    return df


def _get_obs_rep(adata, *, use_raw=False, layer=None, obsm=None, obsp=None):
    """
    Choose array aligned with obs annotation.
    """
    # https://github.com/theislab/scanpy/issues/1546
    if not isinstance(use_raw, bool):
        raise TypeError(f"use_raw expected to be bool, was {type(use_raw)}.")

    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made <= 1
    if choices_made == 0:
        return adata.X
    elif is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    elif is_obsm:
        return adata.obsm[obsm]
    elif is_obsp:
        return adata.obsp[obsp]
    else:
        assert False, (
            "That was unexpected. Please report this bug at:\n\n\t"
            " https://github.com/theislab/scanpy/issues"
        )


def _set_obs_rep(adata, val, *, use_raw=False, layer=None, obsm=None, obsp=None):
    """
    Set value for observation rep.
    """
    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made <= 1
    if choices_made == 0:
        adata.X = val
    elif is_layer:
        adata.layers[layer] = val
    elif use_raw:
        adata.raw.X = val
    elif is_obsm:
        adata.obsm[obsm] = val
    elif is_obsp:
        adata.obsp[obsp] = val
    else:
        assert False, (
            "That was unexpected. Please report this bug at:\n\n\t"
            " https://github.com/theislab/scanpy/issues"
        )

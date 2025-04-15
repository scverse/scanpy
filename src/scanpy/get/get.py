"""Helper functions for accessing data."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.typing import NDArray
from packaging.version import Version

from .._compat import CSBase

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable
    from typing import Any, Literal

    from anndata._core.sparse_dataset import BaseCompressedSparseDataset
    from anndata._core.views import ArrayView

    from .._compat import DaskArray


# --------------------------------------------------------------------------------
# Plotting data helpers
# --------------------------------------------------------------------------------


# TODO: implement diffxpy method, make singledispatch
def rank_genes_groups_df(  # noqa: PLR0912
    adata: AnnData,
    group: str | Iterable[str] | None,
    *,
    key: str = "rank_genes_groups",
    pval_cutoff: float | None = None,
    log2fc_min: float | None = None,
    log2fc_max: float | None = None,
    gene_symbols: str | None = None,
) -> pd.DataFrame:
    """Get :func:`scanpy.tl.rank_genes_groups` results in the form of a :class:`~pandas.DataFrame`.

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
        Return only adjusted p-values below the  cutoff.
    log2fc_min
        Minimum logfc to return.
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
        group = list(adata.uns[key]["names"].dtype.names)
    method = adata.uns[key]["params"]["method"]
    if method == "logreg":
        colnames = ["names", "scores"]
    else:
        colnames = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]

    d = [pd.DataFrame(adata.uns[key][c])[group] for c in colnames]
    d = pd.concat(d, axis=1, names=[None, "group"], keys=colnames)
    if Version(pd.__version__) >= Version("2.1"):
        d = d.stack(level=1, future_stack=True).reset_index()
    else:
        d = d.stack(level=1).reset_index()
    d["group"] = pd.Categorical(d["group"], categories=group)
    d = d.sort_values(["group", "level_0"]).drop(columns="level_0")

    if method != "logreg":
        if pval_cutoff is not None:
            d = d[d["pvals_adj"] < pval_cutoff]
        if log2fc_min is not None:
            d = d[d["logfoldchanges"] > log2fc_min]
        if log2fc_max is not None:
            d = d[d["logfoldchanges"] < log2fc_max]
    if gene_symbols is not None:
        d = d.join(adata.var[gene_symbols], on="names")

    for pts, name in {"pts": "pct_nz_group", "pts_rest": "pct_nz_reference"}.items():
        if pts in adata.uns[key]:
            pts_df = (
                adata.uns[key][pts][group]
                .rename_axis(index="names")
                .reset_index()
                .melt(id_vars="names", var_name="group", value_name=name)
            )
            d = d.merge(pts_df)

    # remove group column for backward compat if len(group) == 1
    if len(group) == 1:
        d.drop(columns="group", inplace=True)

    return d.reset_index(drop=True)


def _check_indices(
    dim_df: pd.DataFrame,
    alt_index: pd.Index,
    *,
    dim: Literal["obs", "var"],
    keys: Iterable[str],
    alias_index: pd.Index | None = None,
    use_raw: bool = False,
) -> tuple[list[str], list[str], list[str]]:
    """Check indices for `obs_df` and `var_df`."""
    alt_repr = "adata.raw" if use_raw else "adata"

    alt_dim = ("obs", "var")[dim == "obs"]

    alias_name = None
    if alias_index is not None:
        alt_names = pd.Series(alt_index, index=alias_index)
        alias_name = alias_index.name
        alt_search_repr = f"{alt_dim}['{alias_name}']"
    else:
        alt_names = pd.Series(alt_index, index=alt_index)
        alt_search_repr = f"{alt_dim}_names"

    col_keys = []
    index_keys = []
    index_aliases = []
    not_found = []

    # check that adata.obs does not contain duplicated columns
    # if duplicated columns names are present, they will
    # be further duplicated when selecting them.
    if not dim_df.columns.is_unique:
        dup_cols = dim_df.columns[dim_df.columns.duplicated()].tolist()
        msg = (
            f"adata.{dim} contains duplicated columns. Please rename or remove "
            "these columns first.\n`"
            f"Duplicated columns {dup_cols}"
        )
        raise ValueError(msg)

    if not alt_index.is_unique:
        msg = (
            f"{alt_repr}.{alt_dim}_names contains duplicated items\n"
            f"Please rename these {alt_dim} names first for example using "
            f"`adata.{alt_dim}_names_make_unique()`"
        )
        raise ValueError(msg)

    # use only unique keys, otherwise duplicated keys will
    # further duplicate when reordering the keys later in the function
    for key in dict.fromkeys(keys):
        if key in dim_df.columns:
            col_keys.append(key)
            if key in alt_names.index:
                msg = f"The key {key!r} is found in both adata.{dim} and {alt_repr}.{alt_search_repr}."
                raise KeyError(msg)
        elif key in alt_names.index:
            val = alt_names[key]
            if isinstance(val, pd.Series):
                # while var_names must be unique, adata.var[gene_symbols] does not
                # It's still ambiguous to refer to a duplicated entry though.
                assert alias_index is not None
                msg = f"Found duplicate entries for {key!r} in {alt_repr}.{alt_search_repr}."
                raise KeyError(msg)
            index_keys.append(val)
            index_aliases.append(key)
        else:
            not_found.append(key)
    if len(not_found) > 0:
        msg = (
            f"Could not find keys {not_found!r} in columns of `adata.{dim}` or in"
            f" {alt_repr}.{alt_search_repr}."
        )
        raise KeyError(msg)

    return col_keys, index_keys, index_aliases


def _get_array_values(
    X,
    dim_names: pd.Index,
    keys: Iterable[str],
    *,
    axis: Literal[0, 1],
    backed: bool,
):
    # TODO: This should be made easier on the anndata side
    mutable_idxer = [slice(None), slice(None)]
    idx = dim_names.get_indexer(keys)

    # for backed AnnData is important that the indices are ordered
    if backed:
        idx_order = np.argsort(idx)
        rev_idxer = mutable_idxer.copy()
        mutable_idxer[axis] = idx[idx_order]
        rev_idxer[axis] = np.argsort(idx_order)
        matrix = X[tuple(mutable_idxer)][tuple(rev_idxer)]
    else:
        mutable_idxer[axis] = idx
        matrix = X[tuple(mutable_idxer)]

    if isinstance(matrix, CSBase):
        matrix = matrix.toarray()

    return matrix


def obs_df(
    adata: AnnData,
    keys: Collection[str] = (),
    obsm_keys: Iterable[tuple[str, int]] = (),
    *,
    layer: str | None = None,
    gene_symbols: str | None = None,
    use_raw: bool = False,
) -> pd.DataFrame:
    """Return values for observations in adata.

    Params
    ------
    adata
        AnnData object to get values from.
    keys
        Keys from either `.var_names`, `.var[gene_symbols]`, or `.obs.columns`.
    obsm_keys
        Tuples of `(key from obsm, column index of obsm[key])`.
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

    >>> import scanpy as sc
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> plotdf = sc.get.obs_df(
    ...     pbmc, keys=["CD8B", "n_genes"], obsm_keys=[("X_umap", 0), ("X_umap", 1)]
    ... )
    >>> plotdf.columns
    Index(['CD8B', 'n_genes', 'X_umap-0', 'X_umap-1'], dtype='object')
    >>> plotdf.plot.scatter("X_umap-0", "X_umap-1", c="CD8B")  # doctest: +SKIP
    <Axes: xlabel='X_umap-0', ylabel='X_umap-1'>

    Calculating mean expression for marker genes by cluster:

    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> marker_genes = ["CD79A", "MS4A1", "CD8A", "CD8B", "LYZ"]
    >>> genedf = sc.get.obs_df(pbmc, keys=["louvain", *marker_genes])
    >>> grouped = genedf.groupby("louvain", observed=True)
    >>> mean, var = grouped.mean(), grouped.var()

    """
    if isinstance(keys, str):
        keys = [keys]
    if use_raw:
        assert layer is None, (
            "Cannot specify use_raw=True and a layer at the same time."
        )
        var = adata.raw.var
    else:
        var = adata.var
    alias_index = pd.Index(var[gene_symbols]) if gene_symbols is not None else None

    obs_cols, var_idx_keys, var_symbols = _check_indices(
        adata.obs,
        var.index,
        dim="obs",
        keys=keys,
        alias_index=alias_index,
        use_raw=use_raw,
    )

    # Make df
    df = pd.DataFrame(index=adata.obs_names)

    # add var values
    if len(var_idx_keys) > 0:
        matrix = _get_array_values(
            _get_obs_rep(adata, layer=layer, use_raw=use_raw),
            var.index,
            var_idx_keys,
            axis=1,
            backed=adata.isbacked,
        )
        df = pd.concat(
            [df, pd.DataFrame(matrix, columns=var_symbols, index=adata.obs_names)],
            axis=1,
        )

    # add obs values
    if len(obs_cols) > 0:
        df = pd.concat([df, adata.obs[obs_cols]], axis=1)

    # reorder columns to given order (including duplicates keys if present)
    if keys:
        df = df[keys]

    for k, idx in obsm_keys:
        added_k = f"{k}-{idx}"
        val = adata.obsm[k]
        if isinstance(val, np.ndarray):
            df[added_k] = np.ravel(val[:, idx])
        elif isinstance(val, CSBase):
            df[added_k] = np.ravel(val[:, idx].toarray())
        elif isinstance(val, pd.DataFrame):
            df[added_k] = val.loc[:, idx]

    return df


def var_df(
    adata: AnnData,
    keys: Collection[str] = (),
    varm_keys: Iterable[tuple[str, int]] = (),
    *,
    layer: str | None = None,
) -> pd.DataFrame:
    """Return values for observations in adata.

    Params
    ------
    adata
        AnnData object to get values from.
    keys
        Keys from either `.obs_names`, or `.var.columns`.
    varm_keys
        Tuples of `(key from varm, column index of varm[key])`.
    layer
        Layer of `adata` to use as expression values.

    Returns
    -------
    A dataframe with `adata.var_names` as index, and values specified by `keys`
    and `varm_keys`.

    """
    # Argument handling
    if isinstance(keys, str):
        keys = [keys]
    var_cols, obs_idx_keys, _ = _check_indices(
        adata.var, adata.obs_names, dim="var", keys=keys
    )

    # initialize df
    df = pd.DataFrame(index=adata.var.index)

    if len(obs_idx_keys) > 0:
        matrix = _get_array_values(
            _get_obs_rep(adata, layer=layer),
            adata.obs_names,
            obs_idx_keys,
            axis=0,
            backed=adata.isbacked,
        ).T
        df = pd.concat(
            [df, pd.DataFrame(matrix, columns=obs_idx_keys, index=adata.var_names)],
            axis=1,
        )

    # add obs values
    if len(var_cols) > 0:
        df = pd.concat([df, adata.var[var_cols]], axis=1)

    # reorder columns to given order
    if keys:
        df = df[keys]

    for k, idx in varm_keys:
        added_k = f"{k}-{idx}"
        val = adata.varm[k]
        if isinstance(val, np.ndarray):
            df[added_k] = np.ravel(val[:, idx])
        elif isinstance(val, CSBase):
            df[added_k] = np.ravel(val[:, idx].toarray())
        elif isinstance(val, pd.DataFrame):
            df[added_k] = val.loc[:, idx]
    return df


def _get_obs_rep(
    adata: AnnData,
    *,
    use_raw: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    obsp: str | None = None,
) -> (
    np.ndarray | CSBase | pd.DataFrame | ArrayView | BaseCompressedSparseDataset | None
):
    """Choose array aligned with obs annotation."""
    # https://github.com/scverse/scanpy/issues/1546
    if not isinstance(use_raw, bool):
        msg = f"use_raw expected to be bool, was {type(use_raw)}."
        raise TypeError(msg)

    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made in {0, 1}
    if choices_made == 0:
        return adata.X
    if is_layer:
        return adata.layers[layer]
    if use_raw:
        return adata.raw.X
    if is_obsm:
        return adata.obsm[obsm]
    if is_obsp:
        return adata.obsp[obsp]
    msg = (
        "That was unexpected. Please report this bug at:\n\n\t"
        "https://github.com/scverse/scanpy/issues"
    )
    raise AssertionError(msg)


def _set_obs_rep(
    adata: AnnData,
    val: Any,
    *,
    use_raw: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    obsp: str | None = None,
):
    """Set value for observation rep."""
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
        msg = (
            "That was unexpected. Please report this bug at:\n\n"
            "\thttps://github.com/scverse/scanpy/issues"
        )
        raise AssertionError(msg)


M = TypeVar("M", bound=NDArray[np.bool_] | NDArray[np.floating] | pd.Series | None)


def _check_mask(
    data: AnnData | np.ndarray | CSBase | DaskArray,
    mask: str | M,
    dim: Literal["obs", "var"],
    *,
    allow_probabilities: bool = False,
) -> M:  # Could also be a series, but should be one or the other
    """Validate mask argument.

    Params
    ------
    data
        Annotated data matrix or numpy array.
    mask
        Mask (or probabilities if `allow_probabilities=True`).
        Either an appropriatley sized array, or name of a column.
    dim
        The dimension being masked.
    allow_probabilities
        Whether to allow probabilities as `mask`
    """
    if mask is None:
        return mask
    desc = "mask/probabilities" if allow_probabilities else "mask"

    if isinstance(mask, str):
        if not isinstance(data, AnnData):
            msg = f"Cannot refer to {desc} with string without providing anndata object as argument"
            raise ValueError(msg)

        annot: pd.DataFrame = getattr(data, dim)
        if mask not in annot.columns:
            msg = (
                f"Did not find `adata.{dim}[{mask!r}]`. "
                f"Either add the {desc} first to `adata.{dim}`"
                f"or consider using the {desc} argument with an array."
            )
            raise ValueError(msg)
        mask_array = annot[mask].to_numpy()
    else:
        if len(mask) != data.shape[0 if dim == "obs" else 1]:
            msg = f"The shape of the {desc} do not match the data."
            raise ValueError(msg)
        mask_array = mask

    is_bool = pd.api.types.is_bool_dtype(mask_array.dtype)
    if not allow_probabilities and not is_bool:
        msg = "Mask array must be boolean."
        raise ValueError(msg)
    elif allow_probabilities and not (
        is_bool or pd.api.types.is_float_dtype(mask_array.dtype)
    ):
        msg = f"{desc} array must be boolean or floating point."
        raise ValueError(msg)

    return mask_array

from __future__ import annotations

import warnings
from functools import partial
from math import sqrt
from typing import TYPE_CHECKING

import numba
import numpy as np
import pandas as pd
from anndata import AnnData

from ... import logging as logg
from ..._compat import CSBase, njit
from ..._settings import Verbosity, settings
from ..._utils import _doc_params, check_nonnegative_integers, view_to_actual
from ...experimental._docs import (
    doc_adata,
    doc_check_values,
    doc_dist_params,
    doc_genes_batch_chunk,
    doc_inplace,
    doc_layer,
)
from ...get import _get_obs_rep
from ...preprocessing._distributed import materialize_as_ndarray
from ...preprocessing._utils import _get_mean_var

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray


@njit
def _calculate_res_sparse(
    indptr: NDArray[np.integer],
    index: NDArray[np.integer],
    data: NDArray[np.float64],
    *,
    sums_genes: NDArray[np.float64],
    sums_cells: NDArray[np.float64],
    sum_total: np.float64,
    clip: np.float64,
    theta: np.float64,
    n_genes: int,
    n_cells: int,
) -> NDArray[np.float64]:
    def get_value(cell: int, sparse_idx: int, stop_idx: int) -> np.float64:
        """Return the value at the specified cell location if it exists, or zero otherwise."""
        if sparse_idx < stop_idx and index[sparse_idx] == cell:
            return data[sparse_idx]
        else:
            return np.float64(0.0)

    def clac_clipped_res_sparse(gene: int, cell: int, value: np.float64) -> np.float64:
        mu = sums_genes[gene] * sums_cells[cell] / sum_total
        mu_sum = value - mu
        pre_res = mu_sum / sqrt(mu + mu * mu / theta)
        res = np.float64(min(max(pre_res, -clip), clip))
        return res

    residuals = np.zeros(n_genes, dtype=np.float64)
    for gene in numba.prange(n_genes):
        start_idx = indptr[gene]
        stop_idx = indptr[gene + 1]

        sparse_idx = start_idx
        var_sum = np.float64(0.0)
        sum_clipped_res = np.float64(0.0)
        for cell in range(n_cells):
            value = get_value(cell, sparse_idx, stop_idx)
            clipped_res = clac_clipped_res_sparse(gene, cell, value)
            if value > 0:
                sparse_idx += 1
            sum_clipped_res += clipped_res

        mean_clipped_res = sum_clipped_res / n_cells
        sparse_idx = start_idx
        for cell in range(n_cells):
            value = get_value(cell, sparse_idx, stop_idx)
            clipped_res = clac_clipped_res_sparse(gene, cell, value)
            if value > 0:
                sparse_idx += 1
            diff = clipped_res - mean_clipped_res
            var_sum += diff * diff

        residuals[gene] = var_sum / n_cells
    return residuals


@njit
def _calculate_res_dense(
    matrix,
    *,
    sums_genes: NDArray[np.float64],
    sums_cells: NDArray[np.float64],
    sum_total: np.float64,
    clip: np.float64,
    theta: np.float64,
    n_genes: int,
    n_cells: int,
) -> NDArray[np.float64]:
    def clac_clipped_res_dense(gene: int, cell: int) -> np.float64:
        mu = sums_genes[gene] * sums_cells[cell] / sum_total
        value = matrix[cell, gene]

        mu_sum = value - mu
        pre_res = mu_sum / sqrt(mu + mu * mu / theta)
        res = np.float64(min(max(pre_res, -clip), clip))
        return res

    residuals = np.zeros(n_genes, dtype=np.float64)

    for gene in numba.prange(n_genes):
        sum_clipped_res = np.float64(0.0)
        for cell in range(n_cells):
            sum_clipped_res += clac_clipped_res_dense(gene, cell)
        mean_clipped_res = sum_clipped_res / n_cells

        var_sum = np.float64(0.0)
        for cell in range(n_cells):
            clipped_res = clac_clipped_res_dense(gene, cell)
            diff = clipped_res - mean_clipped_res
            var_sum += diff * diff

        residuals[gene] = var_sum / n_cells
    return residuals


def _highly_variable_pearson_residuals(  # noqa: PLR0912, PLR0915
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_top_genes: int = 1000,
    batch_key: str | None = None,
    chunksize: int = 1000,
    check_values: bool = True,
    layer: str | None = None,
    subset: bool = False,
    inplace: bool = True,
) -> pd.DataFrame | None:
    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer)
    computed_on = layer if layer else "adata.X"

    # Check for raw counts
    if check_values and not check_nonnegative_integers(X):
        warnings.warn(
            "`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
            UserWarning,
            stacklevel=3,
        )
    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        msg = "Pearson residuals require theta > 0"
        raise ValueError(msg)
    # prepare clipping

    if batch_key is None:
        batch_info = np.zeros(adata.shape[0], dtype=int)
    else:
        batch_info = adata.obs[batch_key].values
    n_batches = len(np.unique(batch_info))

    # Get pearson residuals for each batch separately
    residual_gene_vars = []
    for batch in np.unique(batch_info):
        adata_subset_prefilter = adata[batch_info == batch]
        X_batch_prefilter = _get_obs_rep(adata_subset_prefilter, layer=layer)

        # Filter out zero genes
        with settings.verbosity.override(Verbosity.error):
            nonzero_genes = np.ravel(X_batch_prefilter.sum(axis=0)) != 0
        adata_subset = adata_subset_prefilter[:, nonzero_genes]
        X_batch = _get_obs_rep(adata_subset, layer=layer)

        # Prepare clipping
        if clip is None:
            n = X_batch.shape[0]
            clip = np.sqrt(n)
        if clip < 0:
            msg = "Pearson residuals require `clip>=0` or `clip=None`."
            raise ValueError(msg)

        if isinstance(X_batch, CSBase):
            X_batch = X_batch.tocsc()
            X_batch.eliminate_zeros()
            calculate_res = partial(
                _calculate_res_sparse,
                X_batch.indptr,
                X_batch.indices,
                X_batch.data.astype(np.float64),
            )
        else:
            X_batch = np.array(X_batch, dtype=np.float64, order="F")
            calculate_res = partial(_calculate_res_dense, X_batch)

        sums_genes = np.array(X_batch.sum(axis=0)).ravel()
        sums_cells = np.array(X_batch.sum(axis=1)).ravel()
        sum_total = np.sum(sums_genes)

        residual_gene_var = calculate_res(
            sums_genes=sums_genes,
            sums_cells=sums_cells,
            sum_total=np.float64(sum_total),
            clip=np.float64(clip),
            theta=np.float64(theta),
            n_genes=X_batch.shape[1],
            n_cells=X_batch.shape[0],
        )

        # Add 0 values for genes that were filtered out
        unmasked_residual_gene_var = np.zeros(len(nonzero_genes))
        unmasked_residual_gene_var[nonzero_genes] = residual_gene_var
        residual_gene_vars.append(unmasked_residual_gene_var.reshape(1, -1))

    residual_gene_vars = np.concatenate(residual_gene_vars, axis=0)

    # Get rank per gene within each batch
    # argsort twice gives ranks, small rank means most variable
    ranks_residual_var = np.argsort(np.argsort(-residual_gene_vars, axis=1), axis=1)
    ranks_residual_var = ranks_residual_var.astype(np.float32)
    # count in how many batches a genes was among the n_top_genes
    highly_variable_nbatches = np.sum(
        (ranks_residual_var < n_top_genes).astype(int), axis=0
    )
    # set non-top genes within each batch to nan
    ranks_residual_var[ranks_residual_var >= n_top_genes] = np.nan
    ranks_masked_array = np.ma.masked_invalid(ranks_residual_var)
    # Median rank across batches, ignoring batches in which gene was not selected
    medianrank_residual_var = np.ma.median(ranks_masked_array, axis=0).filled(np.nan)

    means, variances = materialize_as_ndarray(_get_mean_var(X))
    df = pd.DataFrame.from_dict(
        dict(
            means=means,
            variances=variances,
            residual_variances=np.mean(residual_gene_vars, axis=0),
            highly_variable_rank=medianrank_residual_var,
            highly_variable_nbatches=highly_variable_nbatches.astype(np.int64),
            highly_variable_intersection=highly_variable_nbatches == n_batches,
        )
    )
    df = df.set_index(adata.var_names)

    # Sort genes by how often they selected as hvg within each batch and
    # break ties with median rank of residual variance across batches
    df.sort_values(
        ["highly_variable_nbatches", "highly_variable_rank"],
        ascending=[False, True],
        na_position="last",
        inplace=True,
    )

    high_var = np.zeros(df.shape[0], dtype=bool)
    high_var[:n_top_genes] = True
    df["highly_variable"] = high_var
    df = df.loc[adata.var_names, :]

    if inplace:
        adata.uns["hvg"] = {"flavor": "pearson_residuals", "computed_on": computed_on}
        logg.hint(
            "added\n"
            "    'highly_variable', boolean vector (adata.var)\n"
            "    'highly_variable_rank', float vector (adata.var)\n"
            "    'highly_variable_nbatches', int vector (adata.var)\n"
            "    'highly_variable_intersection', boolean vector (adata.var)\n"
            "    'means', float vector (adata.var)\n"
            "    'variances', float vector (adata.var)\n"
            "    'residual_variances', float vector (adata.var)"
        )
        adata.var["means"] = df["means"].values
        adata.var["variances"] = df["variances"].values
        adata.var["residual_variances"] = df["residual_variances"]
        adata.var["highly_variable_rank"] = df["highly_variable_rank"].values
        if batch_key is not None:
            adata.var["highly_variable_nbatches"] = df[
                "highly_variable_nbatches"
            ].values
            adata.var["highly_variable_intersection"] = df[
                "highly_variable_intersection"
            ].values
        adata.var["highly_variable"] = df["highly_variable"].values

        if subset:
            adata._inplace_subset_var(df["highly_variable"].values)

    else:
        if batch_key is None:
            df = df.drop(
                ["highly_variable_nbatches", "highly_variable_intersection"], axis=1
            )
        if subset:
            df = df.iloc[df.highly_variable.values, :]

        return df


@_doc_params(
    adata=doc_adata,
    dist_params=doc_dist_params,
    genes_batch_chunk=doc_genes_batch_chunk,
    check_values=doc_check_values,
    layer=doc_layer,
    inplace=doc_inplace,
)
def highly_variable_genes(  # noqa: PLR0913
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_top_genes: int | None = None,
    batch_key: str | None = None,
    chunksize: int = 1000,
    flavor: Literal["pearson_residuals"] = "pearson_residuals",
    check_values: bool = True,
    layer: str | None = None,
    subset: bool = False,
    inplace: bool = True,
) -> pd.DataFrame | None:
    """Select highly variable genes using analytic Pearson residuals :cite:p:`Lause2021`.

    In :cite:t:`Lause2021`, Pearson residuals of a negative binomial offset model are computed
    (with overdispersion `theta` shared across genes). By default, overdispersion
    `theta=100` is used and residuals are clipped to `sqrt(n_obs)`. Finally, genes
    are ranked by residual variance.

    Expects raw count input.

    Parameters
    ----------
    {adata}
    {dist_params}
    {genes_batch_chunk}
    flavor
        Choose the flavor for identifying highly variable genes. In this experimental
        version, only 'pearson_residuals' is functional.
    {check_values}
    {layer}
    subset
        If `True`, subset the data to highly-variable genes after finding them.
        Otherwise merely indicate highly variable genes in `adata.var` (see below).
    {inplace}

    Returns
    -------
    If `inplace=True`, `adata.var` is updated with the following fields. Otherwise,
    returns the same fields as :class:`~pandas.DataFrame`.

    highly_variable : :class:`bool`
        boolean indicator of highly-variable genes.
    means : :class:`float`
        means per gene.
    variances : :class:`float`
        variance per gene.
    residual_variances : :class:`float`
        For `flavor='pearson_residuals'`, residual variance per gene. Averaged in the
        case of multiple batches.
    highly_variable_rank : :class:`float`
        For `flavor='pearson_residuals'`, rank of the gene according to residual.
        variance, median rank in the case of multiple batches.
    highly_variable_nbatches : :class:`int`
        If `batch_key` given, denotes in how many batches genes are detected as HVG.
    highly_variable_intersection : :class:`bool`
        If `batch_key` given, denotes the genes that are highly variable in all batches.

    Notes
    -----
    Experimental version of `sc.pp.highly_variable_genes()`

    """
    logg.info("extracting highly variable genes")

    if not isinstance(adata, AnnData):
        msg = (
            "`pp.highly_variable_genes` expects an `AnnData` argument, "
            "pass `inplace=False` if you want to return a `pd.DataFrame`."
        )
        raise ValueError(msg)

    if flavor == "pearson_residuals":
        if n_top_genes is None:
            msg = (
                "`pp.highly_variable_genes` requires the argument `n_top_genes`"
                " for `flavor='pearson_residuals'`"
            )
            raise ValueError(msg)
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
    else:
        msg = "This is an experimental API and only `flavor=pearson_residuals` is available."
        raise ValueError(msg)

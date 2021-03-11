"""Moran's I global spatial autocorrelation."""
from typing import Any, Union, Optional, Sequence
from anndata import AnnData

import numpy as np
import pandas as pd
import numba.types as nt

from numba import njit
from statsmodels.stats.multitest import multipletests

it = nt.int64
ft = nt.float64
ip = np.int64
fp = np.float64


def morans_i(
    adata: AnnData,
    connectivity_key: str = "connectivities",
    genes: Optional[Union[str, Sequence[str]]] = None,
    n_perms: int = 100,
    corr_method: Optional[str] = "fdr_bh",
    layer: Optional[str] = None,
    copy: bool = False,
) -> Optional[pd.DataFrame]:
    """
    Calculate Moran’s I Global Autocorrelation Statistic.

    The statistics is described here TODO.

    Parameters
    ----------
    adata
    connectivity_key
    genes
        List of gene names, as stored in :attr:`anndata.AnnData.var_names`, used to compute Moran's I statistics.
        If `None`, it's computed :attr:`anndata.AnnData.var` ``['highly_variable']``, if present. Otherwise,
        it's computed for all genes.
    n_perms
        Number of permutations.
    corr_method
        FDR correction method.
    layer
        Layer in :attr:`anndata.AnnData.layers` to use. If `None`, use :attr:`anndata.AnnData.X`.
    %(seed)s
    %(copy)s
    %(parallelize)s

    Returns
    -------
    If ``copy = True``, returns a :class:`pandas.DataFrame` with the following keys:

        - `'I'` - Moran's I statistic.
        - `'pval_sim'` - p-value based on permutations.
        - `'VI_sim'` - variance of `'I'` from permutations.
        - `'pval_sim_{{corr_method}}'` - the corrected p-values if ``corr_method != None`` .

    Otherwise, modifies the ``adata`` with the following key:

        - :attr:`anndata.AnnData.uns` ``['moranI']`` - the above mentioned dataframe.
    """

    if genes is None:
        if "highly_variable" in adata.var.columns:
            genes = adata[:, adata.var.highly_variable.values].var_names.values
        else:
            genes = adata.var_names.values

    adj = adata.obsp[connectivity_key]

    moran_list = []
    indptr = adj.indptr.astype(ip)
    indices = adj.indices.astype(ip)
    data = adj.data.astype(fp)

    for g in genes:
        counts = adata.obs_vector(g, layer=layer).astype(fp, copy=True)
        mi = _moran_score_perms(counts, indptr, indices, data, n_perms)
        moran_list.append(mi)

    df = pd.DataFrame(moran_list, columns=["I", "pval_sim", "VI_sim"], index=genes)

    if corr_method is not None:
        _, pvals_adj, _, _ = multipletests(
            df["pval_sim"].values, alpha=0.05, method=corr_method
        )
        df[f"pval_sim_{corr_method}"] = pvals_adj

    df.sort_values(by="I", ascending=False, inplace=True)

    if copy:
        return df
    adata.uns["moranI"] = df


@njit(
    ft(ft[:], it[:], it[:], ft[:], ft[:], ft, ft),
    parallel=False,
    fastmath=True,
)
def _compute_moran(
    counts: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    z: np.ndarray,
    data_sum: fp,
    z2ss: fp,
) -> Any:

    zl = np.empty(len(z))
    N = len(indptr) - 1
    for i in range(N):
        s = slice(indptr[i], indptr[i + 1])
        i_indices = indices[s]
        i_data = data[s]
        zl[i] = np.sum(i_data * z[i_indices])
    inum = (z * zl).sum()

    return len(counts) / data_sum * inum / z2ss


@njit(
    ft[:](ft[:], it[:], it[:], ft[:], it),
    parallel=False,
    fastmath=True,
)
def _moran_score_perms(
    counts: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    n_perms: ip,
) -> np.ndarray:

    perms = np.empty(n_perms, dtype=ft)
    res = np.empty(3, dtype=ft)

    z = counts - counts.mean()
    data_sum = data.sum()
    z2ss = (z * z).sum()

    res[0] = _compute_moran(counts, indptr, indices, data, z, data_sum, z2ss)

    for p in range(len(perms)):
        np.random.shuffle(z)
        perms[p] = _compute_moran(counts, indptr, indices, data, z, data_sum, z2ss)
    res[1] = (np.sum(perms > res[0]) + 1) / (n_perms + 1)
    res[2] = np.var(perms)
    return res

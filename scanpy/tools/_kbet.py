from typing import Union, Tuple
from warnings import warn

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, OptimizeResult
from scipy.sparse import spmatrix
from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests
from anndata import AnnData

from ..neighbors import neighbors, Neighbors


def kbet(
    adata: AnnData,
    batch_key: str = 'batch',
    *,
    alpha: float = .05,
    adjacency: spmatrix = None,
    copy: bool = False
) -> Union[AnnData, Tuple[float, np.ndarray]]:
    """kBET: k-nearest neighbour batch effect test.

    Use the heuristic :func:`scanpy.api.pp.kbet_neighbors` to use the ideal
    neighborhood size in :func:`scanpy.api.pp.neighbors`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    batch_key
        The column in :attr:`anndata.AnnData.uns` to use as batch ID.
    alpha
        family-wise error rate. The p-values for the χ²-test that are
        > ``alpha`` are used for the rejection rate.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        ``adata.uns['neighbors']['connectivities']``.
    copy
        Copy instance before computation and return a copy.
        Otherwise, perform computation in-place and return acceptance
        rate and p-values.

    Returns
    -------
    ``adata``
        If ``copy == True``, a copy of the input ``adata`` will be returned.
        Else, ``adata.uns['kbet']`` will be a tuple, see below:
    ``acceptance``
        If ``copy == False``, the acceptance rate (1-rejection) is returned.
    ``p_values``
        If ``copy == False``, the second returned value is the per-cell corrected p-values.
    """
    if adjacency is None:
        nbs = adata.uns.get('neighbors')
        if nbs is None:
            raise ValueError('No neighbors found. Provide the `adjacency` parameter or run `sc.pp.neighbors(adata)`')
        adjacency = nbs['connectivities']  # type: spmatrix
    adjacency = adjacency.tocsr()

    n_obs = adata.n_obs
    batch_ids = pd.Categorical(adata.obs[batch_key])
    # dof = len(batch_ids.unique()) - 1

    freqs_all = batch_ids.value_counts().sort_index() / len(batch_ids)
    freqs_neighbors = np.ndarray((n_obs, len(batch_ids.categories)))

    mask = adjacency != 0
    # cat_2d = np.tile(batch_ids, (n_obs, 1))
    # mapped = np.where(mask.A, cat_2d, None)
    for obs_i in range(n_obs):
        row_idx = mask[:, obs_i].A.flatten()
        freqs_obs = batch_ids[row_idx].value_counts().sort_index() / row_idx.sum()
        freqs_neighbors[obs_i, :] = freqs_obs

    _, p_vals_uncor = chisquare(freqs_neighbors, freqs_all, axis=1)
    rejected, p_vals, *_ = multipletests(p_vals_uncor, alpha)
    rate_rej = rejected.sum() / len(p_vals)

    rate_acc = 1 - rate_rej
    if copy:
        ad_ret = adata.copy()
        ad_ret.uns['kbet'] = rate_acc
        ad_ret.obs['kbet'] = p_vals
        return ad_ret
    else:
        return rate_acc, p_vals


def kbet_neighbors(
    adata: AnnData,
    batch_key: str = 'batch',
    *,
    return_k: bool = False,
    copy: bool = False,
    **neighbors_args
) -> Union[int, AnnData, None]:
    """A heuristic for kBET neighborhood based on batch sizes.

    Parameters
    ----------
    adata
        Annotated data object
    batch_key
        The column in :attr:`anndata.AnnData.uns` to use as batch ID.
    return_k
        Return the neighborhood size instead of calling :func:`scanpy.api.pp.neighbors`.
    copy
        If ``return_n == False``, return a copy of ``adata`` with the new neighbors set or
        just set ``adata.uns['neighbors']`` and return ``None``?
    neighbors_args
        Arguments passed to :func:`scanpy.api.pp.neighbors`.

    Returns
    -------
    ``n_neighbors``
        If ``return_n=True`` is set, return an :class:`int` with the
        heuristically determined neighborhood size.
    ``adata`` or ``None``
        Depends on ``copy``. See :func:`scanpy.api.pp.neighbors`.
    """
    if return_k and copy:
        raise ValueError(
            'You passed both return_n=True and copy=True to `kbet_neighbors`. '
            'I don’t know what to make of that.'
        )

    batch_ids = pd.Categorical(adata.obs[batch_key])
    batch_sizes = batch_ids.value_counts()

    k0 = np.floor(np.mean(batch_sizes) * (3/4))
    res = minimize_scalar(
        _kbet_step,
        bounds=[10, k0],
        args=(adata, batch_key),
        method='Bounded',
        # tol=.5,  # IDK how to minimize step size
    )  # type: OptimizeResult
    if res.success:
        k = int(res.x)
    else:
        warn('Could not find a good neighborhood size: %s', res)
        k = np.floor(np.mean(batch_sizes) / 4)

    if return_k:
        return k

    if 'neighbors' in adata.uns:
        warn('`adata.uns` already contains neighbors. `kbet_neighbors` will overwrite them.')

    return neighbors(adata, k, **neighbors_args, copy=copy)


def _kbet_step(k: int, adata: AnnData, batch_key: str):
    neighs = Neighbors(adata)
    neighs.compute_neighbors(n_neighbors=int(k))
    acc, _ = kbet(adata, batch_key, adjacency=neighs.connectivities)
    return acc

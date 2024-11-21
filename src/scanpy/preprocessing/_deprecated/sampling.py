from __future__ import annotations

from typing import TYPE_CHECKING

from ..._compat import _legacy_numpy_gen, old_positionals
from .._simple import sample

if TYPE_CHECKING:
    import numpy as np
    from anndata import AnnData
    from numpy.typing import NDArray
    from scipy.sparse import csc_matrix, csr_matrix

    from ..._compat import _LegacyRandom

    CSMatrix = csr_matrix | csc_matrix


@old_positionals("n_obs", "random_state", "copy")
def subsample(
    data: AnnData | np.ndarray | CSMatrix,
    fraction: float | None = None,
    *,
    n_obs: int | None = None,
    random_state: _LegacyRandom = 0,
    copy: bool = False,
) -> AnnData | tuple[np.ndarray | CSMatrix, NDArray[np.int64]] | None:
    """\
    Subsample to a fraction of the number of observations.

    .. deprecated:: 1.11.0

       Use :func:`~scanpy.pp.sample` instead.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    fraction
        Subsample to this `fraction` of the number of observations.
    n_obs
        Subsample to this number of observations.
    random_state
        Random seed to change subsampling.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.

    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~anndata.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """

    rng = _legacy_numpy_gen(random_state)
    return sample(
        data=data, fraction=fraction, n=n_obs, rng=rng, copy=copy, replace=False, axis=0
    )

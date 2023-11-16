from __future__ import annotations

from collections.abc import Callable

import pytest
import numpy as np
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

from scanpy._utils.compute.is_constant import is_constant
from scanpy.neighbors._common import (
    _ind_dist_shortcut,
    _get_sparse_matrix_from_indices_distances,
)


def mk_knn_matrix(
    n_obs: int,
    n_neighbors: int,
    *,
    plus_one: bool = False,
    includes_self: bool = False,
    duplicates: bool = False,
) -> sparse.csr_matrix:
    if plus_one:
        n_neighbors += 1
    dists = np.random.randn(n_obs, n_neighbors)
    if includes_self:
        dists[:, 0] = 0.0
    if duplicates:
        # Don’t use the first column, as that can be the cell itself
        dists[n_obs // 4 : n_obs, 2] = 0.0
    idxs = np.arange(n_obs * n_neighbors).reshape((n_obs, n_neighbors))
    mat = _get_sparse_matrix_from_indices_distances(
        idxs, dists, n_obs=n_obs, n_neighbors=n_neighbors
    )
    if duplicates:
        # Make sure the actual matrix has a regular sparsity pattern
        assert is_constant(mat.getnnz(axis=1))
        # Make sure implicit zeros for duplicates would change the sparsity pattern
        mat_sparsified = mat.copy()
        mat_sparsified.eliminate_zeros()
        assert not is_constant(mat_sparsified.getnnz(axis=1))
    return mat


@pytest.mark.parametrize("n_neighbors", [3, pytest.param(None, id="all")])
@pytest.mark.parametrize("plus_one", [True, False])
@pytest.mark.parametrize("includes_self", [True, False])
@pytest.mark.parametrize("duplicates", [True, False])
def test_ind_dist_shortcut_manual(
    n_neighbors: int | None, plus_one: bool, includes_self: bool, duplicates: bool
):
    n_obs = 500
    if n_neighbors is None:
        n_neighbors = n_obs - 1 if plus_one else n_obs

    mat = mk_knn_matrix(
        n_obs,
        n_neighbors,
        plus_one=plus_one,
        includes_self=includes_self,
        duplicates=duplicates,
    )

    assert (mat.nnz / n_obs) in {n_neighbors, n_neighbors + 1}
    assert _ind_dist_shortcut(mat, n_neighbors) is not None


@pytest.mark.parametrize("n_neighbors", [3, pytest.param(None, id="all")])
@pytest.mark.parametrize(
    "mk_mat",
    [
        pytest.param(
            lambda n_obs, n_neighbors: KNeighborsTransformer(
                n_neighbors=n_neighbors
            ).fit_transform(np.random.randn(n_obs, n_obs // 4)),
            id="sklearn_auto",
        )
    ],
)
def test_ind_dist_shortcut_premade(
    n_neighbors: int | None, mk_mat: Callable[[int, int], sparse.csr_matrix]
):
    n_obs = 500
    if n_neighbors is None:
        n_neighbors = n_obs - 1  # KNeighborsTransformer will add 1
    mat = mk_mat(n_obs, n_neighbors)

    assert (mat.nnz / n_obs) in {n_neighbors, n_neighbors + 1}
    assert _ind_dist_shortcut(mat, n_neighbors) is not None

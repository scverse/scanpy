from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from sklearn.neighbors import KNeighborsTransformer

from scanpy._utils.compute.is_constant import is_constant
from scanpy.neighbors._common import (
    _get_sparse_matrix_from_indices_distances,
    _has_self_column,
    _ind_dist_shortcut,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from scipy import sparse

    from scanpy._compat import CSRBase


def mk_knn_matrix(
    n_obs: int,
    n_neighbors: int,
    *,
    style: Literal["basic", "rapids", "sklearn"],
    duplicates: bool = False,
) -> CSRBase:
    n_col = n_neighbors + (1 if style == "sklearn" else 0)
    dists = np.abs(np.random.randn(n_obs, n_col)) + 1e-8
    idxs = np.arange(n_obs * n_col).reshape((n_col, n_obs)).T
    if style == "rapids":
        idxs[:, 0] += 1  # does not include cell itself
    else:
        dists[:, 0] = 0.0  # includes cell itself
    if duplicates:
        # Don’t use the first column, as that might be the cell itself
        dists[n_obs // 4 : n_obs, 2] = 0.0
    # keep self column to simulate output from kNN transformers
    mat = _get_sparse_matrix_from_indices_distances(idxs, dists, keep_self=True)

    # check if out helper here works as expected
    assert _has_self_column(idxs, dists) == (style != "rapids")
    if duplicates:
        # Make sure the actual matrix has a regular sparsity pattern
        assert is_constant(mat.getnnz(axis=1))
        # Make sure implicit zeros for duplicates would change the sparsity pattern
        mat_sparsified = mat.copy()
        mat_sparsified.eliminate_zeros()
        assert not is_constant(mat_sparsified.getnnz(axis=1))

    return mat


@pytest.mark.parametrize("n_neighbors", [3, pytest.param(None, id="all")])
@pytest.mark.parametrize("style", ["basic", "rapids", "sklearn"])
@pytest.mark.parametrize("duplicates", [True, False], ids=["duplicates", "unique"])
def test_ind_dist_shortcut_manual(
    *,
    n_neighbors: int | None,
    style: Literal["basic", "rapids", "sklearn"],
    duplicates: bool,
):
    n_obs = 10
    if n_neighbors is None:
        n_neighbors = n_obs
    mat = mk_knn_matrix(n_obs, n_neighbors, style=style, duplicates=duplicates)

    assert (mat.nnz / n_obs) == n_neighbors + (1 if style == "sklearn" else 0)
    assert _ind_dist_shortcut(mat) is not None


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
    n_neighbors: int | None,
    mk_mat: Callable[[int, int], sparse.csr_matrix],  # noqa: TID251
):
    n_obs = 10
    if n_neighbors is None:
        # KNeighborsTransformer interprets this as “number of neighbors excluding cell itself”
        # so it can be at most n_obs - 1
        n_neighbors = n_obs - 1
    mat = mk_mat(n_obs, n_neighbors)

    assert (mat.nnz / n_obs) == n_neighbors + 1
    assert _ind_dist_shortcut(mat) is not None

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pytest
from fast_array_utils.stats import is_constant
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

from scanpy.neighbors._common import (
    _get_indices_distances_from_sparse_matrix,
    _get_sparse_matrix_from_indices_distances,
    _has_self_column,
    _ind_dist_shortcut,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

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
    assert _has_self_column(idxs) == (style != "rapids")
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


def mk_variable_knn_matrix_self_inclusive(n_obs: int) -> CSRBase:
    """Create matrix with variable neighbors per row, including self with distance 0."""
    indices_list = [[0, 1], [1, 0, 2], [2, 1], [3, 0, 1]]
    distances_list = [[0.0, 0.5], [0.0, 0.3, 0.4], [0.0, 0.7], [0.0, 0.8, 0.6]]

    row_indices = []
    col_indices = []
    data = []

    if n_obs > len(indices_list):
        msg = (
            f"n_obs={n_obs} is too large for the predefined variable neighbor structure "
            f"with {len(indices_list)} rows."
        )
        raise ValueError(msg)

    for i in range(n_obs):
        for idx, dist in zip(indices_list[i], distances_list[i], strict=True):
            row_indices.append(i)
            col_indices.append(idx)
            data.append(dist)

    mat = sparse.csr_matrix(  # noqa: TID251
        (np.array(data), (np.array(row_indices), np.array(col_indices))),
        shape=(n_obs, n_obs),
    )
    return mat


def mk_variable_knn_matrix_self_exclusive(n_obs: int) -> CSRBase:
    """Create matrix with variable neighbors per row, excluding self (no distance 0)."""
    indices_list = [[1], [0, 2], [1], [0, 1]]
    distances_list = [[0.5], [0.3, 0.4], [0.7], [0.8, 0.6]]

    row_indices = []
    col_indices = []
    data = []

    if n_obs > len(indices_list):
        msg = (
            f"n_obs={n_obs} is too large for the predefined variable neighbor structure "
            f"with {len(indices_list)} rows."
        )
        raise ValueError(msg)

    for i in range(n_obs):
        for idx, dist in zip(indices_list[i], distances_list[i], strict=True):
            row_indices.append(i)
            col_indices.append(idx)
            data.append(dist)

    mat = sparse.csr_matrix(  # noqa: TID251
        (np.array(data), (np.array(row_indices), np.array(col_indices))),
        shape=(n_obs, n_obs),
    )
    return mat


@pytest.mark.parametrize(
    ("matrix_func", "self_inclusive"),
    [
        pytest.param(mk_variable_knn_matrix_self_inclusive, True, id="self_inclusive"),
        # pytest.param(mk_variable_knn_matrix_self_exclusive, False, id="self_exclusive"),
    ],
)
def test_variable_neighbors_uses_slow_path(matrix_func, self_inclusive):
    """Test variable neighbor counts trigger slow path with warning.

    This test mocks `sc.pp.neighbors`, which always validates the neighbors returned from
    the nearest neighbor search with `_get_indices_distances_from_sparse_matrix` and then
    converts them back to a sparse.

    Here we use `_get_sparse_matrix_from_indices_distances`, as is used by `method='binary'`
    in `sc.pp.neighbors`.
    """
    n_obs = 4
    n_neighbors = 3  # int(self_inclusive) + 2

    mat = matrix_func(n_obs)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        indices, distances = _get_indices_distances_from_sparse_matrix(
            mat,
            n_neighbors=n_neighbors,
        )
        assert len(w) == 1
        assert "no constant number of neighbors" in str(w[0].message)

    assert indices.shape == (n_obs, n_neighbors)

    # first column should always be self (distance 0)
    assert np.array_equal(indices[:, 0], np.arange(n_obs))
    assert np.allclose(distances[:, 0], 0.0)

    connectivities = _get_sparse_matrix_from_indices_distances(
        indices,
        distances,
        keep_self=False,
    )
    assert connectivities.shape == (n_obs, n_obs)

    # TODO: This is a hack, as we prefill distances with zeros
    # `_get_sparse_matrix_from_indices_distances` does not eliminate zeros itself
    # but rather searches for the self column and removes it if requested
    # We can also not rely on `_get_sparse_matrix_from_indices_distances` being
    # called, as it is specific to `method='binary'` in `sc.pp.neighbors`.
    connectivities.eliminate_zeros()

    # check different rows have different numbers of connections
    nnz_per_row = np.diff(connectivities.indptr)
    assert not is_constant(nnz_per_row), (
        f"Expected variable connectivity, got {nnz_per_row}"
    )

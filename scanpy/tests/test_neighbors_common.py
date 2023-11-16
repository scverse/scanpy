from functools import partial

import pytest
import numpy as np
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

from scanpy._utils.compute.is_constant import is_constant
from scanpy.neighbors._common import _ind_dist_shortcut


def mk_knn_matrix(
    n_obs: int,
    n_neighbors: int,
    *,
    plus_one: bool = False,
    duplicates: bool = False,
):
    if plus_one:
        n_neighbors += 1
    n_nonzero = n_obs * n_neighbors
    dists = np.random.randn(n_obs * n_neighbors)
    if duplicates:
        dists[n_obs // 4 : n_obs] = 0.0
    idxs = np.arange(n_obs * n_neighbors)
    indptr = np.arange(0, n_nonzero + 1, n_neighbors)
    mat = sparse.csr_matrix((dists, idxs, indptr), shape=(n_obs, n_obs))
    if duplicates:
        # Make sure the actual matrix has a regular sparsity pattern
        assert is_constant(mat.getnnz(axis=1))
        # Make sure implicit zeros for duplicates would change the sparsity pattern
        mat_sparsified = mat.copy()
        mat_sparsified.eliminate_zeros()
        assert not is_constant(mat_sparsified.getnnz(axis=1))
    return mat


@pytest.mark.parametrize("n_neighbors", [3, pytest.param(None, id="all")])
@pytest.mark.parametrize(
    "mk_mat",
    [
        pytest.param(partial(mk_knn_matrix), id="n"),
        pytest.param(partial(mk_knn_matrix, plus_one=True), id="n+1"),
        pytest.param(partial(mk_knn_matrix, duplicates=True), id="dupes"),
        pytest.param(
            lambda n_obs, n_neighbors: KNeighborsTransformer(
                n_neighbors=n_neighbors
            ).fit_transform(np.random.randn(n_obs, n_obs // 4)),
            id="sklearn_auto",
        ),
    ],
)
def test_ind_dist_shortcut(n_neighbors, mk_mat):
    n_obs = 500
    if n_neighbors is None:
        n_neighbors = n_obs - 1  # KNeighborsTransformer will add 1
    mat = mk_mat(n_obs, n_neighbors)

    assert (mat.nnz / n_obs) in {n_neighbors, n_neighbors + 1}
    assert _ind_dist_shortcut(mat, n_neighbors) is not None

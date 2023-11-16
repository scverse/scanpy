import numpy as np
import pytest
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

from scanpy.neighbors._common import _ind_dist_shortcut


@pytest.mark.parametrize(
    "mk_mat",
    [
        pytest.param(
            lambda d, n: KNeighborsTransformer(n_neighbors=n).fit_transform(d),
            id="sklearn_auto",
        ),
    ],
)
def test_ind_dist_shortcut(mk_mat):
    n_cells = 500
    d = sparse.rand(n_cells, 100, format="csr", density=0.2)
    n = 3
    mat = mk_mat(d, n)
    assert (mat.nnz / n_cells) in {n, n + 1}
    assert _ind_dist_shortcut(mat, n) is not None

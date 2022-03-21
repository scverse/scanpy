import scanpy as sc

import pytest


def test_deprecate_multicore_tsne():
    pbmc = sc.datasets.pbmc68k_reduced()

    with pytest.warns(
        UserWarning, match="calling tsne with n_jobs > 1 would use MulticoreTSNE"
    ):
        sc.tl.tsne(pbmc, n_jobs=2)

    with pytest.warns(FutureWarning, match="Argument `use_fast_tsne` is deprecated"):
        sc.tl.tsne(pbmc, use_fast_tsne=True)

    with pytest.warns(UserWarning, match="Falling back to scikit-learn"):
        sc.tl.tsne(pbmc, use_fast_tsne=True)

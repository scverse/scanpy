import pytest

import scanpy as sc
import scanpy.external as sce

pytest.importorskip("palantir")


def test_palantir_core():
    adata = sc.datasets.pbmc3k_processed()

    sce.tl.palantir(adata=adata, n_components=5, knn=30)
    assert adata.layers['palantir_imp'].shape[0], "palantir_imp matrix Error!"

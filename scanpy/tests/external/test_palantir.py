import scanpy.external as sce
from scanpy.testing._helpers.data import pbmc3k_processed
from scanpy.testing._pytest.marks import needs


pytestmark = [needs('palantir')]


def test_palantir_core():
    adata = pbmc3k_processed()

    sce.tl.palantir(adata=adata, n_components=5, knn=30)
    assert adata.layers['palantir_imp'].shape[0], "palantir_imp matrix Error!"

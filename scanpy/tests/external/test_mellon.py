import pytest

import scanpy as sc
import scanpy.external as sce
from scanpy.testing._helpers.data import pbmc3k_processed
from scanpy.testing._pytest.marks import needs


@needs("mellon")
def test_mellon_core():
    adata = pbmc3k_processed()

    sce.tl.mellon(adata=adata, repr_key="X_pca")
    assert adata.obs['mellon_log_density'].shape[0], "mellon cell-state density Error!"

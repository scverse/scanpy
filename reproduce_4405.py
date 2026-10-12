from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

import scanpy as sc

rng = np.random.default_rng(0)
X = rng.negative_binomial(3, 0.3, (300, 120)).astype(np.float32)

duplicate_names = [f"g{i % 100}" for i in range(120)]
unique_names = [f"g{i}" for i in range(120)]

for names in (unique_names, duplicate_names):
    adata = AnnData(X.copy(), var=pd.DataFrame(index=names))

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=30,
    )

    print("Selected genes:", adata.var["highly_variable"].sum())

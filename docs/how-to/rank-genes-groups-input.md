# Trace the input to a differential expression test

In Scanpy 1.13, `sc.tl.rank_genes_groups(..., method="t-test")` uses one expression matrix for both Welch’s t-test and the reported log2 fold change. If you omit `use`, Scanpy selects `adata.raw.X` when present and `adata.X` otherwise. You can select `X` or a named layer explicitly. The example uses log1p expression in all three places, including the raw snapshot.

Give one gene opposite patterns in the snapshot and the current matrix. This makes it easy to see which matrix Scanpy read.

```python
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

def log_expression(values):
    return np.log1p(values)[:, None]

adata = AnnData(log_expression([8, 10, 12, 1, 2, 3]))
adata.obs["condition"] = pd.Categorical(["A"] * 3 + ["B"] * 3)
adata.raw = adata.copy()
adata.X = log_expression([1, 2, 3, 8, 10, 12])
adata.layers["alt"] = log_expression([1, 5, 15, 2, 3, 4])

def compare(use=None, *, mean_in_log_space=False):
    sc.tl.rank_genes_groups(
        adata, "condition", groups=["A"], reference="B", method="t-test",
        use=use, mean_in_log_space=mean_in_log_space,
    )
    row = sc.get.rank_genes_groups_df(adata, group="A").iloc[0]
    return row[["scores", "pvals", "logfoldchanges"]].astype(float)

print(pd.DataFrame({
    "default (raw.X)": compare(),
    "X": compare("X"),
    "layers.alt": compare("layers.alt"),
    "layers.alt, mean in log space": compare("layers.alt", mean_in_log_space=True),
}).T.round(3))
```

| Input | Welch’s t | p-value | Log2 fold change |
|---|---:|---:|---:|
| Default with raw snapshot | 5.836 | 0.010 | 2.322 |
| X | −5.836 | 0.010 | −2.322 |
| alt layer | 0.627 | 0.589 | 1.222 |
| alt layer, mean in log space | 0.627 | 0.589 | 0.710 |

With `mean_in_log_space=False`, Scanpy reverses log1p before averaging for the fold change. In the snapshot, group A averages 10 and group B averages 2, giving log2(10/2) ≈ 2.322. Selecting `X` swaps those means and flips the signs of both the fold change and the test statistic.

The last two rows use the same layer. Changing `mean_in_log_space` leaves the test statistic and p-value unchanged because Welch’s test still uses the log1p values. It changes the fold change because, with `True`, Scanpy averages the log values first and then reverses log1p. Scanpy adds a small offset to each group mean before taking the ratio.

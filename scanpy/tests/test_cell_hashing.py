from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np

def test_cell_hashing():
    X = np.ones((10, 10))
    for x in range(10):
        X[x, x] = 1000
    X[0, 1] = 1000
    test_data = AnnData(X)
    sc.pp.demultiplex_hashing(test_data)

    ground_truth = [["Doublet", "0_1"]]
    for x in range(1, 10):
        ground_truth.append([str(x), str(x)])
    ground_truth_df = pd.DataFrame(ground_truth, columns=["ID", "CLASSIFICATION"])
    assert all(test_data.obs["ID"].values == ground_truth_df["ID"].values)
    assert all(test_data.obs["CLASSIFICATION"].values == ground_truth_df["CLASSIFICATION"].values)

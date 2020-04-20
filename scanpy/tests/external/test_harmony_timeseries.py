from itertools import product

import pytest
from anndata import AnnData

import scanpy as sc
import scanpy.external as sce

pytest.importorskip("harmony")


def test_load_timepoints_from_anndata_list():
    adata_ref = sc.datasets.pbmc3k()
    start = [596, 615, 1682, 1663, 1409, 1432]
    adata = AnnData.concatenate(
        *(adata_ref[i : i + 1000] for i in start),
        join="outer",
        batch_key="sample",
        batch_categories=[f"sa{i}_Rep{j}" for i, j in product((1, 2, 3), (1, 2))],
    )
    adata.obs["time_points"] = adata.obs["sample"].str.split("_", expand=True)[0]
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    sce.tl.harmony_timeseries(adata=adata, tp="time_points", n_components=None)
    assert all(
        [adata.obsp['harmony_aff'].shape[0], adata.obsp['harmony_aff_aug'].shape[0]]
    ), "harmony_timeseries augmented affinity matrix Error!"

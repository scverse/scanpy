from itertools import product

from anndata import AnnData

import scanpy as sc
import scanpy.external as sce
from scanpy.testing._helpers.data import pbmc3k
from scanpy.testing._pytest.marks import needs


@needs("harmonypy")
def test_load_timepoints_from_anndata_list():
    adata_ref = pbmc3k()
    start = [596, 615, 1682, 1663, 1409, 1432]
    adata = AnnData.concatenate(
        *(adata_ref[i : i + 1000] for i in start),
        join="outer",
        batch_key="sample",
        batch_categories=[f"sa{i}_Rep{j}" for i, j in product((1, 2, 3), (1, 2))],
    )
    adata.obs["time_points"] = adata.obs["sample"].str.split("_", expand=True)[0]
    adata.obs["time_points"] = adata.obs["time_points"].astype("category")
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)

    sce.tl.harmony_timeseries(adata=adata, tp="time_points", n_components=None)
    assert all(
        [adata.obsp['harmony_aff'].shape[0], adata.obsp['harmony_aff_aug'].shape[0]]
    ), "harmony_timeseries augmented affinity matrix Error!"


@needs("harmonypy")
def test_harmony_integrate():
    """
    Test that Harmony integrate works.

    This is a very simple test that just checks to see if the Harmony
    integrate wrapper succesfully added a new field to ``adata.obsm``
    and makes sure it has the same dimensions as the original PCA table.
    """
    adata = pbmc3k()
    sc.pp.recipe_zheng17(adata)
    sc.tl.pca(adata)
    adata.obs['batch'] = 1350 * ['a'] + 1350 * ['b']
    sce.pp.harmony_integrate(adata, 'batch')
    assert adata.obsm['X_pca_harmony'].shape == adata.obsm['X_pca'].shape

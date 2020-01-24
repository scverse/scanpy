import scanpy as sc
import scanpy.external as sce

sample_names = ['sa1_Rep1', 'sa1_Rep2', 'sa2_Rep1', 'sa2_Rep2']
adata_ref = sc.datasets.pbmc3k()
n = int(adata_ref.shape[0] / 4)
datalist = [adata_ref[n * i : n * (i + 1)] for i in range(4)]
timepoints = [i.split("_")[0] for i in sample_names]
for ad, sn, tp in zip(datalist, sample_names, timepoints):
    ad.obs["time_points"] = tp
    ad.obs["sample_name"] = sn
adata = datalist[0].concatenate(*datalist[1:], join="outer")
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)


def load_timepoints_from_anndata_list():

    sce.tl.harmony_timeseries(adata=adata, tp="time_points", n_components=None)
    assert all(
        [adata.uns['aff_aug_harmony'].shape[0], adata.uns['aff_harmony'].shape[0]]
    ), "harmony_timeseries augmented affinity matrix Error!"

import scanpy as sc
import scanpy.external as sce

sample_names = ['sa1_Rep1', 'sa1_Rep2', 'sa2_Rep1', 'sa2_Rep2']
adata_ref = sc.datasets.pbmc3k()
n = int(adata_ref.shape[0] / 4)
datalist = [adata_ref[n * i : n * (i + 1)] for i in range(4)]


def load_timepoints_from_anndata_list():

    d = sce.tl.harmony_timeseries(data=datalist, sample_names=sample_names)
    d.log_transform = True
    d.process()
    assert isinstance(
        d.harmony_adata, sc.AnnData
    ), "harmony_timeseries augmented affinity matrix Error!"

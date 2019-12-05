import scanpy.external as sce
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'harmony_ext')
COUNTS = os.listdir(ROOT)
COUNTS = [os.path.join(ROOT, i) for i in COUNTS]
sample_names = ['sa1_Rep1', 'sa1_Rep2', 'sa2_Rep1', 'sa2_Rep2']


def load_timepoints():

    d = sce.tl.harmony(adata=COUNTS, sample_names=sample_names)
    d.log_transform = True
    d.process()

    assert d


def load_timepoints_from_anndata_list():

    adatalist = [sc.read(i) for i in COUNTS]
    d = sce.tl.harmony(adata=adatalist, sample_names=sample_names)
    d.log_transform = True
    d.process()

    assert d

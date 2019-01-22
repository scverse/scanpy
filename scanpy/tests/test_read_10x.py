import scanpy as sc
import numpy as np
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(ROOT, '_data', '10x_data')

def assert_anndata_equal(a1, a2):
    assert a1.shape == a2.shape
    assert all(a1.obs == a2.obs)
    assert all(a1.var == a2.var)
    assert np.allclose(a1.X.todense(), a2.X.todense())

def test_read_10x_mtx():
    sc.read_10x_mtx(os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices', 'hg19_chr21'),
                    var_names='gene_symbols', cache=True)
    sc.read_10x_mtx(os.path.join(ROOT, '3.0.0', 'filtered_feature_bc_matrix'),
                    var_names='gene_symbols', cache=True)

def test_read_10x_h5_v1():
    spec_genome_v1 = sc.read_10x_h5(os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices_h5.h5'),
                                    genome='hg19_chr21')
    nospec_genome_v1 = sc.read_10x_h5(os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices_h5.h5'))
    assert_anndata_equal(spec_genome_v1, nospec_genome_v1)

def test_read_10x_h5():
    spec_genome_v3 = sc.read_10x_h5(os.path.join(ROOT, '3.0.0', 'filtered_feature_bc_matrix.h5'), 
                                    genome='GRCh38_chr21')
    nospec_genome_v3 = sc.read_10x_h5(os.path.join(ROOT, '3.0.0', 'filtered_feature_bc_matrix.h5'))
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)

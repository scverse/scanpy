import scanpy.api as sc
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(ROOT, '_data', '10x_data')


def test_read_10x_mtx():
    sc.read_10x_mtx(os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices', 'hg19_chr21'),
                    var_names='gene_symbols', cache=True)
    sc.read_10x_mtx(os.path.join(ROOT, '3.0.0', 'filtered_feature_bc_matrix'),
                    var_names='gene_symbols', cache=True)


def test_read_10x_h5():
    sc.read_10x_h5(os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices_h5.h5'),
                   genome='hg19_chr21')
    sc.read_10x_h5(os.path.join(ROOT, '3.0.0',  'filtered_feature_bc_matrix.h5'),
                   genome='GRCh38_chr21')

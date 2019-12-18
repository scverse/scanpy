<<<<<<< HEAD
import h5py
import numpy as np
import os
=======
from pathlib import Path

import h5py
import numpy as np
>>>>>>> upstream/master
import pytest
import scanpy as sc


<<<<<<< HEAD
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

def test_error_10x_h5_legacy(tmp_path):
    onepth = os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices_h5.h5')
    twopth = str(tmp_path / "two_genomes.h5")
=======
ROOT = Path(__file__).parent
ROOT = ROOT / '_data' / '10x_data'


def assert_anndata_equal(a1, a2):
    assert a1.shape == a2.shape
    assert (a1.obs == a2.obs).all(axis=None)
    assert (a1.var == a2.var).all(axis=None)
    assert np.allclose(a1.X.todense(), a2.X.todense())


def test_read_10x_v1():
    v1_mtx = sc.read_10x_mtx(
        ROOT / '1.2.0' / 'filtered_gene_bc_matrices' / 'hg19_chr21',
        var_names='gene_symbols'
    )
    v1_h5 = sc.read_10x_h5(ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5')
    assert_anndata_equal(v1_mtx, v1_h5)


def test_read_10x_v3():
    v3_mtx = sc.read_10x_mtx(
        ROOT / '3.0.0' / 'filtered_feature_bc_matrix',
        var_names='gene_symbols',
    )
    v3_h5 = sc.read_10x_h5(ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5')
    v3_h5.var.drop(columns="genome", inplace=True)
    assert_anndata_equal(v3_mtx, v3_h5)


def test_read_10x_h5_v1():
    spec_genome_v1 = sc.read_10x_h5(
        ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5',
        genome='hg19_chr21',
    )
    nospec_genome_v1 = sc.read_10x_h5(ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5')
    assert_anndata_equal(spec_genome_v1, nospec_genome_v1)


def test_read_10x_h5():
    spec_genome_v3 = sc.read_10x_h5(
        ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5',
        genome='GRCh38_chr21',
    )
    nospec_genome_v3 = sc.read_10x_h5(ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5')
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)


def test_error_10x_h5_legacy(tmp_path):
    onepth = ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5'
    twopth = tmp_path / "two_genomes.h5"
>>>>>>> upstream/master
    with h5py.File(onepth, "r") as one, h5py.File(twopth, "w") as two:
        one.copy("hg19_chr21", two)
        one.copy("hg19_chr21", two, name="hg19_chr21_copy")
    with pytest.raises(ValueError):
        sc.read_10x_h5(twopth)
    sc.read_10x_h5(twopth, genome="hg19_chr21_copy")

<<<<<<< HEAD
def test_error_missing_genome():
    legacy_pth = os.path.join(ROOT, '1.2.0', 'filtered_gene_bc_matrices_h5.h5')
    v3_pth = os.path.join(ROOT, '3.0.0', 'filtered_feature_bc_matrix.h5')
=======

def test_error_missing_genome():
    legacy_pth = ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5'
    v3_pth = ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5'
>>>>>>> upstream/master
    with pytest.raises(ValueError, match=r".*hg19_chr21.*"):
        sc.read_10x_h5(legacy_pth, genome="not a genome")
    with pytest.raises(ValueError, match=r".*GRCh38_chr21.*"):
        sc.read_10x_h5(v3_pth, genome="not a genome")

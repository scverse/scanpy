from pathlib import Path

import h5py
import numpy as np
import pytest
import scanpy as sc


ROOT = Path(__file__).parent
ROOT = ROOT / '_data' / '10x_data'


def assert_anndata_equal(a1, a2):
    assert a1.shape == a2.shape
    assert (a1.obs == a2.obs).all(axis=None)
    assert (a1.var == a2.var).all(axis=None)
    assert np.allclose(a1.X.todense(), a2.X.todense())


@pytest.mark.parametrize(
    ['mtx_path', 'h5_path'],
    [
        pytest.param(
            ROOT / '1.2.0' / 'filtered_gene_bc_matrices' / 'hg19_chr21',
            ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5',
        ),
        pytest.param(
            ROOT / '3.0.0' / 'filtered_feature_bc_matrix',
            ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5',
        ),
    ],
)
def test_read_10x(tmp_path, mtx_path, h5_path):
    mtx = sc.read_10x_mtx(mtx_path, var_names="gene_symbols")
    h5 = sc.read_10x_h5(h5_path)

    # Drop genome column for comparing v3
    if "3.0.0" in str(h5_path):
        h5.var.drop(columns="genome", inplace=True)

    # Check equivalence
    assert_anndata_equal(mtx, h5)

    # Test that it can be written:
    from_mtx_pth = tmp_path / "from_mtx.h5ad"
    from_h5_pth = tmp_path / "from_h5.h5ad"

    mtx.write(from_mtx_pth)
    h5.write(from_h5_pth)

    assert_anndata_equal(sc.read_h5ad(from_mtx_pth), sc.read_h5ad(from_h5_pth))


def test_read_10x_h5_v1():
    spec_genome_v1 = sc.read_10x_h5(
        ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5', genome='hg19_chr21',
    )
    nospec_genome_v1 = sc.read_10x_h5(
        ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5'
    )
    assert_anndata_equal(spec_genome_v1, nospec_genome_v1)


def test_read_10x_h5():
    spec_genome_v3 = sc.read_10x_h5(
        ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5', genome='GRCh38_chr21',
    )
    nospec_genome_v3 = sc.read_10x_h5(ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5')
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)


def test_error_10x_h5_legacy(tmp_path):
    onepth = ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5'
    twopth = tmp_path / "two_genomes.h5"
    with h5py.File(onepth, "r") as one, h5py.File(twopth, "w") as two:
        one.copy("hg19_chr21", two)
        one.copy("hg19_chr21", two, name="hg19_chr21_copy")
    with pytest.raises(ValueError):
        sc.read_10x_h5(twopth)
    sc.read_10x_h5(twopth, genome="hg19_chr21_copy")


def test_error_missing_genome():
    legacy_pth = ROOT / '1.2.0' / 'filtered_gene_bc_matrices_h5.h5'
    v3_pth = ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5'
    with pytest.raises(ValueError, match=r".*hg19_chr21.*"):
        sc.read_10x_h5(legacy_pth, genome="not a genome")
    with pytest.raises(ValueError, match=r".*GRCh38_chr21.*"):
        sc.read_10x_h5(v3_pth, genome="not a genome")


def test_read_visium_counts():
    # TODO: What is the purpose of this test?
    h5_pth = ROOT / '../visium_data/1.0.0/filtered_feature_bc_matrix.h5'
    spec_genome_v3 = sc.read_10x_h5(h5_pth, genome='GRCh38')
    nospec_genome_v3 = sc.read_10x_h5(h5_pth)
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)


def test_10x_h5_gex():
    # Tests that gex option doesn't, say, make the function return None
    h5_pth = ROOT / '3.0.0' / 'filtered_feature_bc_matrix.h5'
    assert_anndata_equal(
        sc.read_10x_h5(h5_pth, gex_only=True), sc.read_10x_h5(h5_pth, gex_only=False)
    )

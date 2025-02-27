from __future__ import annotations

import shutil
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import pytest

import scanpy as sc

ROOT = Path(__file__).parent
ROOT = ROOT / "_data" / "10x_data"
VISIUM_ROOT = Path(__file__).parent / "_data" / "visium_data"


def assert_anndata_equal(a1, a2):
    assert a1.shape == a2.shape
    assert (a1.obs == a2.obs).all(axis=None)
    assert (a1.var == a2.var).all(axis=None)
    assert np.allclose(a1.X.todense(), a2.X.todense())


@pytest.mark.parametrize(
    ("mtx_path", "h5_path"),
    [
        pytest.param(
            ROOT / "1.2.0" / "filtered_gene_bc_matrices" / "hg19_chr21",
            ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5",
        ),
        pytest.param(
            ROOT / "3.0.0" / "filtered_feature_bc_matrix",
            ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5",
        ),
    ],
)
@pytest.mark.parametrize("prefix", [None, "prefix_"])
def test_read_10x(tmp_path, mtx_path, h5_path, prefix):
    if prefix is not None:
        # Build files named "prefix_XXX.xxx" in a temporary directory.
        mtx_path_orig = mtx_path
        mtx_path = tmp_path / "filtered_gene_bc_matrices_prefix"
        mtx_path.mkdir()
        for item in mtx_path_orig.iterdir():
            if item.is_file():
                shutil.copyfile(item, mtx_path / f"{prefix}{item.name}")

    mtx = sc.read_10x_mtx(mtx_path, var_names="gene_symbols", prefix=prefix)
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
        ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5",
        genome="hg19_chr21",
    )
    nospec_genome_v1 = sc.read_10x_h5(
        ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    )
    assert_anndata_equal(spec_genome_v1, nospec_genome_v1)


def test_read_10x_h5_v2_multiple_genomes():
    genome1_v1 = sc.read_10x_h5(
        ROOT / "1.2.0" / "multiple_genomes.h5",
        genome="hg19_chr21",
    )
    genome2_v1 = sc.read_10x_h5(
        ROOT / "1.2.0" / "multiple_genomes.h5",
        genome="another_genome",
    )
    # the test data are such that X is the same shape for both "genomes",
    # but the values are different
    assert (genome1_v1.X != genome2_v1.X).sum() > 0, (
        "loading data from two different genomes in 10x v2 format. "
        "should be different, but is the same. "
    )


def test_read_10x_h5():
    spec_genome_v3 = sc.read_10x_h5(
        ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5",
        genome="GRCh38_chr21",
    )
    nospec_genome_v3 = sc.read_10x_h5(ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5")
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)


def test_error_10x_h5_legacy(tmp_path):
    onepth = ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    twopth = tmp_path / "two_genomes.h5"
    with h5py.File(onepth, "r") as one, h5py.File(twopth, "w") as two:
        one.copy("hg19_chr21", two)
        one.copy("hg19_chr21", two, name="hg19_chr21_copy")
    with pytest.raises(ValueError, match=r"contains more than one genome"):
        sc.read_10x_h5(twopth)
    sc.read_10x_h5(twopth, genome="hg19_chr21_copy")


def test_error_missing_genome():
    legacy_pth = ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    v3_pth = ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5"
    with pytest.raises(ValueError, match=r".*hg19_chr21.*"):
        sc.read_10x_h5(legacy_pth, genome="not a genome")
    with pytest.raises(ValueError, match=r".*GRCh38_chr21.*"):
        sc.read_10x_h5(v3_pth, genome="not a genome")


@pytest.fixture(params=[1, 2])
def visium_pth(request, tmp_path) -> Path:
    visium1_pth = VISIUM_ROOT / "1.0.0"
    if request.param == 1:
        return visium1_pth
    elif request.param == 2:
        visium2_pth = tmp_path / "visium2"
        with patch.object(shutil, "copystat"):
            # copy only data, not file metadata
            shutil.copytree(visium1_pth, visium2_pth)
        header = "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"
        orig = visium2_pth / "spatial" / "tissue_positions_list.csv"
        csv = f"{header}\n{orig.read_text()}"
        orig.unlink()
        (orig.parent / "tissue_positions.csv").write_text(csv)
        return visium2_pth
    else:
        pytest.fail("add branch for new visium version")


@pytest.mark.filterwarnings("ignore:Use `squidpy.*` instead:FutureWarning")
def test_read_visium_counts(visium_pth):
    """Test checking that read_visium reads the right genome."""
    spec_genome_v3 = sc.read_visium(visium_pth, genome="GRCh38")
    nospec_genome_v3 = sc.read_visium(visium_pth)
    assert_anndata_equal(spec_genome_v3, nospec_genome_v3)


def test_10x_h5_gex():
    # Tests that gex option doesn't, say, make the function return None
    h5_pth = ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5"
    assert_anndata_equal(
        sc.read_10x_h5(h5_pth, gex_only=True), sc.read_10x_h5(h5_pth, gex_only=False)
    )


def test_10x_probe_barcode_read():
    # Tests the 10x probe barcode matrix is read correctly
    h5_pth = VISIUM_ROOT / "2.1.0" / "raw_probe_bc_matrix.h5"
    probe_anndata = sc.read_10x_h5(h5_pth)
    assert set(probe_anndata.var.columns) == {
        "feature_types",
        "filtered_probes",
        "gene_ids",
        "gene_name",
        "genome",
        "probe_ids",
        "probe_region",
    }
    assert set(probe_anndata.obs.columns) == {"filtered_barcodes"}
    assert probe_anndata.shape == (4987, 1000)
    assert probe_anndata.X.nnz == 858

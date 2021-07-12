import scanpy as sc
import numpy as np


def test_paga_positions_reproducible():
    """Check exact reproducibility and effect of random_state on paga positions"""
    # https://github.com/theislab/scanpy/issues/1859
    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.paga(pbmc, "bulk_labels")

    a = pbmc.copy()
    b = pbmc.copy()
    c = pbmc.copy()

    sc.pl.paga(a, show=False, random_state=42)
    sc.pl.paga(b, show=False, random_state=42)
    sc.pl.paga(c, show=False, random_state=13)

    np.testing.assert_array_equal(a.uns["paga"]["pos"], b.uns["paga"]["pos"])
    assert a.uns["paga"]["pos"].tolist() != c.uns["paga"]["pos"].tolist()

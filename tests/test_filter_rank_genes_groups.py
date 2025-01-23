from __future__ import annotations

import numpy as np
import pytest

from scanpy.tools import filter_rank_genes_groups, rank_genes_groups
from testing.scanpy._helpers.data import pbmc68k_reduced

NAMES_NO_REF = [
    ["CD3D", "ITM2A", "CD3D", "CCL5", "CD7", "nan", "CD79A", "nan", "NKG7", "LYZ"],
    ["CD3E", "CD3D", "nan", "NKG7", "CD3D", "AIF1", "CD79B", "nan", "GNLY", "CST3"],
    ["IL32", "RPL39", "nan", "CST7", "nan", "nan", "nan", "SNHG7", "CD7", "nan"],
    ["nan", "SRSF7", "IL32", "GZMA", "nan", "LST1", "IGJ", "nan", "CTSW", "nan"],
    ["nan", "nan", "CD2", "CTSW", "CD8B", "TYROBP", "ISG20", "SNHG8", "GZMB", "nan"],
]

NAMES_REF = [
    ["CD3D", "ITM2A", "CD3D", "nan", "CD3D", "nan", "CD79A", "nan", "CD7"],
    ["nan", "nan", "nan", "CD3D", "nan", "AIF1", "nan", "nan", "NKG7"],
    ["nan", "nan", "nan", "NKG7", "nan", "FCGR3A", "ISG20", "SNHG7", "CTSW"],
    ["nan", "CD3D", "nan", "CCL5", "CD7", "nan", "CD79B", "nan", "GNLY"],
    ["CD3E", "IL32", "nan", "IL32", "CD27", "FCER1G", "nan", "nan", "nan"],
]

NAMES_NO_REF_COMPARE_ABS = [
    [
        *("CD3D", "ITM2A", "HLA-DRB1", "CCL5", "HLA-DPA1"),
        *("nan", "CD79A", "nan", "NKG7", "LYZ"),
    ],
    [
        *("HLA-DPA1", "nan", "CD3D", "NKG7", "HLA-DRB1"),
        *("AIF1", "CD79B", "nan", "GNLY", "CST3"),
    ],
    [
        *("nan", "PSAP", "CD74", "CST7", "CD74"),
        *("PSAP", "FCER1G", "SNHG7", "CD7", "HLA-DRA"),
    ],
    [
        *("IL32", "nan", "HLA-DRB5", "GZMA", "HLA-DRB5"),
        *("LST1", "nan", "nan", "CTSW", "HLA-DRB1"),
    ],
    [
        *("nan", "FCER1G", "HLA-DPB1", "CTSW", "HLA-DPB1"),
        *("TYROBP", "TYROBP", "S100A10", "GZMB", "HLA-DPA1"),
    ],
]


EXPECTED = {
    ("Dendritic", False): np.array(NAMES_REF),
    ("rest", False): np.array(NAMES_NO_REF),
    ("rest", True): np.array(NAMES_NO_REF_COMPARE_ABS),
}


@pytest.mark.parametrize(
    ("reference", "pts", "abs"),
    [
        pytest.param("Dendritic", False, False, id="ref-no_pts-no_abs"),
        pytest.param("rest", False, False, id="rest-no_pts-no_abs"),
        pytest.param("rest", True, False, id="rest-pts-no_abs"),
        pytest.param("rest", True, True, id="rest-pts-abs"),
    ],
)
def test_filter_rank_genes_groups(reference, pts, abs):
    adata = pbmc68k_reduced()

    rank_genes_groups(
        adata,
        "bulk_labels",
        reference=reference,
        pts=pts,
        method="wilcoxon",
        rankby_abs=abs,
        n_genes=5,
    )
    if abs:
        filter_rank_genes_groups(
            adata,
            compare_abs=True,
            min_in_group_fraction=-1,
            max_out_group_fraction=1,
            min_fold_change=3.1,
        )
    else:
        filter_rank_genes_groups(
            adata,
            min_in_group_fraction=0.25,
            min_fold_change=1,
            max_out_group_fraction=0.5,
        )

    assert np.array_equal(
        EXPECTED[reference, abs],
        np.array(adata.uns["rank_genes_groups_filtered"]["names"].tolist()),
    )

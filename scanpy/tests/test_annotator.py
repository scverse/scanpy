from importlib.util import find_spec

import pandas as pd
import numpy as np
from scanpy import AnnData
from scanpy.external.tl import annotator

import pytest


@pytest.fixture
def markers():
    return pd.DataFrame(
        [
            ["Type 1", "111"],
            ["Type 1", "112"],
            ["Type 1", "113"],
            ["Type 1", "114"],
            ["Type 2", "211"],
            ["Type 2", "212"],
            ["Type 2", "213"],
            ["Type 2", "214"],
        ],
        columns=["Cell Type", "Gene"],
    )


@pytest.fixture
def data():
    genes = ["111", "112", "113", "114", "211", "212", "213", "214"]
    return pd.DataFrame(
        np.array(
            [
                [1, 1, 1, 1.1, 0, 0, 0, 0],
                [1, 0.8, 0.9, 1, 0, 0, 0, 0],
                [0.7, 1.1, 1, 1.2, 0, 0, 0, 0],
                [0.8, 0.7, 1.1, 1, 0, 0.1, 0, 0],
                [0, 0, 0, 0, 1.05, 1.05, 1.1, 1],
                [0, 0, 0, 0, 1.1, 1.0, 1.05, 1.1],
                [0, 0, 0, 0, 1.05, 0.9, 1.1, 1.1],
                [0, 0, 0, 0, 0.9, 0.9, 1.2, 1],
            ]
        ),
        columns=genes,
    )


@pytest.fixture
def adata(data):
    return AnnData(data.values, var=data.columns.values)


def basic_check(annotations, adata):
    assert type(annotations) == AnnData
    assert len(annotations) == len(adata)
    assert annotations.shape == (8, 2)  # two types in the data
    assert np.nansum(annotations.X) > 0
    assert np.nanmax(annotations.X) <= 1
    assert np.nanmin(annotations.X) >= 0


@pytest.mark.skipif(
    find_spec('pointannotator') is None, reason="point-annotator not installed"
)
def test_annotator(adata, markers):
    annotations = annotator(
        adata, markers, normalize=False, num_genes=15
    )

    basic_check(annotations, adata)


@pytest.mark.skipif(
    find_spec('pointannotator') is None, reason="point-annotator not installed"
)
def test_remove_empty_column(adata, markers):
    """
    Type 3 column must be removed here, since this cell type does not
    belong to any cell.
    """
    additinal_markers = pd.DataFrame(
        [["Type 3", "311"], ["Type 3", "312"], ["Type 3", "313"]],
        columns=["Cell Type", "Gene"],
    )
    markers = markers.append(additinal_markers)

    annotations = annotator(adata, markers, num_genes=20)

    basic_check(annotations, adata)

    annotations = annotator(
        adata,
        markers,
        num_genes=20,
        return_nonzero_annotations=False,
    )
    assert len(annotations) == len(adata)
    assert annotations.shape == (8, 3)  # three types in the data
    assert np.nansum(annotations.X) > 0
    assert np.nanmax(annotations.X) <= 1
    assert np.nanmin(annotations.X) >= 0


@pytest.mark.skipif(
    find_spec('pointannotator') is None, reason="point-annotator not installed"
)
def test_sf(adata, markers):
    """
    Test annotations with hypergeom.sf
    """
    annotations = annotator(
        adata, markers, num_genes=15, p_value_fun="hypergeom"
    )

    basic_check(annotations, adata)


@pytest.mark.skipif(
    find_spec('pointannotator') is None, reason="point-annotator not installed"
)
def test_scoring(adata, markers, data):
    # scoring SCORING_EXP_RATIO
    annotations = annotator(
        adata, markers, num_genes=15, scoring="exp_ratio"
    )

    basic_check(annotations, adata)

    # scoring SCORING_MARKERS_SUM
    annotations = annotator(
        adata,
        markers,
        num_genes=15,
        scoring="sum_of_expressed_markers",
    )

    assert type(annotations) == AnnData
    assert len(annotations) == len(adata)
    assert annotations.shape == (8, 2)  # two types in the data

    # based on provided data it should match
    # the third row is skipped, since it is special
    assert pytest.approx(annotations.X[0, 0]) == data.iloc[0].sum()
    assert pytest.approx(annotations.X[5, 1]) == data.iloc[5].sum()

    # scoring SCORING_LOG_FDR
    annotations = annotator(
        adata, markers, num_genes=15, scoring="log_fdr"
    )

    assert type(annotations) == AnnData
    assert len(annotations) == len(adata)
    assert annotations.shape == (8, 2)  # two types in the data

    # scoring SCORING_LOG_PVALUE
    annotations = annotator(
        adata, markers, num_genes=15, scoring="log_p_value"
    )

    assert type(annotations) == AnnData
    assert len(annotations) == len(adata)
    assert annotations.shape == (8, 2)  # two types in the data

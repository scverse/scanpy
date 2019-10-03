import unittest
from importlib.util import find_spec

import pandas as pd
import numpy as np
from scanpy import AnnData
from scanpy.external.tl import annotator

import pytest


@pytest.mark.skipif(
    find_spec('pointannotator') is None, reason="point-annotator not installed"
)
class AnnotatorTests(unittest.TestCase):
    def setUp(self):
        self.markers = pd.DataFrame(
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

        genes = ["111", "112", "113", "114", "211", "212", "213", "214"]
        self.data = pd.DataFrame(
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

        # transform data to AnnData
        self.anndata = AnnData(self.data.values, var=self.data.columns.values)

    def basic_check(self, annotations):
        self.assertEqual(type(annotations), AnnData)
        self.assertEqual(len(annotations), len(self.anndata))
        self.assertTupleEqual(
            annotations.shape, (8, 2)
        )  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_annotator(self):
        annotations = annotator(
            self.anndata, self.markers, normalize=False, num_genes=15
        )

        self.basic_check(annotations)

    def test_remove_empty_column(self):
        """
        Type 3 column must be removed here, since this cell type does not
        belong to any cell.
        """
        additinal_markers = pd.DataFrame(
            [["Type 3", "311"], ["Type 3", "312"], ["Type 3", "313"]],
            columns=["Cell Type", "Gene"],
        )
        markers = self.markers.append(additinal_markers)

        annotations = annotator(self.anndata, markers, num_genes=20)

        self.basic_check(annotations)

        annotations = annotator(
            self.anndata,
            markers,
            num_genes=20,
            return_nonzero_annotations=False,
        )
        self.assertEqual(len(annotations), len(self.anndata))
        self.assertTupleEqual(
            annotations.shape, (8, 3)
        )  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_sf(self):
        """
        Test annotations with hypergeom.sf
        """
        annotations = annotator(
            self.anndata, self.markers, num_genes=15, p_value_fun="hypergeom"
        )

        self.basic_check(annotations)

    def test_scoring(self):
        # scoring SCORING_EXP_RATIO
        annotations = annotator(
            self.anndata, self.markers, num_genes=15, scoring="exp_ratio"
        )

        self.basic_check(annotations)

        # scoring SCORING_MARKERS_SUM
        annotations = annotator(
            self.anndata,
            self.markers,
            num_genes=15,
            scoring="sum_of_expressed_markers",
        )

        self.assertEqual(type(annotations), AnnData)
        self.assertEqual(len(annotations), len(self.anndata))
        self.assertTupleEqual(
            annotations.shape, (8, 2)
        )  # two types in the data

        # based on provided data it should match
        # the third row is skipped, since it is special
        self.assertAlmostEqual(
            annotations.X[0, 0], self.data.iloc[0].sum(), places=6
        )
        self.assertAlmostEqual(
            annotations.X[5, 1], self.data.iloc[5].sum(), places=6
        )

        # scoring SCORING_LOG_FDR
        annotations = annotator(
            self.anndata, self.markers, num_genes=15, scoring="log_fdr"
        )

        self.assertEqual(type(annotations), AnnData)
        self.assertEqual(len(annotations), len(self.anndata))
        self.assertTupleEqual(
            annotations.shape, (8, 2)
        )  # two types in the data

        # scoring SCORING_LOG_PVALUE
        annotations = annotator(
            self.anndata, self.markers, num_genes=15, scoring="log_p_value"
        )

        self.assertEqual(type(annotations), AnnData)
        self.assertEqual(len(annotations), len(self.anndata))
        self.assertTupleEqual(
            annotations.shape, (8, 2)
        )  # two types in the data

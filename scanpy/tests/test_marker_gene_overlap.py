from __future__ import annotations

import numpy as np
from anndata import AnnData

import scanpy as sc


def generate_test_data():
    # Create an artificial data set
    test_data = AnnData(X=np.ones((9, 10)))
    test_data.uns["rank_genes_groups"] = dict()
    test_data.uns["rank_genes_groups"]["names"] = np.rec.fromarrays(
        [["a", "b", "c", "d", "e"], ["a", "f", "g", "h", "i"]], names="c0,c1"
    )
    test_data.uns["rank_genes_groups"]["pvals_adj"] = np.rec.fromarrays(
        [[0.001, 0.01, 0.02, 0.05, 0.6], [0.001, 0.01, 0.02, 0.05, 0.6]], names="c0,c1"
    )

    marker_genes = {"type 1": {"a", "b", "c"}, "type 2": {"a", "f", "g"}}

    return test_data, marker_genes


def test_marker_overlap_base():
    # Test all overlap calculations on artificial data
    test_data, marker_genes = generate_test_data()

    t1 = sc.tl.marker_gene_overlap(test_data, marker_genes)

    assert t1["c0"]["type 1"] == 3.0
    assert t1["c1"]["type 2"] == 3.0


def test_marker_overlap_normalization():
    test_data, marker_genes = generate_test_data()

    t2 = sc.tl.marker_gene_overlap(test_data, marker_genes, normalize="reference")
    t3 = sc.tl.marker_gene_overlap(test_data, marker_genes, normalize="data")

    assert t2["c0"]["type 1"] == 1.0
    assert t3["c1"]["type 2"] == 0.6


def test_marker_overlap_methods():
    test_data, marker_genes = generate_test_data()

    t4 = sc.tl.marker_gene_overlap(test_data, marker_genes, method="overlap_coef")
    t5 = sc.tl.marker_gene_overlap(test_data, marker_genes, method="jaccard")

    assert t4["c0"]["type 1"] == 1.0
    assert t5["c0"]["type 1"] == 0.6


def test_marker_overlap_subsetting():
    test_data, marker_genes = generate_test_data()

    t6 = sc.tl.marker_gene_overlap(test_data, marker_genes, top_n_markers=2)
    t7 = sc.tl.marker_gene_overlap(test_data, marker_genes, adj_pval_threshold=0.01)

    assert t6["c0"]["type 1"] == 2.0
    assert t7["c0"]["type 1"] == 1.0

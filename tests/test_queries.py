from __future__ import annotations

import pandas as pd
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


@pytest.mark.internet
@needs.gprofiler
def test_enrich():
    pbmc = pbmc68k_reduced()
    sc.tl.rank_genes_groups(pbmc, "louvain", n_genes=pbmc.shape[1])
    enrich_anndata = sc.queries.enrich(pbmc, "1")
    de = pd.DataFrame()
    for k in ["pvals_adj", "names"]:
        de[k] = pbmc.uns["rank_genes_groups"][k]["1"]
    de_genes = de.loc[lambda x: x["pvals_adj"] < 0.05, "names"]
    enrich_list = sc.queries.enrich(list(de_genes))
    assert (enrich_anndata == enrich_list).all().all()

    # scverse/scanpy/#1043
    sc.tl.filter_rank_genes_groups(pbmc, min_fold_change=1)
    sc.queries.enrich(pbmc, "1")

    gene_dict = {"set1": ["KLF4", "PAX5"], "set2": ["SOX2", "NANOG"]}
    enrich_list = sc.queries.enrich(
        gene_dict, org="hsapiens", gprofiler_kwargs=dict(sources=["GO:BP"])
    )
    assert "set1" in enrich_list["query"].unique()
    assert "set2" in enrich_list["query"].unique()


@pytest.mark.internet
@needs.pybiomart
def test_mito_genes():
    pbmc = pbmc68k_reduced()
    mt_genes = sc.queries.mitochondrial_genes("hsapiens")
    assert (
        pbmc.var_names.isin(mt_genes["external_gene_name"]).sum() == 1
    )  # Should only be MT-ND3

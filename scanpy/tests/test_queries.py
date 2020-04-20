import pandas as pd
import pytest
import scanpy as sc


@pytest.mark.internet
def test_enrich():
    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.rank_genes_groups(pbmc, "louvain", n_genes=pbmc.shape[1])
    enrich_anndata = sc.queries.enrich(pbmc, "1")
    de = pd.DataFrame()
    for k in ["pvals_adj", "names"]:
        de[k] = pbmc.uns["rank_genes_groups"][k]["1"]
    de_genes = de.loc[lambda x: x["pvals_adj"] < 0.05, "names"]
    enrich_list = sc.queries.enrich(list(de_genes))
    assert (enrich_anndata == enrich_list).all().all()

    # theislab/scanpy/#1043
    sc.tl.filter_rank_genes_groups(pbmc, min_fold_change=1)
    sc.queries.enrich(pbmc, "1")


@pytest.mark.internet
def test_mito_genes():
    pbmc = sc.datasets.pbmc68k_reduced()
    mt_genes = sc.queries.mitochondrial_genes("hsapiens")
    assert (
        pbmc.var_names.isin(mt_genes["external_gene_name"]).sum() == 1
    )  # Should only be MT-ND3

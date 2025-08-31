from __future__ import annotations

from functools import partial

import pytest
from anndata import read_h5ad

import scanpy as sc


@pytest.mark.flaky(reruns=5, reruns_delay=2)
@pytest.mark.parametrize(
    ("name", "func", "msg"),
    [
        pytest.param("PCA", sc.pp.pca, " with chunked as False", id="pca"),
        pytest.param(
            "PCA", partial(sc.pp.pca, layer="X_copy"), " from layers", id="pca_layer"
        ),
        pytest.param(
            "regress_out",
            partial(sc.pp.regress_out, keys=["n_counts", "percent_mito"]),
            "",
            id="regress_out",
        ),
        pytest.param(
            "dendrogram", partial(sc.tl.dendrogram, groupby="cat"), "", id="dendrogram"
        ),
        pytest.param("tsne", sc.tl.tsne, "", id="tsne"),
        pytest.param("scale", sc.pp.scale, "", id="scale"),
        pytest.param(
            "downsample_counts",
            partial(sc.pp.downsample_counts, counts_per_cell=1000),
            "",
            id="downsample_counts",
        ),
        pytest.param(
            "filter_genes",
            partial(sc.pp.filter_genes, max_cells=1000),
            "",
            id="filter_genes",
        ),
        pytest.param(
            "filter_cells",
            partial(sc.pp.filter_cells, max_genes=1000),
            "",
            id="filter_cells",
        ),
        pytest.param(
            "rank_genes_groups",
            partial(sc.tl.rank_genes_groups, groupby="cat"),
            "",
            id="rank_genes_groups",
        ),
        pytest.param(
            "score_genes",
            partial(sc.tl.score_genes, gene_list=map(str, range(100))),
            "",
            id="score_genes",
        ),
    ],
)
def test_backed_error(backed_adata, name, func, msg):
    with pytest.raises(
        NotImplementedError,
        match=f"{name} is not implemented for matrices of type {type(backed_adata.X)}{msg}",
    ):
        func(backed_adata)


def test_log1p_backed_errors(backed_adata):
    with pytest.raises(
        NotImplementedError,
        match="log1p is not implemented for backed AnnData with backed mode not r+",
    ):
        sc.pp.log1p(backed_adata, chunked=True)
    backed_adata.file.close()
    backed_adata = read_h5ad(backed_adata.filename, backed="r+")
    with pytest.raises(
        NotImplementedError,
        match=f"log1p is not implemented for matrices of type {type(backed_adata.X)} without `chunked=True`",
    ):
        sc.pp.log1p(backed_adata)
    backed_adata.layers["X_copy"] = backed_adata.X
    layer_type = type(backed_adata.layers["X_copy"])
    with pytest.raises(
        NotImplementedError,
        match=f"log1p is not implemented for matrices of type {layer_type} from layers",
    ):
        sc.pp.log1p(backed_adata, layer="X_copy")
    backed_adata.file.close()


def test_scatter_backed(backed_adata):
    sc.pp.pca(backed_adata, chunked=True)
    sc.pl.scatter(backed_adata, color="0", basis="pca", show=False)


def test_dotplot_backed(backed_adata):
    sc.pl.dotplot(backed_adata, ["0", "1", "2", "3"], groupby="cat", show=False)

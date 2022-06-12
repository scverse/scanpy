import scanpy as sc
from sklearn import preprocessing
from scanpy.tests._data._cached_datasets import pbmc68k_reduced
import pytest


def test_deprecate_multicore_tsne():
    pbmc = pbmc68k_reduced()

    with pytest.warns(
        UserWarning, match="calling tsne with n_jobs > 1 would use MulticoreTSNE"
    ):
        sc.tl.tsne(pbmc, n_jobs=2)

    with pytest.warns(FutureWarning, match="Argument `use_fast_tsne` is deprecated"):
        sc.tl.tsne(pbmc, use_fast_tsne=True)

    with pytest.warns(UserWarning, match="Falling back to scikit-learn"):
        sc.tl.tsne(pbmc, use_fast_tsne=True)


def test_deprecate_use_highly_variable_genes():    
    pbmc = pbmc68k_reduced()

    with pytest.warns(UserWarning, match="Argument `use_highly_variable` is deprecated"):
        sc.pp.pca(pbmc, use_highly_variable=True)

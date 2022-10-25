from functools import partial
from importlib.util import find_spec
from pathlib import Path

from matplotlib import cm
import numpy as np

import scanpy as sc
from scanpy.tests._data._cached_datasets import pbmc68k_reduced, pbmc3k_processed
import pytest

HERE: Path = Path(__file__).parent
ROOT = HERE / '_images'
FIGS = HERE / 'figures'


needs_igraph = pytest.mark.skipif(
    not find_spec("igraph"),
    reason="needs module `igraph` (`pip install python-igraph`)",
)


@pytest.fixture(scope="module")
def pbmc():
    pbmc = pbmc68k_reduced()
    sc.tl.paga(pbmc, groups='bulk_labels')
    pbmc.obs['cool_feature'] = pbmc[:, 'CST3'].X.squeeze()
    return pbmc


@needs_igraph
@pytest.mark.parametrize(
    "test_id,func",
    [
        ("master_paga", sc.pl.paga),
        ("master_paga_continuous", partial(sc.pl.paga, color="CST3")),
        ("master_paga_continuous_obs", partial(sc.pl.paga, color="cool_feature")),
        (
            "master_paga_continuous_multiple",
            partial(sc.pl.paga, color=['CST3', 'GATA2']),
        ),
        ("master_paga_compare", partial(sc.pl.paga_compare, legend_fontoutline=2)),
        (
            "master_paga_compare_continuous",
            partial(sc.pl.paga_compare, color='CST3', legend_fontsize=5),
        ),
        (
            "master_paga_compare_pca",
            partial(sc.pl.paga_compare, basis='X_pca', legend_fontweight='normal'),
        ),
    ],
)
def test_paga_plots(image_comparer, pbmc, test_id, func):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=30)
    common = dict(threshold=0.5, max_edge_width=1.0, random_state=0, show=False)

    func(pbmc, **common)
    save_and_compare_images(test_id)


@needs_igraph
def test_paga_pie(image_comparer, pbmc):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=30)

    colors = {
        c: {cm.Set1(_): 0.33 for _ in range(3)}
        for c in pbmc.obs["bulk_labels"].cat.categories
    }
    colors["Dendritic"] = {cm.Set2(_): 0.25 for _ in range(4)}

    sc.pl.paga(pbmc, color=colors, colorbar=False)
    save_and_compare_images('master_paga_pie')


@needs_igraph
def test_paga_path(image_comparer, pbmc):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc.uns['iroot'] = 0
    sc.tl.dpt(pbmc)
    sc.pl.paga_path(
        pbmc,
        nodes=['Dendritic'],
        keys=['HES4', 'SRM', 'CSTB'],
        show=False,
    )
    save_and_compare_images('master_paga_path')


@needs_igraph
def test_paga_compare(image_comparer):
    # Tests that https://github.com/scverse/scanpy/issues/1887 is fixed
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = pbmc3k_processed()
    sc.tl.paga(pbmc, groups="louvain")

    sc.pl.paga_compare(pbmc, basis="umap", show=False)

    save_and_compare_images('master_paga_compare_pbmc3k')


@needs_igraph
def test_paga_positions_reproducible():
    """Check exact reproducibility and effect of random_state on paga positions"""
    # https://github.com/scverse/scanpy/issues/1859
    pbmc = pbmc68k_reduced()
    sc.tl.paga(pbmc, "bulk_labels")

    a = pbmc.copy()
    b = pbmc.copy()
    c = pbmc.copy()

    sc.pl.paga(a, show=False, random_state=42)
    sc.pl.paga(b, show=False, random_state=42)
    sc.pl.paga(c, show=False, random_state=13)

    np.testing.assert_array_equal(a.uns["paga"]["pos"], b.uns["paga"]["pos"])
    assert a.uns["paga"]["pos"].tolist() != c.uns["paga"]["pos"].tolist()

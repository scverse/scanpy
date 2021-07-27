from functools import partial
from pathlib import Path

from matplotlib import cm

import scanpy as sc

import pytest

HERE: Path = Path(__file__).parent
ROOT = HERE / '_images'
FIGS = HERE / 'figures'


@pytest.fixture(scope="module")
def pbmc():
    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.paga(pbmc, groups='bulk_labels')
    pbmc.obs['cool_feature'] = pbmc[:, 'CST3'].X.squeeze()
    return pbmc


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


def test_paga_pie(image_comparer, pbmc):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=30)

    colors = {
        c: {cm.Set1(_): 0.33 for _ in range(3)}
        for c in pbmc.obs["bulk_labels"].cat.categories
    }
    colors["Dendritic"] = {cm.Set2(_): 0.25 for _ in range(4)}

    sc.pl.paga(pbmc, color=colors, colorbar=False)
    save_and_compare_images('master_paga_pie')


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


def test_paga_compare(image_comparer):
    # Tests that https://github.com/theislab/scanpy/issues/1887 is fixed
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc3k_processed()
    sc.tl.paga(pbmc, groups="louvain")

    sc.pl.paga_compare(pbmc, basis="umap", show=False)

    save_and_compare_images('master_paga_compare_pbmc3k')

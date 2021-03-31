import pytest

import scanpy as sc
import scanpy.external as sce
import numpy as np
import os
import os.path as osp


def test_scnym():
    """Test that scNym runs and returns plausible classifications
    on a demonstration dataset.
    """
    pytest.importorskip("scnym")
    TEST_URL = "https://storage.googleapis.com/calico-website-mca-storage/kang_2017_stim_pbmc.h5ad"

    np.random.seed(1)
    # load a testing dataset
    adata = sc.datasets._datasets.read(
        sc.settings.datasetdir / 'kang17.h5ad', backup_url=TEST_URL
    )
    target_bidx = adata.obs['stim'] == 'stim'
    adata.obs['cell'] = np.array(adata.obs['cell'])
    adata.obs['ground_truth'] = np.array(adata.obs['cell'])
    adata.obs.loc[target_bidx, 'cell'] = 'Unlabeled'

    # downsample to speed up testing
    ridx = np.random.choice(
        adata.shape[0],
        size=2048,
        replace=False,
    )
    adata = adata[ridx, :].copy()

    # train an scNym model
    print("training...")
    config = {'n_epochs': 1}
    sce.tl.scnym(
        adata=adata,
        task='train',
        groupby='cell',
        out_path=str(sc.settings.datasetdir),
        config=config,
    )

    assert 'scNym_train_results' in adata.uns.keys()
    assert osp.exists(
        osp.join(str(sc.settings.datasetdir), '00_best_model_weights.pkl')
    )

    # predict cell types with an scNym model
    print("predicting...")
    sce.tl.scnym(
        adata=adata,
        task='predict',
        key_added='scNym',
        out_path=str(sc.settings.datasetdir),
        trained_model=str(sc.settings.datasetdir),
        config=config,
    )

    assert 'X_scnym' in adata.obsm.keys()
    assert 'scNym' in adata.obs.columns

    # check accuracy
    target_gt = np.array(adata.obs.loc[target_bidx, 'ground_truth'])
    pred = np.array(adata.obs.loc[target_bidx, 'scNym'])
    acc = np.sum(target_gt == pred) / len(pred)
    print(f"acc: {acc:.6f}")
    assert acc > 0.5, "low accuracy suggests a performance issue"

    return


if __name__ == "__main__":
    test_scnym()

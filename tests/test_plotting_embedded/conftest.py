from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scanpy as sc

HERE: Path = Path(__file__).parent


@pytest.fixture(scope="module")
def adata():
    # A bit cute.
    from matplotlib.image import imread
    from sklearn.cluster import DBSCAN
    from sklearn.datasets import make_blobs

    empty_pixel = np.array([1.0, 1.0, 1.0, 0]).reshape(1, 1, -1)
    image = imread(HERE.parent.parent / "docs/_static/img/Scanpy_Logo_RGB.png")
    x, y = np.where(np.logical_and.reduce(~np.equal(image, empty_pixel), axis=2))

    # Just using to calculate the hex coords
    hexes = plt.hexbin(x, y, gridsize=(44, 100))
    counts = hexes.get_array()
    pixels = hexes.get_offsets()[counts != 0]
    plt.close()

    labels = DBSCAN(eps=20, min_samples=2).fit(pixels).labels_
    order = np.argsort(labels)
    adata = sc.AnnData(
        make_blobs(
            pd.Series(labels[order]).value_counts().values,
            n_features=20,
            shuffle=False,
            random_state=42,
        )[0],
        obs={"label": pd.Categorical(labels[order].astype(str))},
        obsm={"spatial": pixels[order, ::-1]},
        uns={
            "spatial": {
                "scanpy_img": {
                    "images": {"hires": image},
                    "scalefactors": {
                        "tissue_hires_scalef": 1,
                        "spot_diameter_fullres": 10,
                    },
                }
            }
        },
    )
    sc.pp.pca(adata)

    # Adding some missing values
    adata.obs["label_missing"] = adata.obs["label"].copy()
    adata.obs.loc[::2, "label_missing"] = np.nan

    adata.obs["1_missing"] = adata.obs_vector("1")
    adata.obs.loc[
        adata.obsm["spatial"][:, 0] < adata.obsm["spatial"][:, 0].mean(), "1_missing"
    ] = np.nan

    return adata
